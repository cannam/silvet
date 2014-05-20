/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  Silvet

  A Vamp plugin for note transcription.
  Centre for Digital Music, Queen Mary University of London.
    
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of the
  License, or (at your option) any later version.  See the file
  COPYING included with this distribution for more information.
*/

#include "EM.h"

#include <cstdlib>
#include <cmath>

#include <iostream>

#include "VectorOps.h"
#include "Allocators.h"
#include "Instruments.h"

using std::vector;
using std::cerr;
using std::endl;

using namespace breakfastquay;

static float epsilon = 1e-10;

EM::EM(const InstrumentPack *pack, bool useShifts) :
    m_pack(pack),
    m_noteCount(pack->templateNoteCount),
    m_shiftCount(useShifts ? pack->templateMaxShift * 2 + 1 : 1),
    m_binCount(pack->templateHeight),
    m_sourceCount(pack->templates.size()),
    m_pitchSparsity(1.1),
    //!!! note: slightly less source sparsity might help; also
    //!!! consider a modest shift sparsity e.g. 1.1
    m_sourceSparsity(1.3),
    m_lowestPitch(pack->lowestNote),
    m_highestPitch(pack->highestNote)
{
    m_pitches = allocate<float>(m_noteCount);
    m_updatePitches = allocate<float>(m_noteCount);
    for (int n = 0; n < m_noteCount; ++n) {
        m_pitches[n] = drand48();
    }

    if (useShifts) {
        m_shifts = allocate_channels<float>(m_shiftCount, m_noteCount);
        m_updateShifts = allocate_channels<float>(m_shiftCount, m_noteCount);
        for (int f = 0; f < m_shiftCount; ++f) {
            for (int n = 0; n < m_noteCount; ++n) {
                m_shifts[f][n] = drand48();
            }
        }
    } else {
        m_shifts = 0;
        m_updateShifts = 0;
    }
    
    m_sources = allocate_channels<float>(m_sourceCount, m_noteCount);
    m_updateSources = allocate_channels<float>(m_sourceCount, m_noteCount);
    for (int i = 0; i < m_sourceCount; ++i) {
        for (int n = 0; n < m_noteCount; ++n) {
            m_sources[i][n] = (inRange(i, n) ? 1.0 : 0.0);
        }
    }

    m_estimate = allocate<float>(m_binCount);
    m_q = allocate<float>(m_binCount);
}

EM::~EM()
{
    deallocate(m_q);
    deallocate(m_estimate);
    deallocate_channels(m_sources, m_sourceCount);
    deallocate_channels(m_updateSources, m_sourceCount);
    deallocate_channels(m_shifts, m_shiftCount);
    deallocate_channels(m_updateShifts, m_shiftCount);
    deallocate(m_pitches);
    deallocate(m_updatePitches);
}

void
EM::rangeFor(int instrument, int &minPitch, int &maxPitch)
{
    minPitch = m_pack->templates[instrument].lowestNote;
    maxPitch = m_pack->templates[instrument].highestNote;
}

bool
EM::inRange(int instrument, int pitch)
{
    int minPitch, maxPitch;
    rangeFor(instrument, minPitch, maxPitch);
    return (pitch >= minPitch && pitch <= maxPitch);
}

void
EM::normaliseColumn(float *column, int size)
{
    float sum = v_sum(column, size);
    v_scale(column, 1.0 / sum, size);
}

void
EM::normaliseGrid(float **grid, int size1, int size2)
{
    float *denominators = allocate_and_zero<float>(size2);

    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            denominators[j] += grid[i][j];
        }
    }

    for (int i = 0; i < size1; ++i) {
        v_divide(grid[i], denominators, size2);
    }

    deallocate(denominators);
}

void
EM::iterate(const double *column)
{
    float *norm = allocate<float>(m_binCount);
    v_convert(norm, column, m_binCount);
    normaliseColumn(norm, m_binCount);
    expectation(norm);
    maximisation(norm);
    deallocate(norm);
}

const float *
EM::templateFor(int instrument, int note, int shift)
{
    const float *base = m_pack->templates.at(instrument).data.at(note).data();
    if (m_shifts) {
        return base + shift;
    } else {
        return base + m_pack->templateMaxShift;
    }
}

void
EM::expectation(const float *column)
{
//    cerr << ".";

    v_set(m_estimate, epsilon, m_binCount);

    for (int f = 0; f < m_shiftCount; ++f) {

        const float *shiftIn = m_shifts ? m_shifts[f] : 0;

        for (int i = 0; i < m_sourceCount; ++i) {

            const float *sourceIn = m_sources[i];

            int lowest, highest;
            rangeFor(i, lowest, highest);

            for (int n = lowest; n <= highest; ++n) {

                const float source = sourceIn[n];
                const float shift = shiftIn ? shiftIn[n] : 1.0;
                const float pitch = m_pitches[n];

                const float factor = pitch * source * shift;
                const float *w = templateFor(i, n, f);

                v_add_with_gain(m_estimate, w, factor, m_binCount);
            }
        }
    }

    //!!! try normalising so as to sum to unity and then taking L2
    //!!! norm of the two (for measure of how close they are)

    for (int i = 0; i < m_binCount; ++i) {
        m_q[i] = column[i] / m_estimate[i];
    }
}

void
EM::maximisation(const float *column)
{
    v_set(m_updatePitches, epsilon, m_noteCount);

    for (int i = 0; i < m_sourceCount; ++i) {
        v_set(m_updateSources[i], epsilon, m_noteCount);
    }

    if (m_shifts) {
        for (int i = 0; i < m_shiftCount; ++i) {
            v_set(m_updateShifts[i], epsilon, m_noteCount);
        }
    }

    float *contributions = allocate<float>(m_binCount);

    for (int f = 0; f < m_shiftCount; ++f) {

        const float *shiftIn = m_shifts ? m_shifts[f] : 0;
        float *shiftOut = m_shifts ? m_updateShifts[f] : 0;

        for (int i = 0; i < m_sourceCount; ++i) {

            const float *sourceIn = m_sources[i];
            float *sourceOut = m_updateSources[i];

            int lowest, highest;
            rangeFor(i, lowest, highest);

            for (int n = lowest; n <= highest; ++n) {

                const float shift = shiftIn ? shiftIn[n] : 1.0;
                const float source = sourceIn[n];
                const float pitch = m_pitches[n];

                const float factor = pitch * source * shift;
                const float *w = templateFor(i, n, f);

                v_copy(contributions, w, m_binCount);
                v_multiply(contributions, m_q, m_binCount);

                float total = factor * v_sum(contributions, m_binCount);

                m_updatePitches[n] += total;
                sourceOut[n] += total;

                if (shiftOut) {
                    shiftOut[n] += total;
                }
            }
        }
    }

    if (m_pitchSparsity != 1.0) {
        for (int n = 0; n < m_noteCount; ++n) {
            m_updatePitches[n] = 
                pow(m_updatePitches[n], m_pitchSparsity);
        }
    }

    if (m_sourceSparsity != 1.0) {
        for (int i = 0; i < m_sourceCount; ++i) {
            for (int n = 0; n < m_noteCount; ++n) {
                m_updateSources[i][n] =
                    pow(m_updateSources[i][n], m_sourceSparsity);
            }
        }
    }

    normaliseColumn(m_updatePitches, m_noteCount);
    std::swap(m_pitches, m_updatePitches);

    normaliseGrid(m_updateSources, m_sourceCount, m_noteCount);
    std::swap(m_sources, m_updateSources);

    if (m_shifts) {
        normaliseGrid(m_updateShifts, m_shiftCount, m_noteCount);
        std::swap(m_shifts, m_updateShifts);
    }
}


