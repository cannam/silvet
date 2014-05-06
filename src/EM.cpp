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

#include "data/include/templates.h"

#include <cstdlib>
#include <cmath>

#include <iostream>

#include "VectorOps.h"
#include "Allocators.h"

using std::vector;
using std::cerr;
using std::endl;

using namespace breakfastquay;

static double epsilon = 1e-16;

EM::EM() :
    m_noteCount(SILVET_TEMPLATE_NOTE_COUNT),
    m_shiftCount(SILVET_TEMPLATE_MAX_SHIFT * 2 + 1),
    m_binCount(SILVET_TEMPLATE_HEIGHT),
    m_sourceCount(SILVET_TEMPLATE_COUNT),
    m_pitchSparsity(1.1),
    m_sourceSparsity(1.3),
    m_lowestPitch(silvet_templates_lowest_note),
    m_highestPitch(silvet_templates_highest_note)
{
    m_pitches = allocate<double>(m_noteCount);
    for (int n = 0; n < m_noteCount; ++n) {
        m_pitches[n] = drand48();
    }

    m_shifts = allocate_channels<double>(m_shiftCount, m_noteCount);
    for (int f = 0; f < m_shiftCount; ++f) {
        for (int n = 0; n < m_noteCount; ++n) {
            m_shifts[f][n] = drand48();
        }
    }
    
    m_sources = allocate_channels<double>(m_sourceCount, m_noteCount);
    for (int i = 0; i < m_sourceCount; ++i) {
        for (int n = 0; n < m_noteCount; ++n) {
            m_sources[i][n] = (inRange(i, n) ? 1.0 : 0.0);
        }
    }

    m_estimate = allocate<double>(m_binCount);
    m_q = allocate<double>(m_binCount);
}

EM::~EM()
{
    deallocate(m_q);
    deallocate(m_estimate);
    deallocate_channels(m_sources, m_sourceCount);
    deallocate_channels(m_shifts, m_shiftCount);
    deallocate(m_pitches);
}

void
EM::rangeFor(int instrument, int &minPitch, int &maxPitch)
{
    minPitch = silvet_templates[instrument].lowest;
    maxPitch = silvet_templates[instrument].highest;
}

bool
EM::inRange(int instrument, int pitch)
{
    int minPitch, maxPitch;
    rangeFor(instrument, minPitch, maxPitch);
    return (pitch >= minPitch && pitch <= maxPitch);
}

void
EM::normaliseColumn(double *column, int size)
{
    double sum = v_sum(column, size);
    v_scale(column, 1.0 / sum, size);
}

void
EM::normaliseGrid(double **grid, int size1, int size2)
{
    double *denominators = allocate_and_zero<double>(size2);

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
    double *norm = allocate<double>(m_binCount);
    v_copy(norm, column, m_binCount);
    normaliseColumn(norm, m_binCount);
    expectation(norm);
    maximisation(norm);
}

const double *
EM::templateFor(int instrument, int note, int shift)
{
    return silvet_templates[instrument].data[note] + shift;
}

void
EM::expectation(const double *column)
{
//    cerr << ".";

    for (int i = 0; i < m_binCount; ++i) {
        m_estimate[i] = epsilon;
    }

    for (int i = 0; i < m_sourceCount; ++i) {
        for (int n = 0; n < m_noteCount; ++n) {
            const double pitch = m_pitches[n];
            const double source = m_sources[i][n];
            for (int f = 0; f < m_shiftCount; ++f) {
                const double *w = templateFor(i, n, f);
                const double shift = m_shifts[f][n];
                const double factor = pitch * source * shift;
                for (int j = 0; j < m_binCount; ++j) {
                    m_estimate[j] += w[j] * factor;
                }
            }
        }
    }

    for (int i = 0; i < m_binCount; ++i) {
        m_q[i] = column[i] / m_estimate[i];
    }
}

void
EM::maximisation(const double *column)
{
    double *newPitches = allocate<double>(m_noteCount);
    v_set(newPitches, epsilon, m_noteCount);

    double **newShifts = allocate_channels<double>(m_shiftCount, m_noteCount);
    for (int i = 0; i < m_shiftCount; ++i) {
        v_set(newShifts[i], epsilon, m_noteCount);
    }

    double **newSources = allocate_channels<double>(m_sourceCount, m_noteCount);
    for (int i = 0; i < m_sourceCount; ++i) {
        v_set(newSources[i], epsilon, m_noteCount);
    }

    double *contributions = allocate<double>(m_binCount);

    for (int n = 0; n < m_noteCount; ++n) {

        const double pitch = m_pitches[n];

        for (int f = 0; f < m_shiftCount; ++f) {

            const double shift = m_shifts[f][n];

            for (int i = 0; i < m_sourceCount; ++i) {

                const double source = m_sources[i][n];
                const double factor = pitch * source * shift;
                const double *w = templateFor(i, n, f);

                v_copy(contributions, w, m_binCount);
                v_add(contributions, m_q, m_binCount);
                v_scale(contributions, factor, m_binCount);

                double total = v_sum(contributions, m_binCount);

                if (n >= m_lowestPitch && n <= m_highestPitch) {

                    newPitches[n] += total;

                    if (inRange(i, n)) {
                        newSources[i][n] += total;
                    }
                }

                newShifts[f][n] += total;
            }
        }
    }

    for (int n = 0; n < m_noteCount; ++n) {
        if (m_pitchSparsity != 1.0) {
            newPitches[n] = pow(newPitches[n], m_pitchSparsity);
        }
        if (m_sourceSparsity != 1.0) {
            for (int i = 0; i < m_sourceCount; ++i) {
                newSources[i][n] = pow(newSources[i][n], m_sourceSparsity);
            }
        }
    }

    normaliseColumn(newPitches, m_noteCount);
    normaliseGrid(newShifts, m_shiftCount, m_noteCount);
    normaliseGrid(newSources, m_sourceCount, m_noteCount);

    deallocate(m_pitches);
    deallocate_channels(m_shifts, m_shiftCount);
    deallocate_channels(m_sources, m_sourceCount);

    m_pitches = newPitches;
    m_shifts = newShifts;
    m_sources = newSources;
}


