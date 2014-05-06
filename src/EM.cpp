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

#include <vector>

using std::vector;
using std::cerr;
using std::endl;

static double epsilon = 1e-16;

EM::EM() :
    m_noteCount(SILVET_TEMPLATE_NOTE_COUNT),
    m_shiftCount(SILVET_TEMPLATE_MAX_SHIFT * 2 + 1),
    m_binCount(SILVET_TEMPLATE_HEIGHT),
    m_instrumentCount(SILVET_TEMPLATE_COUNT),
    m_pitchSparsity(1.1),
    m_sourceSparsity(1.3),
    m_lowestPitch(silvet_templates_lowest_note),
    m_highestPitch(silvet_templates_highest_note)
{
    m_pitches = V(m_noteCount);
    for (int n = 0; n < m_noteCount; ++n) {
        m_pitches[n] = drand48();
    }

    m_shifts = Grid(m_shiftCount);
    for (int f = 0; f < m_shiftCount; ++f) {
        m_shifts[f] = V(m_noteCount);
        for (int n = 0; n < m_noteCount; ++n) {
            m_shifts[f][n] = drand48();
        }
    }
    
    m_sources = Grid(m_instrumentCount);
        for (int i = 0; i < m_instrumentCount; ++i) {
        m_sources[i] = V(m_noteCount);
        for (int n = 0; n < m_noteCount; ++n) {
            m_sources[i][n] = (inRange(i, n) ? 1.0 : 0.0);
        }
    }

    m_estimate = V(m_binCount);
    m_q = V(m_binCount);
}

EM::~EM()
{
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
EM::normaliseColumn(V &column)
{
    double sum = 0.0;
    for (int i = 0; i < (int)column.size(); ++i) {
        sum += column[i];
    }
    for (int i = 0; i < (int)column.size(); ++i) {
        column[i] /= sum;
    }
}

void
EM::normaliseGrid(Grid &grid)
{
    V denominators(grid[0].size());

    for (int i = 0; i < (int)grid.size(); ++i) {
        for (int j = 0; j < (int)grid[i].size(); ++j) {
            denominators[j] += grid[i][j];
        }
    }

    for (int i = 0; i < (int)grid.size(); ++i) {
        for (int j = 0; j < (int)grid[i].size(); ++j) {
            grid[i][j] /= denominators[j];
        }
    }
}

void
EM::iterate(V column)
{
    normaliseColumn(column);
    expectation(column);
    maximisation(column);
}

const float *
EM::templateFor(int instrument, int note, int shift)
{
    return silvet_templates[instrument].data[note] + shift;
}

void
EM::expectation(const V &column)
{
//    cerr << ".";

    for (int i = 0; i < m_binCount; ++i) {
        m_estimate[i] = epsilon;
    }

    for (int i = 0; i < m_instrumentCount; ++i) {
        for (int n = 0; n < m_noteCount; ++n) {
            const double pitch = m_pitches[n];
            const double source = m_sources[i][n];
            for (int f = 0; f < m_shiftCount; ++f) {
                const float *w = templateFor(i, n, f);
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
EM::maximisation(const V &column)
{
    V newPitches = m_pitches;
    Grid newShifts = m_shifts;
    Grid newSources = m_sources;

    for (int n = 0; n < m_noteCount; ++n) {

        const double pitch = m_pitches[n];
        newPitches[n] = epsilon;

        for (int f = 0; f < m_shiftCount; ++f) {

            const double shift = m_shifts[f][n];
            newShifts[f][n] = epsilon;

            for (int i = 0; i < m_instrumentCount; ++i) {

                const double source = m_sources[i][n];
                newSources[i][n] = epsilon;

                const float *w = templateFor(i, n, f);
                const double factor = pitch * source * shift;

                if (n >= m_lowestPitch && n <= m_highestPitch) {

                    for (int j = 0; j < m_binCount; ++j) {
                        newPitches[n] += w[j] * m_q[j] * factor;
                    }

                    if (inRange(i, n)) {
                        for (int j = 0; j < m_binCount; ++j) {
                            newSources[i][n] += w[j] * m_q[j] * factor;
                        }
                    }
                }

                for (int j = 0; j < m_binCount; ++j) {
                    newShifts[f][n] += w[j] * m_q[j] * factor;
                }
            }
        }
    }

    for (int n = 0; n < m_noteCount; ++n) {
        if (m_pitchSparsity != 1.0) {
            newPitches[n] = pow(newPitches[n], m_pitchSparsity);
        }
        if (m_sourceSparsity != 1.0) {
            for (int i = 0; i < m_instrumentCount; ++i) {
                newSources[i][n] = pow(newSources[i][n], m_sourceSparsity);
            }
        }
    }

    normaliseColumn(newPitches);
    normaliseGrid(newShifts);
    normaliseGrid(newSources);

    m_pitches = newPitches;
    m_shifts = newShifts;
    m_sources = newSources;
}


