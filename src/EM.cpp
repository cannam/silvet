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
    m_pitchCount(m_noteCount * m_shiftCount),
    m_binCount(SILVET_TEMPLATE_HEIGHT),
    m_instrumentCount(SILVET_TEMPLATE_COUNT),
    m_pitchSparsity(1.1),
    m_sourceSparsity(1.3)
{
    m_lowestPitch = 
        silvet_templates_lowest_note * m_shiftCount;
    m_highestPitch =
        silvet_templates_highest_note * m_shiftCount + m_shiftCount - 1;

    m_pitches = V(m_pitchCount);

    for (int n = 0; n < m_pitchCount; ++n) {
        m_pitches[n] = drand48();
    }
    
    m_sources = Grid(m_instrumentCount);
    
    for (int i = 0; i < m_instrumentCount; ++i) {
        m_sources[i] = V(m_pitchCount);
        for (int n = 0; n < m_pitchCount; ++n) {
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
    minPitch = silvet_templates[instrument].lowest * m_shiftCount;
    maxPitch = silvet_templates[instrument].highest * m_shiftCount
        + m_shiftCount - 1;
}

bool
EM::inRange(int instrument, int pitch)
{
    int minPitch, maxPitch;
    rangeFor(instrument, minPitch, maxPitch);
    return (pitch >= minPitch && pitch <= maxPitch);
}

void
EM::normalise(V &column)
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
EM::iterate(V column)
{
    normalise(column);
    expectation(column);
    maximisation(column);
}

const float *
EM::templateFor(int instrument, int pitch)
{
    int note = pitch / m_shiftCount;
    int shift = pitch % m_shiftCount;
    return silvet_templates[instrument].data[note] + shift;
}

void
EM::expectation(const V &column)
{
    cerr << ".";

    for (int i = 0; i < m_binCount; ++i) {
        m_estimate[i] = epsilon;
    }

    for (int i = 0; i < m_instrumentCount; ++i) {
        for (int n = 0; n < m_pitchCount; ++n) {
            const float *w = templateFor(i, n);
            double pitch = m_pitches[n];
            double source = m_sources[i][n];
            for (int j = 0; j < m_binCount; ++j) {
                m_estimate[j] += w[j] * pitch * source;
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

    for (int n = 0; n < m_pitchCount; ++n) {
        newPitches[n] = epsilon;
        if (n >= m_lowestPitch && n <= m_highestPitch) {
            for (int i = 0; i < m_instrumentCount; ++i) {
                const float *w = templateFor(i, n);
                double pitch = m_pitches[n];
                double source = m_sources[i][n];
                for (int j = 0; j < m_binCount; ++j) {
                    newPitches[n] += w[j] * m_q[j] * pitch * source;
                }
            }
        }
        if (m_pitchSparsity != 1.0) {
            newPitches[n] = pow(newPitches[n], m_pitchSparsity);
        }
    }
    normalise(newPitches);

    Grid newSources = m_sources;

    for (int i = 0; i < m_instrumentCount; ++i) {
        for (int n = 0; n < m_pitchCount; ++n) {
            newSources[i][n] = epsilon;
            if (inRange(i, n)) {
                const float *w = templateFor(i, n);
                double pitch = m_pitches[n];
                double source = m_sources[i][n];
                for (int j = 0; j < m_binCount; ++j) {
                    newSources[i][n] += w[j] * m_q[j] * pitch * source;
                }
            }
            if (m_sourceSparsity != 1.0) {
                newSources[i][n] = pow(newSources[i][n], m_sourceSparsity);
            }
        }
        normalise(newSources[i]);
    }

    m_pitches = newPitches;
    m_sources = newSources;
}

void
EM::report()
{
    vector<int> sounding;
    for (int n = 0; n < m_pitchCount; ++n) {
        if (m_pitches[n] > 0.05) {
            sounding.push_back(n);
        }
    }
    cerr << " sounding: ";
    for (int i = 0; i < (int)sounding.size(); ++i) {
        cerr << sounding[i] << " ";
        int maxj = -1;
        double maxs = 0.0;
        for (int j = 0; j < m_instrumentCount; ++j) {
            if (j == 0 || m_sources[j][sounding[i]] > maxs) {
                maxj = j;
                maxs = m_sources[j][sounding[i]];
            }
        }
        cerr << silvet_templates[maxj].name << " ";
    }
    cerr << endl;
}

