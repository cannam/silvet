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

#include <iostream>

#include <vector>

using std::vector;
using std::cerr;
using std::endl;

static double epsilon = 1e-16;

EM::EM() :
    m_notes(SILVET_TEMPLATE_NOTE_COUNT),
    m_bins(SILVET_TEMPLATE_HEIGHT),
    m_instruments(SILVET_TEMPLATE_COUNT)
{
    m_lowest = 0;
    m_highest = m_notes - 1;

    for (int i = 0; i < m_instruments; ++i) {
        if (i == 0 || silvet_templates[i].lowest < m_lowest) {
            m_lowest = silvet_templates[i].lowest;
        }
        if (i == 0 || silvet_templates[i].highest > m_highest) {
            m_highest = silvet_templates[i].highest;
        }
    }

    m_pitches = V(m_notes);

    for (int n = 0; n < m_notes; ++n) {
        m_pitches[n] = drand48();
    }
    
    m_sources = Grid(m_instruments);
    
    for (int i = 0; i < m_instruments; ++i) {
        m_sources[i] = V(m_notes);
        for (int n = 0; n < m_notes; ++n) {
            m_sources[i][n] = (inRange(i, n) ? 1.0 : 0.0);
        }
    }

    m_estimate = V(m_bins);
    m_q = V(m_bins);
}

EM::~EM()
{
}

bool
EM::inRange(int instrument, int note)
{
    return (note >= silvet_templates[instrument].lowest &&
            note <= silvet_templates[instrument].highest);
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

void
EM::expectation(const V &column)
{
    cerr << ".";

    for (int i = 0; i < m_bins; ++i) {
        m_estimate[i] = epsilon;
    }

    for (int i = 0; i < m_instruments; ++i) {
        for (int n = 0; n < m_notes; ++n) {
            float *w = silvet_templates[i].data[n];
            double pitch = m_pitches[n];
            double source = m_sources[i][n];
            for (int j = 0; j < m_bins; ++j) {
                m_estimate[j] += w[j] * pitch * source;
            }
        }
    }

    for (int i = 0; i < m_bins; ++i) {
        m_q[i] = column[i] / m_estimate[i];
    }
}

void
EM::maximisation(const V &column)
{
    V newPitches = m_pitches;

    for (int n = 0; n < m_notes; ++n) {
        newPitches[n] = epsilon;
        if (n >= m_lowest && n <= m_highest) {
            for (int i = 0; i < m_instruments; ++i) {
                float *w = silvet_templates[i].data[n];
                double pitch = m_pitches[n];
                double source = m_sources[i][n];
                for (int j = 0; j < m_bins; ++j) {
                    newPitches[n] += w[j] * m_q[j] * pitch * source;
                }
            }
        }
    }
    normalise(newPitches);

    Grid newSources = m_sources;

    for (int i = 0; i < m_instruments; ++i) {
        for (int n = 0; n < m_notes; ++n) {
            newSources[i][n] = epsilon;
            if (inRange(i, n)) {
                float *w = silvet_templates[i].data[n];
                double pitch = m_pitches[n];
                double source = m_sources[i][n];
                for (int j = 0; j < m_bins; ++j) {
                    newSources[i][n] += w[j] * m_q[j] * pitch * source;
                }
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
    for (int n = 0; n < m_notes; ++n) {
        if (m_pitches[n] > 0.05) {
            sounding.push_back(n);
        }
    }
    cerr << " sounding: ";
    for (int i = 0; i < (int)sounding.size(); ++i) {
        cerr << sounding[i] << " ";
        int maxj = -1;
        double maxs = 0.0;
        for (int j = 0; j < m_instruments; ++j) {
            if (j == 0 || m_sources[j][sounding[i]] > maxs) {
                maxj = j;
                maxs = m_sources[j][sounding[i]];
            }
        }
        cerr << silvet_templates[maxj].name << " ";
    }
    cerr << endl;
}

