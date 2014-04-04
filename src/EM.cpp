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

    for (int n = 0; n < m_notes; ++i) {
        m_pitches[n] = drand48();
    }
    
    m_sources = Grid(m_instruments);
    
    for (int i = 0; i < m_instruments; ++i) {
        m_sources[i] = V(m_notes);
        for (int n = 0; n < m_notes; ++n) {
            m_sources[i][n] = (inRange(i, n) ? 1.0 : 0.0);
        }
    }

    m_q = V(m_bins);
    
    for (int w = 0; w < m_bins; ++w) {
        m_q[w] = epsilon;
    }
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

