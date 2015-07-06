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

#ifndef SILVET_LIVE_INSTRUMENTS_H
#define SILVET_LIVE_INSTRUMENTS_H

#include "Instruments.h"

/**
 * Adapt an instrument pack into a "live" version, with fewer bins per
 * octave and so lower CQ latency.
 */
class LiveAdapter
{
public:
    static InstrumentPack adapt(const InstrumentPack &original);
    static std::vector<InstrumentPack> adaptAll(const std::vector<InstrumentPack> &);
};

#endif
