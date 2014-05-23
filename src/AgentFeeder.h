/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  Silvet

  A Vamp plugin for note transcription.
  Centre for Digital Music, Queen Mary University of London.
  This file Copyright 2012 Chris Cannam.
    
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of the
  License, or (at your option) any later version.  See the file
  COPYING included with this distribution for more information.
*/

#ifndef AGENT_FEEDER_H
#define AGENT_FEEDER_H

#include "AgentHypothesis.h"

class AgentFeeder
{
public:
    virtual void feed(AgentHypothesis::Observation) = 0;
    virtual void finish() = 0;

    virtual ~AgentFeeder() { }
};

#endif
