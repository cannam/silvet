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

#ifndef SILVET_EM_H
#define SILVET_EM_H

#include <vector>

class EM
{
public:
    EM();
    ~EM();

    void iterate(const std::vector<double> &column);

private:
    typedef std::vector<double> V;
    typedef std::vector<std::vector<double> > Grid;

    V m_pitches;
    Grid m_sources;
    Grid m_q;
    
    int m_lowest;
    int m_highest;
};

#endif
