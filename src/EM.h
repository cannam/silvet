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

    void iterate(std::vector<double> column);

    const std::vector<double> &getEstimate() const { 
	return m_estimate;
    }
    const std::vector<double> &getPitchDistribution() const {
	return m_pitches;
    }
    const std::vector<std::vector<double> > &getSources() const {
	return m_sources; 
    }

private:
    typedef std::vector<double> V;
    typedef std::vector<std::vector<double> > Grid;

    V m_pitches;
    Grid m_shifts;
    Grid m_sources;

    V m_estimate;
    V m_q;
    
    const int m_noteCount;
    const int m_shiftCount; // 1 + 2 * max template shift
    const int m_binCount;
    const int m_instrumentCount;
    
    const double m_pitchSparsity;
    const double m_sourceSparsity;

    const int m_lowestPitch;
    const int m_highestPitch;

    void normaliseColumn(V &column);
    void normaliseGrid(Grid &grid);
    void expectation(const V &column);
    void maximisation(const V &column);

    const float *templateFor(int instrument, int note, int shift);
    void rangeFor(int instrument, int &minPitch, int &maxPitch);
    bool inRange(int instrument, int pitch);
};

#endif
