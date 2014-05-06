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

    int getBinCount() const { return m_binCount; } // size of input column
    int getNoteCount() const { return m_noteCount; } // size of pitch column
    int getSourceCount() const { return m_sourceCount; }

    void iterate(double *column);

    const double *getEstimate() const { // bin count
	return m_estimate;
    }
    const double *getPitchDistribution() const { // note count
	return m_pitches;
    }
    const double **getSources() const { // source count * note count
	return m_sources; 
    }

private:
    double *m_pitches;
    double **m_shifts;
    double **m_sources;

    double *m_estimate;
    double *m_q;
    
    const int m_noteCount;
    const int m_shiftCount; // 1 + 2 * max template shift
    const int m_binCount;
    const int m_sourceCount;
    
    const double m_pitchSparsity;
    const double m_sourceSparsity;

    const int m_lowestPitch;
    const int m_highestPitch;

    void normaliseColumn(double *column, int size);
    void normaliseGrid(double **grid, int width, int height);

    void expectation(double *column); // size is m_binCount
    void maximisation(double *column); // size is m_binCount

    const double *templateFor(int instrument, int note, int shift);
    void rangeFor(int instrument, int &minPitch, int &maxPitch);
    bool inRange(int instrument, int pitch);
};

#endif
