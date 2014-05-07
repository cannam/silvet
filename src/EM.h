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
    EM(bool useShifts);
    ~EM();

    int getBinCount() const { return m_binCount; }
    int getNoteCount() const { return m_noteCount; }
    int getSourceCount() const { return m_sourceCount; }

    /**
     * Carry out one iteration using the given column as input. The
     * column must have getBinCount() values.
     */
    void iterate(const double *column);

    /**
     * Return the estimated distribution after the current iteration.
     * Like the input, this will have getBinCount() values.
     */
    const double *getEstimate() const {
	return m_estimate;
    }

    /**
     * Return the pitch distribution for the current estimate.  The
     * returned array has getNoteCount() values.
     */
    const double *getPitchDistribution() const {
	return m_pitches;
    }
    
    /** 
     * Return the source distribution for the current estimate. The
     * returned pointer refers to getSourceCount() arrays of
     * getNoteCount() values.
     */
    const double *const *getSources() const {
	return m_sources; 
    }

private:
    double *m_pitches;
    double **m_shifts;
    double **m_sources;

    double *m_updatePitches;
    double **m_updateShifts;
    double **m_updateSources;

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
    void normaliseGrid(double **grid, int size1, int size2);

    void expectation(const double *column); // size is m_binCount
    void maximisation(const double *column); // size is m_binCount

    const double *templateFor(int instrument, int note, int shift);
    void rangeFor(int instrument, int &minPitch, int &maxPitch);
    bool inRange(int instrument, int pitch);
};

#endif
