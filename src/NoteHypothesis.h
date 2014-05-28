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

#ifndef NOTE_HYPOTHESIS_H
#define NOTE_HYPOTHESIS_H

#include "AgentHypothesis.h"

#include <set>
#include <vector>

/**
 * An AgentHypothesis which tests a series of instantaneous pitch
 * estimates to see whether they fit a single-note relationship.
 * Contains rules specific to testing note pitch and timing.
 */

class NoteHypothesis : public AgentHypothesis
{
public:
    /**
     * Construct an empty hypothesis. This will be in New state and
     * will provisionally accept any estimate.
     */
    NoteHypothesis();

    /**
     * Destroy the hypothesis
     */
    ~NoteHypothesis();

    virtual bool accept(Observation);
    virtual State getState() const;
    virtual Observations getAcceptedObservations() const;

    struct Note {
        Note() : freq(0), time(), duration(), confidence(1.0) { }
        Note(double _f, Vamp::RealTime _t, Vamp::RealTime _d, double _i) :
            freq(_f), time(_t), duration(_d), confidence(_i) { }
        bool operator==(const Note &e) const {
            return e.freq == freq && e.time == time && 
                e.duration == duration && e.confidence == confidence;
        }
	double freq;
	Vamp::RealTime time;
	Vamp::RealTime duration;
        double confidence;
    };
    
    /**
     * Return the mean frequency of the accepted observations
     */
    double getMeanFrequency() const;
    
    /**
     * Return the median frequency of the accepted observations
     */
    double getMedianFrequency() const;
    
    /**
     * Return the median confidence of the accepted observations
     */
    double getMedianConfidence() const;

    /**
     * Return a single note roughly matching this hypothesis
     */
    Note getAveragedNote() const;
    
    /**
     * Return the time of the first accepted observation
     */
    Vamp::RealTime getStartTime() const;

    /**
     * Return the difference between the start time and the end of the
     * final accepted observation
     */
    Vamp::RealTime getDuration() const;

    //!!!
    bool operator==(const NoteHypothesis &other) const {
        return m_state == other.m_state && m_pending == other.m_pending;
    }

    bool operator<(const NoteHypothesis &other) const {
        if (getStartTime() != other.getStartTime()) {
            return getStartTime() < other.getStartTime();
        } else if (m_state != other.m_state) {
            return m_state < other.m_state;
        } else if (m_pending.size() != other.m_pending.size()) {
            return m_pending.size() < other.m_pending.size();
        } else {
            Observations::const_iterator i = m_pending.begin();
            Observations::const_iterator j = other.m_pending.begin();
            while (i != m_pending.end()) {
                if (*i == *j) {
                    ++i;
                    ++j;
                } else {
                    return *i < *j;
                }
            }
            return false;
        }
    }

private:
    bool isWithinTolerance(Observation) const;
    bool isOutOfDateFor(Observation) const;
    bool isSatisfied() const;
    
    State m_state;
    Observations m_pending;
};

#endif
