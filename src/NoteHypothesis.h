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
        Note() : freq(0), time(), duration() { }
        Note(double _f, Vamp::RealTime _t, Vamp::RealTime _d) :
            freq(_f), time(_t), duration(_d) { }
        bool operator==(const Note &e) const {
            return e.freq == freq && e.time == time && e.duration == duration;
        }
	double freq;
	Vamp::RealTime time;
	Vamp::RealTime duration;
    };
    
    /**
     * Return the mean frequency of the accepted observations
     */
    double getMeanFrequency() const;

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

    /**
     * Convert the given sequence of accepted hypotheses into a
     * sampled series of pitches (in Hz), returned at regular
     * intervals determined by the given start time, end time, and
     * interval. The range is [start,end], i.e. the end time is
     * included. The interval must be greater than zero.
     *
     * Unvoiced samples are returned as 0Hz.
     */
    static std::vector<double> sample(const std::set<NoteHypothesis> &,
                                      Vamp::RealTime startTime,
                                      Vamp::RealTime endTime,
                                      Vamp::RealTime interval);

    /**
     *!!! No! Not equally spaced, should be able to be anything [ordered]
     
     * Given a series of equally spaced observations, return a series
     * of the same number of pitches (in Hz) calculated by running an
     * AgentFeeder<NoteHypothesis> on the observations and flattening
     * and sampling the resulting accepted hypotheses.
     *
     * The result should contain only pitches that contributed to
     * recognised notes in the input observations, with the remaining
     * (unvoiced) samples returned as 0Hz.
     *
     * If the input observations are not equally spaced, the result is
     * undefined.
     *!!! (what about rounding errors from RealTime to frame and vice versa?)
     *!!! (should provide a Timebase?)
     *!!! update docs for updated api
     */
    static std::vector<double> winnow(const Observations &,
                                      Vamp::RealTime startTime,
                                      Vamp::RealTime endTime,
                                      Vamp::RealTime interval);

    //!!!
    bool operator==(const NoteHypothesis &other) const {
        return m_state == other.m_state && m_pending == other.m_pending;
    }

    bool operator<(const NoteHypothesis &other) const {
        return getStartTime() < other.getStartTime();
    }

private:
    bool isWithinTolerance(Observation) const;
    bool isOutOfDateFor(Observation) const;
    bool isSatisfied() const;
    
    State m_state;
    Observations m_pending;
};

#endif
