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

#ifndef AGENT_HYPOTHESIS_H
#define AGENT_HYPOTHESIS_H

#include "vamp-sdk/RealTime.h"

#include <set>
#include <map>

/**
 * An agent used to test an incoming series of timed observations or
 * estimates to see whether they fit a consistent single-object
 * relationship.
 *
 * A freshly constructed hypothesis object should be in New state and
 * should accept any observation.
 */

class AgentHypothesis
{
public:
    virtual ~AgentHypothesis() { }

    enum State {

	/// Just constructed, will provisionally accept any observation
	New,

	/// Accepted at least one observation, but not enough evidence to satisfy
	Provisional,

	/// Could not find enough consistency in offered observations
	Rejected,

	/// Have accepted enough consistent observations to satisfy hypothesis
	Satisfied,

	/// Have been satisfied, but evidence has now changed: we're done
	Expired
    };

    struct Observation {

        Observation() : value(0), time(), confidence(1) { }

        Observation(double _f, Vamp::RealTime _t, double _c) :
            value(_f), time(_t), confidence(_c) { }

        bool operator==(const Observation &o) const {
            return o.value == value && o.time == time && o.confidence == confidence;
        }
        bool operator<(const Observation &o) const {
            return
                (time < o.time || 
                 (time == o.time && value < o.value) ||
                 (time == o.time && value == o.value && confidence < o.confidence));
        }

	double value;
        Vamp::RealTime time;
	double confidence;
    };
    typedef std::set<Observation> Observations;

    /**
     * Test the given observation to see whether it is consistent with
     * this hypothesis, and adjust the hypothesis' internal state
     * accordingly. If the observation is not inconsistent with the
     * hypothesis, return true.
     *!!! should be called e.g. test?
     */
    virtual bool accept(Observation) = 0;

    /**
     * Return the current state of this hypothesis.
     */
    virtual State getState() const = 0;

    /**
     * If the hypothesis has been satisfied (i.e. is in Satisfied or
     * Expired state), return the set of observations that it
     * accepted. Otherwise return an empty set
     */
    virtual Observations getAcceptedObservations() const = 0;

    /**
     * Convert the given set of accepted hypotheses (of type
     * subclassed from AgentHypothesis) into a flattened set of their
     * accepted observations.
     * 
     * That is, only one is included for at any given moment, so in
     * the case of overlapping hypotheses, the observations for the
     * earlier are taken until the next hypothesis begins and then the
     * latter's observations begin instead.
     *
     * (If there are gaps between hypotheses, the gaps remain in the
     * output.)
     */
    template <typename HypothesisType>
    static Observations flatten(const std::set<HypothesisType> &agents) {

        typedef typename std::set<HypothesisType>::const_iterator Itr;
        Observations flattened;

        if (agents.empty()) return flattened;
        Observations obs = agents.begin()->getAcceptedObservations();

        for (Itr i = agents.begin(); i != agents.end(); ++i) {

            Itr j = i;
            ++j;

            Observations nextObs;
            if (j != agents.end()) nextObs = j->getAcceptedObservations();

            for (Observations::const_iterator i = obs.begin();
                 i != obs.end(); ++i) {
                if (!nextObs.empty() && i->time >= nextObs.begin()->time) {
                    break;
                }
                flattened.insert(*i);
            }

            obs = nextObs;
        }

        return flattened;
    }
};

#endif
