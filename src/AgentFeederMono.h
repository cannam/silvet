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

#ifndef AGENT_FEEDER_MONO_H
#define AGENT_FEEDER_MONO_H

#include "AgentFeeder.h"

//#define DEBUG_FEEDER 1

/**
 * Take a series of observations or estimates (one at a time) and feed
 * them to a set of agent hypotheses, creating a new candidate agent
 * for each observation and also testing the observation against the
 * existing set of hypotheses.
 *
 * One satisfied hypothesis is considered to be "accepted" at any
 * moment (that is, the earliest contemporary hypothesis to have
 * become satisfied). The series of accepted and completed hypotheses
 * from construction to the present time can be queried through
 * retrieveAcceptedHypotheses(), which detaches them from the current
 * feeder (i.e. hypotheses returned from one call to this function
 * will not be returned by any subsequent call).
 *
 * Call feed() to provide a new observation. Call finish() when all
 * observations have been provided. The set of hypotheses returned by
 * getAcceptedHypotheses() will not be complete unless finish() has
 * been called.
 */
template <typename Hypothesis>
class AgentFeederMono : public AgentFeeder
{
private:
    typedef std::vector<Hypothesis> Hypotheses;

public:
    AgentFeederMono() : m_haveCurrent(false) { }

    virtual void feed(AgentHypothesis::Observation o) {

#ifdef DEBUG_FEEDER        
        std::cerr << "feed: have observation [value = " << o.value << ", time = " << o.time << "]" << std::endl;
#endif

        if (m_haveCurrent) {
            if (m_current.accept(o)) {
                return;
            }
            if (m_current.getState() == Hypothesis::Expired) {
                m_accepted.push_back(m_current);
#ifdef DEBUG_FEEDER        
                std::cerr << "current has expired, pushing to accepted" << std::endl;
#endif
                m_haveCurrent = false;
            }
        }

        bool swallowed = false;

#ifdef DEBUG_FEEDER        
        std::cerr << "not swallowed by current" << std::endl;
#endif

        Hypotheses newCandidates;

        for (typename Hypotheses::iterator i = m_candidates.begin();
             i != m_candidates.end(); ++i) {

            Hypothesis h = *i;
                
            if (swallowed) {

                // don't offer: each observation can only belong to one
                // satisfied hypothesis
                newCandidates.push_back(h);

            } else {

                if (h.accept(o)) {
#ifdef DEBUG_FEEDER        
                    std::cerr << "accepted, state is " << h.getState() << std::endl;
#endif
                    if (h.getState() == Hypothesis::Satisfied) {

                        swallowed = true;
        
                        if (!m_haveCurrent ||
                            m_current.getState() == Hypothesis::Expired ||
                            m_current.getState() == Hypothesis::Rejected) {
#ifdef DEBUG_FEEDER        
                            std::cerr << "current has ended, updating from candidate" << std::endl;
#endif
                            m_current = h;
                            m_haveCurrent = true;
                        } else {
                            newCandidates.push_back(h);
                        }

                    } else {
                        newCandidates.push_back(h);
                    }
                }
            }
        }
            
        if (!swallowed) {
            Hypothesis h;
            h.accept(o); // must succeed, as h is new
            newCandidates.push_back(h);
#ifdef DEBUG_FEEDER        
            std::cerr << "not swallowed, creating new hypothesis" << std::endl;
#endif
        }

        // reap rejected/expired hypotheses from candidates list,
        // and assign back to m_candidates
        
        m_candidates.clear();

        for (typename Hypotheses::const_iterator i = newCandidates.begin();
             i != newCandidates.end(); ++i) {
            Hypothesis h = *i;
            if (h.getState() != Hypothesis::Rejected && 
                h.getState() != Hypothesis::Expired) {
                m_candidates.push_back(h);
            } else {
#ifdef DEBUG_FEEDER        
                std::cerr << "reaping a candidate" << std::endl;
#endif
            }
        }  
#ifdef DEBUG_FEEDER        
        std::cerr << "have " << m_candidates.size() << " candidates" << std::endl;
#endif
    }

    virtual void finish() {
        if (m_haveCurrent &&
            (m_current.getState() == Hypothesis::Satisfied)) {
#ifdef DEBUG_FEEDER        
            std::cerr << "finish: current is satisfied, pushing to accepted" << std::endl;
#endif
            m_accepted.push_back(m_current);
        }
    }

    std::set<Hypothesis> retrieveAcceptedHypotheses() {
        std::set<Hypothesis> hs;
#ifdef DEBUG_FEEDER
        std::cerr << "retrieveAcceptedHypotheses: returning " << m_accepted.size() << " accepted" << std::endl;
#endif
        for (typename Hypotheses::const_iterator i = m_accepted.begin();
             i != m_accepted.end(); ++i) {
            hs.insert(*i);
        }
        m_accepted.clear();
        return hs;
    }

private:
    Hypotheses m_candidates;
    Hypothesis m_current;
    bool m_haveCurrent;
    Hypotheses m_accepted;
};

#endif
