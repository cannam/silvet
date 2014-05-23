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

#ifndef AGENT_FEEDER_POLY_H
#define AGENT_FEEDER_POLY_H

#include "AgentFeeder.h"

#include <cassert>

//#define DEBUG_FEEDER 1

/**
 * Take a series of observations or estimates (one at a time) and feed
 * them to a set of agent hypotheses, creating a new candidate agent
 * for each observation and also testing the observation against the
 * existing set of hypotheses.
 *
 *!!! -- todo: document poly-ness of it
 *
 * Call feed() to provide a new observation. Call finish() when all
 * observations have been provided. The set of hypotheses returned by
 * getAcceptedHypotheses() will not be complete unless finish() has
 * been called.
 */
template <typename Hypothesis>
class AgentFeederPoly : public AgentFeeder
{
public:
    typedef std::set<Hypothesis> Hypotheses;
    
private:
    struct State {
        Hypotheses provisional;
        Hypotheses satisfied;
        Hypotheses completed;
    };
    State m_state;

public:
    AgentFeederPoly() { }

    virtual void feed(AgentHypothesis::Observation o) {
#ifdef DEBUG_FEEDER        
        std::cerr << "feed: have observation [value = " << o.value << ", time = " << o.time << "]" << std::endl;
#endif

        update(m_state, o);
    }

    virtual void finish() {
#ifdef DEBUG_FEEDER        
        std::cerr << "finish: satisfied count == " << m_state.satisfied.size() << std::endl;
#endif
        for (typename Hypotheses::const_iterator i = m_state.satisfied.begin();
             i != m_state.satisfied.end(); ++i) {
            m_state.completed.insert(*i);
        }
    }

    Hypotheses getAcceptedHypotheses() const {
        return m_state.completed;
    }

private:
    void update(State &s, AgentHypothesis::Observation o) {

        /*
          An observation can "belong" to any number of provisional
          hypotheses, but only to one satisfied hypothesis.

          A new observation is first offered to the hypotheses that
          have already been satisfied. If one of these accepts it, it
          gets to keep it and no other hypothesis can have it.
         
          Any observation not accepted by a hypothesis in satisfied
          state is then offered to the provisional hypotheses; any
          number of these may accept it. Also, every observation that
          belongs to no satisfied hypothesis is used as the first
          observation in its own new hypothesis (regardless of how
          many other provisional hypotheses have accepted it).

          When a hypothesis subsequently becomes satisfied, all other
          provisional hypotheses containing any of its observations
          must be discarded.
        */

        // We only ever add to the completed hypotheses, never remove
        // anything from them. But we may remove from provisional (if
        // rejected or transferred to satisfied) and satisfied (when
        // completed).

        Hypotheses toCompleted;
        Hypotheses toSatisfied;
        Hypotheses toProvisional;
        
        bool swallowed = false;

        for (typename Hypotheses::iterator i = s.satisfied.begin();
             i != s.satisfied.end(); ++i) {

            Hypothesis h = *i;

            assert(h.getState() == Hypothesis::Satisfied);

            if (swallowed) {

                // An observation that has already been accepted by a
                // hypothesis cannot be offered to any other, because
                // it can only belong to one satisfied hypothesis. Any
                // subsequent satisfied hypotheses are retained
                // unchanged in our updated state. We can't test them
                // for expiry, because the state is only updated when
                // accept() is called.
                //!!! That looks like a limitation in the Hypothesis API

            } else { // !swallowed

                if (h.accept(o)) {
#ifdef DEBUG_FEEDER        
                    std::cerr << "accepted by satisfied hypothesis " << &(*i) << ", state is " << h.getState() << std::endl;
#endif
                    swallowed = true;
                } else if (h.getState() == Hypothesis::Expired) {
                    toCompleted.insert(h);
                }
            }
        }

        if (swallowed) {

#ifdef DEBUG_FEEDER        
            std::cerr << "was swallowed by satisfied hypothesis" << std::endl;
#endif
            // no provisional hypotheses have become satisfied, no new
            // ones have been introduced

        } else {

#ifdef DEBUG_FEEDER        
            std::cerr << "remained unswallowed by " << newState.satisfied.size() << " satisfied hypotheses" << std::endl;
#endif

            // not swallowed by any satisfied hypothesis

            Hypothesis promoted;
        
            for (typename Hypotheses::iterator i = s.provisional.begin();
                 i != s.provisional.end(); ++i) {

                Hypothesis h = *i;

                assert(h.getState() == Hypothesis::Provisional);

                // can only have one satisfied hypothesis for each
                // observation, so try this only if promoted has not been
                // set to something else yet
                if (promoted == Hypothesis() &&
                    h.accept(o) &&
                    h.getState() == Hypothesis::Satisfied) {
                    toSatisfied.insert(h);
#ifdef DEBUG_FEEDER        
                    std::cerr << "promoting a hypothesis to satisfied, have " << newState.satisfied.size() << " satisfied now" << std::endl;
#endif
                    promoted = h;
                } else if (h.getState() != Hypothesis::Rejected) {
                    // leave as provisional
                }
            }

            if (promoted == Hypothesis()) {

                // No provisional hypothesis has become satisfied

                Hypothesis h;
                h.accept(o);

                if (h.getState() == Hypothesis::Provisional) {
                    toProvisional.insert(h);
                } else if (h.getState() == Hypothesis::Satisfied) {
                    toSatisfied.insert(h);
                }

#ifdef DEBUG_FEEDER        
                std::cerr << "update: new hypothesis of state " << h.getState() << ", provisional count -> " << newState.provisional.size() << std::endl;
#endif

            }
        }

        for (typename Hypotheses::const_iterator i = toCompleted.begin();
             i != toCompleted.end(); ++i) {
            s.satisfied.erase(*i);
            s.completed.insert(*i);
        }
        for (typename Hypotheses::const_iterator i = toSatisfied.begin();
             i != toSatisfied.end(); ++i) {
            s.provisional.erase(*i);
            s.satisfied.insert(*i);
        }
        for (typename Hypotheses::const_iterator i = toProvisional.begin();
             i != toProvisional.end(); ++i) {
            s.provisional.insert(*i);
        }

        reap(s);
    }

    void reap(State &s) {

        // "When a hypothesis subsequently becomes satisfied, all
        // other provisional hypotheses containing any of its
        // observations must be discarded."

        if (s.provisional.empty()) return;

        int reaped = 0;

        Hypotheses toRemove = Hypotheses();;

        for (typename Hypotheses::const_iterator hi = s.provisional.begin();
             hi != s.provisional.end(); ++hi) {
                    
            const AgentHypothesis::Observations obs =
                hi->getAcceptedObservations();

            bool keep = true;

            for (AgentHypothesis::Observations::const_iterator oi = obs.begin();
                 oi != obs.end(); ++oi) {

                for (typename Hypotheses::const_iterator si = s.satisfied.end();
                     si != s.satisfied.begin(); ) {
                
                    --si;

                    const AgentHypothesis::Observations sobs = 
                        si->getAcceptedObservations();

                    if (sobs.find(*oi) != sobs.end()) {
                        keep = false;
                        break;
                    }
                }

                if (!keep) {
                    break;
                }
            }

            if (!keep) {
                toRemove.insert(*hi);
                ++reaped;
            }
        }
        
        for (typename Hypotheses::const_iterator i = toRemove.begin();
             i != toRemove.end(); ++i) {
            s.provisional.erase(i);
        }

#ifdef DEBUG_FEEDER
        std::cerr << "reap: have "
                  << s.satisfied.size() << " satisfied, "
                  << s.provisional.size() << " provisional, "
                  << s.completed.size() << " completed, reaped "
                  << reaped << std::endl;
#endif
    }
};

#endif
