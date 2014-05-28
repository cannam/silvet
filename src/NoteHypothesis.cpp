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

#include "NoteHypothesis.h"
#include "AgentFeederPoly.h"

#include <cmath>
#include <cassert>

#include <map>
#include <algorithm>

using Vamp::RealTime;

using std::cerr;
using std::endl;

#define DEBUG_NOTE_HYPOTHESIS 1

NoteHypothesis::NoteHypothesis()
{
    m_state = New;
}

NoteHypothesis::~NoteHypothesis()
{
}

bool
NoteHypothesis::isWithinTolerance(Observation s) const
{
    if (m_pending.empty()) {
        return true;
    }

    // check we are within a relatively close tolerance of the last
    // candidate
    Observations::const_iterator i = m_pending.end();
    --i;
    Observation last = *i;
    double r = s.value / last.value;
    int cents = lrint(1200.0 * (log(r) / log(2.0)));
#ifdef DEBUG_NOTE_HYPOTHESIS
    cerr << "isWithinTolerance: this " << s.value << " is " << cents
              << " cents from prior " << last.value << endl;
#endif
    if (cents < -60 || cents > 60) return false;

    // and within a slightly bigger tolerance of the current mean
    double meanFreq = getMeanFrequency();
    r = s.value / meanFreq;
    cents = lrint(1200.0 * (log(r) / log(2.0)));
#ifdef DEBUG_NOTE_HYPOTHESIS
    cerr << "isWithinTolerance: this " << s.value << " is " << cents
              << " cents from mean " << meanFreq << endl;
#endif
    if (cents < -80 || cents > 80) return false;
    
    return true;
}

bool
NoteHypothesis::isOutOfDateFor(Observation s) const
{
    if (m_pending.empty()) return false;

    Observations::const_iterator i = m_pending.end();
    --i;
    Observation last = *i;

#ifdef DEBUG_NOTE_HYPOTHESIS
    cerr << "isOutOfDateFor: this " << s.time << " is "
              << (s.time - last.time) << " from last " << last.time
              << " (threshold " << RealTime::fromMilliseconds(40) << ")"
              << endl;
#endif

    return ((s.time - last.time) > RealTime::fromMilliseconds(40));
}

bool 
NoteHypothesis::isSatisfied() const
{
    if (m_pending.empty()) return false;
    
    double meanConfidence = 0.0;
    for (Observations::const_iterator i = m_pending.begin();
         i != m_pending.end(); ++i) {
        meanConfidence += i->confidence;
    }
    meanConfidence /= m_pending.size();

    //!!! surely this depends on the hop size?
    int lengthRequired = 100;
    if (meanConfidence > 0.0) {
        lengthRequired = int(2.0 / meanConfidence + 0.5);
    }
    //!!! 
    lengthRequired = lengthRequired / 2;
    if (lengthRequired < 1) lengthRequired = 1;

#ifdef DEBUG_NOTE_HYPOTHESIS
    cerr << "meanConfidence " << meanConfidence << ", lengthRequired " << lengthRequired << endl;
#endif

    return ((int)m_pending.size() > lengthRequired);
}

static void printState(NoteHypothesis::State s)
{
    switch (s) {
    case NoteHypothesis::New: cerr << "New"; break;
    case NoteHypothesis::Provisional: cerr << "Provisional"; break;
    case NoteHypothesis::Rejected: cerr << "Rejected"; break;
    case NoteHypothesis::Satisfied: cerr << "Satisfied"; break;
    case NoteHypothesis::Expired: cerr << "Expired"; break;
    }
}

bool
NoteHypothesis::accept(Observation s)
{
    bool accept = false;

#ifdef DEBUG_NOTE_HYPOTHESIS
    cerr << "NoteHypothesis[" << this << "]::accept (value " << s.value << ", time " << s.time << ", confidence " << s.confidence << "): state ";
    printState(m_state);
    cerr << "..." << endl;
#endif

    static double negligibleConfidence = 0.0001;

    if (s.confidence < negligibleConfidence) {
        // avoid piling up a lengthy sequence of estimates that are
        // all acceptable but are in total not enough to cause us to
        // be satisfied
        if (m_state == New) {
            m_state = Rejected;
        }
        return false;
    }

    switch (m_state) {

    case New:
        m_state = Provisional;
        accept = true;
        break;

    case Provisional:
        if (isOutOfDateFor(s)) {
            m_state = Rejected;
        } else if (isWithinTolerance(s)) {
            accept = true;
        }
        break;
        
    case Satisfied:
        if (isOutOfDateFor(s)) {
            m_state = Expired;
        } else if (isWithinTolerance(s)) {
            accept = true;
        }
        break;

    case Rejected:
        break;

    case Expired:
        break;
    }

    if (accept) {
#ifdef DEBUG_NOTE_HYPOTHESIS
        cerr << "... accepting" << endl;
#endif
        m_pending.insert(s);
        if (m_state == Provisional && isSatisfied()) {
            m_state = Satisfied;
        }
    } else {
#ifdef DEBUG_NOTE_HYPOTHESIS
        cerr << "... not accepting" << endl;
#endif
    }

#ifdef DEBUG_NOTE_HYPOTHESIS
    cerr << "... -> ";
    printState(m_state);
    cerr << " (pending: " << m_pending.size() << ")" << endl;
#endif

    return accept;
}        

NoteHypothesis::State
NoteHypothesis::getState() const
{
    return m_state;
}

NoteHypothesis::Observations
NoteHypothesis::getAcceptedObservations() const
{
    if (m_state == Satisfied || m_state == Expired) {
        return m_pending;
    } else {
        return Observations();
    }
}

double
NoteHypothesis::getMedianFrequency() const
{
    if (m_pending.empty()) return 0.0;
    std::vector<double> freqs;
    for (Observations::const_iterator i = m_pending.begin();
         i != m_pending.end(); ++i) {
        freqs.push_back(i->value);
    }
    std::sort(freqs.begin(), freqs.end());
    return freqs[freqs.size()/2];
}

double
NoteHypothesis::getMeanFrequency() const
{
    double acc = 0.0;
    if (m_pending.empty()) return acc;
    for (Observations::const_iterator i = m_pending.begin();
         i != m_pending.end(); ++i) {
        acc += i->value;
    }
    acc /= m_pending.size();
    return acc;
}

double
NoteHypothesis::getMedianConfidence() const
{
    if (m_pending.empty()) return 0.0;
    std::vector<double> confs;
    for (Observations::const_iterator i = m_pending.begin();
         i != m_pending.end(); ++i) {
        confs.push_back(i->confidence);
    }
    std::sort(confs.begin(), confs.end());
    return confs[confs.size()/2];
}

NoteHypothesis::Note
NoteHypothesis::getAveragedNote() const
{
    Note n;

    n.time = getStartTime();
    n.duration = getDuration();
    n.freq = getMedianFrequency();
    n.confidence = getMedianConfidence();
    
    return n;
}

RealTime
NoteHypothesis::getStartTime() const
{
    if (!(m_state == Satisfied || m_state == Expired)) {
        return RealTime::zeroTime;
    } else {
        return m_pending.begin()->time;
    }
}

RealTime
NoteHypothesis::getDuration() const
{
//!!! test this! it is wrong
    if (!(m_state == Satisfied || m_state == Expired)) {
        return RealTime::zeroTime;
    } else {
        RealTime start = m_pending.begin()->time;
        Observations::const_iterator i = m_pending.end();
        --i;
        return i->time - start;
    }
}
