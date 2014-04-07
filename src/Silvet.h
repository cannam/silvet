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

#ifndef SILVET_H
#define SILVET_H

#include <vamp-sdk/Plugin.h>

#include <vector>
#include <string>
#include <set>

#include "maths/MedianFilter.h"

using std::string;
using std::vector;
using std::set;

class Resampler;
class CQInterpolated;

class Silvet : public Vamp::Plugin
{
public:
    Silvet(float inputSampleRate);
    virtual ~Silvet();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    int getPluginVersion() const;
    string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    Resampler *m_resampler;
    CQInterpolated *m_cq;

    typedef vector<vector<double> > Grid;

    vector<MedianFilter<double> *> m_postFilter;
    vector<set<int> > m_pianoRoll;

    Grid preProcess(const Grid &);
    FeatureList postProcess(const vector<double> &);
    FeatureSet transcribe(const Grid &);

    string noteName(int n) const;
    float noteFrequency(int n) const;

    int m_blockSize;
    int m_columnCount;
    int m_reducedColumnCount;
    Vamp::RealTime m_startTime;

    mutable int m_notesOutputNo;
    mutable int m_cqOutputNo;
    mutable int m_fcqOutputNo;
    mutable int m_pitchOutputNo;
};

#endif
