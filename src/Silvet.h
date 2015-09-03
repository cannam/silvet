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

#include "MedianFilter.h"
#include "Instruments.h"

#ifndef MAX_EM_THREADS
#define MAX_EM_THREADS 8
#endif

using std::string;
using std::vector;
using std::set;
using std::map;

class Resampler;
class CQSpectrogram;
class FlattenDynamics;

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

    enum ProcessingMode {
        LiveMode = 0,
        HighQualityMode = 1,
    };

protected:
    const std::vector<InstrumentPack> m_instruments;
    const std::vector<InstrumentPack> m_liveInstruments;

    const InstrumentPack &getPack(int instrument) const {
        if (m_mode == LiveMode) {
            return m_liveInstruments[instrument];
        } else {
            return m_instruments[instrument];
        }
    }
    
    Resampler *m_resampler;
    FlattenDynamics *m_flattener;
    CQSpectrogram *m_cq;

    ProcessingMode m_mode;
    
    bool m_fineTuning;
    int m_instrument;
    int m_colsPerSec;

    typedef vector<vector<double> > Grid;

    vector<MedianFilter<double> *> m_postFilter;
    vector<map<int, double> > m_pianoRoll;
    vector<map<int, int> > m_pianoRollShifts;
    map<Vamp::RealTime, float> m_inputGains;
    set<int> m_current;

    Grid preProcess(const Grid &);

    std::pair<vector<double>, vector<int> > applyEM(const InstrumentPack &pack,
                                                    const vector<double> &column);
    
    void postProcess(const vector<double> &pitches,
                     const vector<int> &bestShifts); // -> piano roll column

    struct FeatureChunk {
        FeatureList notes;
        FeatureList onsets;
        FeatureList onOffsets;
    };

    int getShiftCount() const;
    
    FeatureChunk noteTrack(); // notes, on/offsets

    void emitNoteAndOffset(int start, int end, int note,
                           FeatureList &noteFeatures,
                           FeatureList &onOffsetFeatures);

    void emitOnset(int start, int note,
                  FeatureList &onOffsetFeatures);

    Vamp::RealTime getColumnTimestamp(int column);
    
    Feature makeNoteFeature(int start, int end, int note, int shift, double strength);
    Feature makeOnsetFeature(int start, int note, int shift, double strength);
    Feature makeOffsetFeature(int col, int note, int shift);

    int getVelocityFor(double strength, int column);
    
    float getInputGainAt(Vamp::RealTime t);

    void insertTemplateFeatures(FeatureSet &);
    
    void transcribe(const Grid &, FeatureSet &);

    string getChromaName(int n) const;
    string getNoteName(int n, int shift) const;
    float getNoteFrequency(int n, int shift) const;

    int m_blockSize;
    int m_columnCount;
    int m_resampledCount;
    Vamp::RealTime m_startTime;
    bool m_haveStartTime;

    mutable int m_notesOutputNo;
    mutable int m_onsetsOutputNo;
    mutable int m_onOffsetsOutputNo;
    mutable int m_fcqOutputNo;
    mutable int m_pitchOutputNo;
    mutable int m_templateOutputNo;
    mutable int m_chromaOutputNo;
};

#endif
