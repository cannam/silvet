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

#include "Silvet.h"
#include "EM.h"

#include <cq/CQSpectrogram.h>

#include "MedianFilter.h"
#include "constant-q-cpp/src/dsp/Resampler.h"

#include <vector>

#include <cstdio>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using Vamp::RealTime;

static int processingSampleRate = 44100;
static int processingBPO = 60;
static int processingHeight = 545;
static int processingNotes = 88;

Silvet::Silvet(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_instruments(InstrumentPack::listInstrumentPacks()),
    m_resampler(0),
    m_cq(0),
    m_hqMode(true),
    m_fineTuning(false),
    m_instrument(0)
{
}

Silvet::~Silvet()
{
    delete m_resampler;
    delete m_cq;
    for (int i = 0; i < (int)m_postFilter.size(); ++i) {
        delete m_postFilter[i];
    }
}

string
Silvet::getIdentifier() const
{
    return "silvet";
}

string
Silvet::getName() const
{
    return "Silvet Note Transcription";
}

string
Silvet::getDescription() const
{
    // Return something helpful here!
    return "";
}

string
Silvet::getMaker() const
{
    // Your name here
    return "";
}

int
Silvet::getPluginVersion() const
{
    return 1;
}

string
Silvet::getCopyright() const
{
    // This function is not ideally named.  It does not necessarily
    // need to say who made the plugin -- getMaker does that -- but it
    // should indicate the terms under which it is distributed.  For
    // example, "Copyright (year). All Rights Reserved", or "GPL"
    return "";
}

Silvet::InputDomain
Silvet::getInputDomain() const
{
    return TimeDomain;
}

size_t
Silvet::getPreferredBlockSize() const
{
    return 0;
}

size_t 
Silvet::getPreferredStepSize() const
{
    return 0;
}

size_t
Silvet::getMinChannelCount() const
{
    return 1;
}

size_t
Silvet::getMaxChannelCount() const
{
    return 1;
}

Silvet::ParameterList
Silvet::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor desc;
    desc.identifier = "mode";
    desc.name = "Processing mode";
    desc.unit = "";
    desc.description = "Determines the tradeoff of processing speed against transcription quality";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 1;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.push_back("Draft (faster)"); 
    desc.valueNames.push_back("Intensive (higher quality)");
    list.push_back(desc);

    desc.identifier = "soloinstrument";
    desc.name = "Solo instrument";
    desc.unit = "";
    desc.description = "The instrument known to be present in the recording, if there is only one";
    desc.minValue = 0;
    desc.maxValue = m_instruments.size()-1;
    desc.defaultValue = 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    for (int i = 0; i < int(m_instruments.size()); ++i) {
        desc.valueNames.push_back(m_instruments[i].name);
    }
    list.push_back(desc);

    desc.identifier = "finetune";
    desc.name = "Return fine pitch estimates";
    desc.unit = "";
    desc.description = "Return pitch estimates at finer than semitone resolution (works only in Intensive mode)";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    list.push_back(desc);

    return list;
}

float
Silvet::getParameter(string identifier) const
{
    if (identifier == "mode") {
        return m_hqMode ? 1.f : 0.f;
    } else if (identifier == "finetune") {
        return m_fineTuning ? 1.f : 0.f;
    } else if (identifier == "soloinstrument") {
        return m_instrument;
    }
    return 0;
}

void
Silvet::setParameter(string identifier, float value) 
{
    if (identifier == "mode") {
        m_hqMode = (value > 0.5);
    } else if (identifier == "finetune") {
        m_fineTuning = (value > 0.5);
    } else if (identifier == "soloinstrument") {
        m_instrument = lrintf(value);
    }
}

Silvet::ProgramList
Silvet::getPrograms() const
{
    ProgramList list;
    return list;
}

string
Silvet::getCurrentProgram() const
{
    return ""; 
}

void
Silvet::selectProgram(string name)
{
}

Silvet::OutputList
Silvet::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "notes";
    d.name = "Note transcription";
    d.description = "Overall note transcription across selected instruments";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 2;
    d.binNames.push_back("Frequency");
    d.binNames.push_back("Velocity");
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = m_inputSampleRate / (m_cq ? m_cq->getColumnHop() : 62);
    d.hasDuration = true;
    m_notesOutputNo = list.size();
    list.push_back(d);

    return list;
}

std::string
Silvet::noteName(int i) const
{
    static const char *names[] = {
        "A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"
    };

    const char *n = names[i % 12];

    int oct = (i + 9) / 12; 
    
    char buf[20];
    sprintf(buf, "%s%d", n, oct);

    return buf;
}

float
Silvet::noteFrequency(int note, int shiftCount) const
{
    float fineNote = float(note) / float(shiftCount);
    return float(27.5 * pow(2.0, fineNote / 12.0));
}

bool
Silvet::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    if (stepSize != blockSize) {
	cerr << "Silvet::initialise: Step size must be the same as block size ("
	     << stepSize << " != " << blockSize << ")" << endl;
	return false;
    }

    m_blockSize = blockSize;

    reset();

    return true;
}

void
Silvet::reset()
{
    delete m_resampler;
    delete m_cq;

    if (m_inputSampleRate != processingSampleRate) {
	m_resampler = new Resampler(m_inputSampleRate, processingSampleRate);
    } else {
	m_resampler = 0;
    }

    CQParameters params(processingSampleRate,
                        27.5, 
                        processingSampleRate / 3,
                        processingBPO);

    params.q = 0.95; // MIREX code uses 0.8, but it seems 0.9 or lower
                     // drops the FFT size to 512 from 1024 and alters
                     // some other processing parameters, making
                     // everything much, much slower. Could be a flaw
                     // in the CQ parameter calculations, must check
    params.atomHopFactor = 0.3;
    params.threshold = 0.0005;
    params.window = CQParameters::Hann; //!!! todo: test whether it makes any difference

    m_cq = new CQSpectrogram(params, CQSpectrogram::InterpolateLinear);

    m_colsPerSec = m_hqMode ? 50 : 25;

    for (int i = 0; i < (int)m_postFilter.size(); ++i) {
        delete m_postFilter[i];
    }
    m_postFilter.clear();
    for (int i = 0; i < processingNotes; ++i) {
        m_postFilter.push_back(new MedianFilter<double>(3));
    }
    m_pianoRoll.clear();
    m_columnCount = 0;
    m_startTime = RealTime::zeroTime;
}

Silvet::FeatureSet
Silvet::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    if (m_columnCount == 0) {
        m_startTime = timestamp;
    }
    
    vector<double> data;
    for (int i = 0; i < m_blockSize; ++i) {
        data.push_back(inputBuffers[0][i]);
    }

    if (m_resampler) {
	data = m_resampler->process(data.data(), data.size());
    }

    Grid cqout = m_cq->process(data);
    FeatureSet fs = transcribe(cqout);
    return fs;
}

Silvet::FeatureSet
Silvet::getRemainingFeatures()
{
    Grid cqout = m_cq->getRemainingOutput();
    FeatureSet fs = transcribe(cqout);
    return fs;
}

Silvet::FeatureSet
Silvet::transcribe(const Grid &cqout)
{
    Grid filtered = preProcess(cqout);

    FeatureSet fs;

    if (filtered.empty()) return fs;

    int width = filtered.size();

    int iterations = m_hqMode ? 20 : 10;

    vector<EM *> em(width, (EM *)0);
    vector<double> sums(width, 0.0);

#pragma omp parallel for
    for (int i = 0; i < width; ++i) {

        for (int j = 0; j < processingHeight; ++j) {
            sums[i] += filtered.at(i).at(j);
        }

        if (sums[i] < 1e-5) continue;

        em[i] = new EM(&m_instruments[m_instrument], m_hqMode);

        for (int j = 0; j < iterations; ++j) {
            em[i]->iterate(filtered.at(i).data());
        }
    }

    int shiftCount = 1;
    if (m_hqMode && m_fineTuning) {
        shiftCount = m_instruments[m_instrument].templateMaxShift * 2 + 1;
    }
        
    for (int i = 0; i < width; ++i) {

        if (!em[i]) {
            noteTrack(map<int, double>(), shiftCount);
            continue;
        }

        map<int, double> active = postProcess(em[i]->getPitchDistribution(),
                                              em[i]->getShifts(),
                                              shiftCount,
                                              sums[i]);
        
        delete em[i];
        
        FeatureList noteFeatures = noteTrack(active, shiftCount);

        for (FeatureList::const_iterator fi = noteFeatures.begin();
             fi != noteFeatures.end(); ++fi) {
            fs[m_notesOutputNo].push_back(*fi);
        }
    }

    return fs;
}

Silvet::Grid
Silvet::preProcess(const Grid &in)
{
    int width = in.size();

    int spacing = processingSampleRate / m_colsPerSec;

    // need to be careful that col spacing is an integer number of samples!
    assert(spacing * m_colsPerSec == processingSampleRate);

    Grid out;

    // We count the CQ latency in terms of processing hops, but
    // actually it probably isn't an exact number of hops so this
    // isn't quite accurate. But the small constant offset is
    // practically irrelevant compared to the jitter from the frame
    // size we reduce to in a moment
    int latentColumns = m_cq->getLatency() / m_cq->getColumnHop();

    for (int i = 0; i < width; ++i) {

        if (m_columnCount < latentColumns) {
            ++m_columnCount;
            continue;
        }

        int prevSampleNo = (m_columnCount - 1) * m_cq->getColumnHop();
        int sampleNo = m_columnCount * m_cq->getColumnHop();

        bool select = (sampleNo / spacing != prevSampleNo / spacing);

        if (select) {
            vector<double> inCol = in[i];
            vector<double> outCol(processingHeight);

            // we reverse the column as we go (the CQ output is
            // "upside-down", with high frequencies at the start of
            // each column, and we want it the other way around) and
            // then ignore the first 55 (lowest-frequency) bins,
            // giving us 545 bins instead of 600

            for (int j = 0; j < processingHeight; ++j) {
                int ix = inCol.size() - j - 55;
                outCol[j] = inCol[ix];
            }

            vector<double> noiseLevel1 = 
                MedianFilter<double>::filter(40, outCol);
            for (int j = 0; j < processingHeight; ++j) {
                noiseLevel1[j] = std::min(outCol[j], noiseLevel1[j]);
            }

            vector<double> noiseLevel2 = 
                MedianFilter<double>::filter(40, noiseLevel1);
            for (int j = 0; j < processingHeight; ++j) {
                outCol[j] = std::max(outCol[j] - noiseLevel2[j], 0.0);
            }

            out.push_back(outCol);
        }

        ++m_columnCount;
    }

    return out;
}
    
map<int, double>
Silvet::postProcess(const float *pitches,
                    const float *const *shifts,
                    int shiftCount,
                    double gain)
{
    vector<double> filtered;

    for (int j = 0; j < processingNotes; ++j) {
        m_postFilter[j]->push(pitches[j] * gain);
        filtered.push_back(m_postFilter[j]->get());
    }

    // Threshold for level and reduce number of candidate pitches

    int polyphony = 5;

    //!!! make this a parameter (was 4.8, try adjusting, compare levels against matlab code)
    double threshold = 6;
//    double threshold = 4.8;

    typedef std::multimap<double, int> ValueIndexMap;

    ValueIndexMap strengths;

    for (int j = 0; j < processingNotes; ++j) {

        double strength = filtered[j];
        if (strength < threshold) continue;

        // convert note number j to a pitch value p. If we are not
        // using fine tuning, p is the same as j; otherwise p is j *
        // shiftCount + preferred shift

        int p = j;

        if (m_hqMode && m_fineTuning && shiftCount > 1) {
            float bestShiftValue = 0.f;
            int bestShift = 0;
            for (int f = 0; f < shiftCount; ++f) {
                if (f == 0 || shifts[f][j] > bestShiftValue) {
                    bestShiftValue = shifts[f][j];
                    bestShift = f;
                }
            }
            //!!! I think our shift array per note is actually upside down, check this
            p = j * shiftCount + bestShift;
        }

        strengths.insert(ValueIndexMap::value_type(strength, p));
    }

    map<int, double> active;
    ValueIndexMap::const_iterator si = strengths.end();
    while (int(active.size()) < polyphony && si != strengths.begin()) {
        --si;
//        cerr << si->second << " : " << si->first << endl;
        active[si->second] = si->first;
        if (si == strengths.begin()) break;
    }

    return active;
}

Vamp::Plugin::FeatureList
Silvet::noteTrack(const map<int, double> &active, int shiftCount)
{        
    // Minimum duration pruning, and conversion to notes. We can only
    // report notes that have just ended (i.e. that are absent in the
    // latest active set but present in the last set in the piano
    // roll) -- any notes that ended earlier will have been reported
    // already, and if they haven't ended, we don't know their
    // duration.

    int postFilterLatency = int(m_postFilter[0]->getSize() / 2);

    int width = m_pianoRoll.size();

    double columnDuration = 1.0 / m_colsPerSec;

    // only keep notes >= 100ms or thereabouts
    int durationThreshold = floor(0.1 / columnDuration); // columns
    if (durationThreshold < 1) durationThreshold = 1;

    FeatureList noteFeatures;

    if (width < durationThreshold + 1) {
        m_pianoRoll.push_back(active);
        return noteFeatures;
    }
    
    //!!! try: repeated note detection? (look for change in first derivative of the pitch matrix)

    for (map<int, double>::const_iterator ni = m_pianoRoll[width-1].begin();
         ni != m_pianoRoll[width-1].end(); ++ni) {

        int note = ni->first;
        
        if (active.find(note) != active.end()) {
            // the note is still playing
            continue;
        }

        // the note was playing but just ended
        int end = width;
        int start = end-1;

        double maxStrength = 0.0;

        while (m_pianoRoll[start].find(note) != m_pianoRoll[start].end()) {
            double strength = m_pianoRoll[start][note];
            if (strength > maxStrength) {
                maxStrength = strength;
            }
            --start;
        }
        ++start;

        int duration = width - start;
//        cerr << "duration " << duration << " for just-ended note " << note << endl;
        if (duration < durationThreshold) {
            // spurious
            continue;
        }

        int velocity = maxStrength * 2;
        if (velocity > 127) velocity = 127;

//        cerr << "Found a genuine note, starting at " << columnDuration * start << " with duration " << columnDuration * duration << endl;

        Feature nf;
        nf.hasTimestamp = true;
        nf.timestamp = RealTime::fromSeconds
            (columnDuration * (start - postFilterLatency) + 0.02);
        nf.hasDuration = true;
        nf.duration = RealTime::fromSeconds
            (columnDuration * duration);
        nf.values.push_back(noteFrequency(note, shiftCount));
        nf.values.push_back(velocity);
        nf.label = noteName(note);
        noteFeatures.push_back(nf);
    }

    m_pianoRoll.push_back(active);

//    cerr << "returning " << noteFeatures.size() << " complete note(s) " << endl;

    return noteFeatures;
}

