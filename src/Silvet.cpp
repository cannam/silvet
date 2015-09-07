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
#include "flattendynamics-ladspa.h"
#include "LiveInstruments.h"

#include <vector>

#include <cstdio>

#if (defined(MAX_EM_THREADS) && (MAX_EM_THREADS > 1))
#include <future>
using std::future;
using std::async;
#endif

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;

using Vamp::RealTime;

static int processingSampleRate = 44100;

static int binsPerSemitoneLive = 1;
static int binsPerSemitoneNormal = 5;

static int minInputSampleRate = 100;
static int maxInputSampleRate = 192000;

static const Silvet::ProcessingMode defaultMode = Silvet::HighQualityMode;

Silvet::Silvet(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_instruments(InstrumentPack::listInstrumentPacks()),
    m_liveInstruments(LiveAdapter::adaptAll(m_instruments)),
    m_resampler(0),
    m_flattener(0),
    m_cq(0),
    m_mode(defaultMode),
    m_fineTuning(false),
    m_instrument(0),
    m_colsPerSec(50),
    m_haveStartTime(false)
{
}

Silvet::~Silvet()
{
    delete m_resampler;
    delete m_flattener;
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
    return "Estimate the note onsets, pitches, and durations that make up a music recording.";
}

string
Silvet::getMaker() const
{
    return "Queen Mary, University of London";
}

int
Silvet::getPluginVersion() const
{
    return 3;
}

string
Silvet::getCopyright() const
{
    return "Method by Emmanouil Benetos and Simon Dixon; plugin by Chris Cannam and Emmanouil Benetos. GPL licence.";
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
    desc.description = "Sets the tradeoff of processing speed against transcription quality. Live mode is much faster and detects notes with relatively low latency; Intensive mode (the default) is slower but will almost always produce better results.";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = int(defaultMode);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.push_back("Live (faster and lower latency)");
    desc.valueNames.push_back("Intensive (higher quality)");
    list.push_back(desc);

    desc.identifier = "instrument";
    desc.name = "Instrument";
    desc.unit = "";
    desc.description = "The instrument or instruments known to be present in the recording. This affects the set of instrument templates used, as well as the expected level of polyphony in the output. Using a more limited set of instruments than the default will also make the plugin run faster.\nNote that this plugin cannot isolate instruments: you can't use this setting to request notes from only one instrument in a recording with several. Instead, use this as a hint to the plugin about which instruments are actually present.";
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
    desc.description = "Return pitch estimates at finer than semitone resolution. This works only in Intensive mode. Notes that appear to drift in pitch will be split up into shorter notes with individually finer pitches.";
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
        return (float)(int)m_mode;
    } else if (identifier == "finetune") {
        return m_fineTuning ? 1.f : 0.f;
    } else if (identifier == "instrument") {
        return m_instrument;
    }
    return 0;
}

void
Silvet::setParameter(string identifier, float value) 
{
    if (identifier == "mode") {
        m_mode = (ProcessingMode)(int)(value + 0.5);
    } else if (identifier == "finetune") {
        m_fineTuning = (value > 0.5);
    } else if (identifier == "instrument") {
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
    d.description = "Overall note transcription. Each note has time, duration, estimated fundamental frequency, and a synthetic MIDI velocity (1-127) estimated from the strength of the pitch in the mixture.";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 2;
    d.binNames.push_back("Frequency");
    d.binNames.push_back("Velocity");
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = processingSampleRate / (m_cq ? m_cq->getColumnHop() : 62);
    d.hasDuration = true;
    m_notesOutputNo = list.size();
    list.push_back(d);

    d.identifier = "onsets";
    d.name = "Note onsets";
    d.description = "Note onsets, without durations. These can be calculated sooner than complete notes, because it isn't necessary to wait for a note to finish before returning its feature. Each event has time, estimated fundamental frequency in Hz, and a synthetic MIDI velocity (1-127) estimated from the strength of the pitch in the mixture.";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 2;
    d.binNames.push_back("Frequency");
    d.binNames.push_back("Velocity");
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = processingSampleRate / (m_cq ? m_cq->getColumnHop() : 62);
    d.hasDuration = false;
    m_onsetsOutputNo = list.size();
    list.push_back(d);

    d.identifier = "onoffsets";
    d.name = "Note onsets and offsets";
    d.description = "Note onsets and offsets as separate events. Each onset event has time, estimated fundamental frequency in Hz, and a synthetic MIDI velocity (1-127) estimated from the strength of the pitch in the mixture. Offsets are represented in the same way but with a velocity of 0.";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 2;
    d.binNames.push_back("Frequency");
    d.binNames.push_back("Velocity");
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = processingSampleRate / (m_cq ? m_cq->getColumnHop() : 62);
    d.hasDuration = false;
    m_onOffsetsOutputNo = list.size();
    list.push_back(d);

    d.identifier = "timefreq";
    d.name = "Time-frequency distribution";
    d.description = "Filtered constant-Q time-frequency distribution as used as input to the expectation-maximisation algorithm.";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = getPack(0).templateHeight;
    d.binNames.clear();
    if (m_cq) {
        char name[50];
        for (int i = 0; i < getPack(0).templateHeight; ++i) {
            // We have a 600-bin (10 oct 60-bin CQ) of which the
            // lowest-frequency 55 bins have been dropped, for a
            // 545-bin template. The native CQ bins go high->low
            // frequency though, so these are still the first 545 bins
            // as reported by getBinFrequency, though in reverse order
            float freq = m_cq->getBinFrequency
                (getPack(0).templateHeight - i - 1);
            sprintf(name, "%.1f Hz", freq);
            d.binNames.push_back(name);
        }
    }
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = m_colsPerSec;
    d.hasDuration = false;
    m_fcqOutputNo = list.size();
    list.push_back(d);

    d.identifier = "pitchactivation";
    d.name = "Pitch activation distribution";
    d.description = "Pitch activation distribution resulting from expectation-maximisation algorithm, prior to note extraction.";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = getPack(0).templateNoteCount;
    d.binNames.clear();
    if (m_cq) {
        for (int i = 0; i < getPack(0).templateNoteCount; ++i) {
            d.binNames.push_back(getNoteName(i, 0));
        }
    }
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = m_colsPerSec;
    d.hasDuration = false;
    m_pitchOutputNo = list.size();
    list.push_back(d);

    d.identifier = "chroma";
    d.name = "Pitch chroma distribution";
    d.description = "Pitch chroma distribution formed by wrapping the un-thresholded pitch activation distribution into a single octave of semitone bins.";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 12;
    d.binNames.clear();
    if (m_cq) {
        for (int i = 0; i < 12; ++i) {
            d.binNames.push_back(getChromaName(i));
        }
    }
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = m_colsPerSec;
    d.hasDuration = false;
    m_chromaOutputNo = list.size();
    list.push_back(d);

    d.identifier = "templates";
    d.name = "Templates";
    d.description = "Constant-Q spectral templates for the selected instrument pack.";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = getPack(0).templateHeight;
    d.binNames.clear();
    if (m_cq) {
        char name[50];
        for (int i = 0; i < getPack(0).templateHeight; ++i) {
            // We have a 600-bin (10 oct 60-bin CQ) of which the
            // lowest-frequency 55 bins have been dropped, for a
            // 545-bin template. The native CQ bins go high->low
            // frequency though, so these are still the first 545 bins
            // as reported by getBinFrequency, though in reverse order
            float freq = m_cq->getBinFrequency
                (getPack(0).templateHeight - i - 1);
            sprintf(name, "%.1f Hz", freq);
            d.binNames.push_back(name);
        }
    }
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = m_colsPerSec;
    d.hasDuration = false;
    m_templateOutputNo = list.size();
    list.push_back(d);

    return list;
}

std::string
Silvet::getChromaName(int pitch) const
{
    static const char *names[] = {
        "A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"
    };

    return names[pitch];
}
    
std::string
Silvet::getNoteName(int note, int shift) const
{
    string n = getChromaName(note % 12);

    int oct = (note + 9) / 12; 
    
    char buf[30];

    float pshift = 0.f;
    int shiftCount = getShiftCount();
    if (shiftCount > 1) {
        // see getNoteFrequency below
        pshift = 
            float((shiftCount - shift) - int(shiftCount / 2) - 1) / shiftCount;
    }

    if (pshift > 0.f) {
        sprintf(buf, "%s%d+%dc", n.c_str(), oct, int(round(pshift * 100)));
    } else if (pshift < 0.f) {
        sprintf(buf, "%s%d-%dc", n.c_str(), oct, int(round((-pshift) * 100)));
    } else {
        sprintf(buf, "%s%d", n.c_str(), oct);
    }

    return buf;
}

float
Silvet::getNoteFrequency(int note, int shift) const
{
    // Convert shift number to a pitch shift. The given shift number
    // is an offset into the template array, which starts with some
    // zeros, followed by the template, then some trailing zeros.
    // 
    // Example: if we have templateMaxShift == 2 and thus shiftCount
    // == 5, then the number will be in the range 0-4 and the template
    // will have 2 zeros at either end. Thus number 2 represents the
    // template "as recorded", for a pitch shift of 0; smaller indices
    // represent moving the template *up* in pitch (by introducing
    // zeros at the start, which is the low-frequency end), for a
    // positive pitch shift; and higher values represent moving it
    // down in pitch, for a negative pitch shift.

    float pshift = 0.f;
    int shiftCount = getShiftCount();
    if (shiftCount > 1) {
        pshift = 
            float((shiftCount - shift) - int(shiftCount / 2) - 1) / shiftCount;
    }

    float freq = float(27.5 * pow(2.0, (note + pshift) / 12.0));

//    cerr << "note = " << note << ", shift = " << shift << ", shiftCount = "
//         << shiftCount << ", obtained freq = " << freq << endl;
    
    return freq;
}

bool
Silvet::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (m_inputSampleRate < minInputSampleRate ||
        m_inputSampleRate > maxInputSampleRate) {
	cerr << "Silvet::initialise: Unsupported input sample rate "
             << m_inputSampleRate << " (supported min " << minInputSampleRate
             << ", max " << maxInputSampleRate << ")" << endl;
        return false;
    }

    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) {
	cerr << "Silvet::initialise: Unsupported channel count " << channels
             << " (supported min " << getMinChannelCount() << ", max "
             << getMaxChannelCount() << ")" << endl;
        return false;
    }

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
    delete m_flattener;
    delete m_cq;

    if (m_inputSampleRate != processingSampleRate) {
	m_resampler = new Resampler(m_inputSampleRate, processingSampleRate);
    } else {
	m_resampler = 0;
    }

    m_flattener = new FlattenDynamics(m_inputSampleRate); // before resampling
    m_flattener->reset();

    // this happens to be processingSampleRate / 3, and is the top
    // freq used for the EM templates:
    double maxFreq = 14700;

    if (m_mode == LiveMode) {
        // We only have 12 bpo rather than 60, so we need the top bin
        // to be the middle one of the top 5, i.e. 2/5 of a semitone
        // lower than 14700
        maxFreq *= powf(2.0, -1.0 / 30.0);
    }
    
    double minFreq = 27.5;

    if (m_mode == LiveMode) {
        // We don't actually return any notes from the bottom octave,
        // so we can just pad with zeros
        minFreq *= 2;
    }

    int bpo = 12 *
        (m_mode == LiveMode ? binsPerSemitoneLive : binsPerSemitoneNormal);

    CQParameters params(processingSampleRate,
                        minFreq, 
                        maxFreq,
                        bpo);

    params.q = 0.8;
    params.atomHopFactor = (m_mode == LiveMode ? 1.0 : 0.3);
    params.threshold = 0.0005;
    params.decimator =
        (m_mode == LiveMode ?
         CQParameters::FasterDecimator : CQParameters::BetterDecimator);
    params.window = CQParameters::Hann;

    m_cq = new CQSpectrogram(params, CQSpectrogram::InterpolateLinear);

//    cerr << "CQ bins = " << m_cq->getTotalBins() << endl;
//    cerr << "CQ min freq = " << m_cq->getMinFrequency() << " (and for confirmation, freq of bin 0 = " << m_cq->getBinFrequency(0) << ")" << endl;
    
    m_colsPerSec = 50;

    for (int i = 0; i < (int)m_postFilter.size(); ++i) {
        delete m_postFilter[i];
    }
    m_postFilter.clear();
    int postFilterLength = 3;
    for (int i = 0; i < getPack(0).templateNoteCount; ++i) {
        m_postFilter.push_back(new MedianFilter<double>(postFilterLength));
    }
    m_pianoRoll.clear();
    m_inputGains.clear();
    m_columnCount = 0;
    m_resampledCount = 0;
    m_startTime = RealTime::zeroTime;
    m_haveStartTime = false;
}

Silvet::FeatureSet
Silvet::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    
    if (!m_haveStartTime) {

        m_startTime = timestamp;
        m_haveStartTime = true;

        insertTemplateFeatures(fs);
    }

    vector<float> flattened(m_blockSize);
    float gain = 1.f;
    m_flattener->connectInputPort
        (FlattenDynamics::AudioInputPort, inputBuffers[0]);
    m_flattener->connectOutputPort
        (FlattenDynamics::AudioOutputPort, &flattened[0]);
    m_flattener->connectOutputPort
        (FlattenDynamics::GainOutputPort, &gain);
    m_flattener->process(m_blockSize);

    m_inputGains[timestamp] = gain;
    
    vector<double> data;
    for (int i = 0; i < m_blockSize; ++i) {
        double d = flattened[i];
        data.push_back(d);
    }

    if (m_resampler) {

	data = m_resampler->process(data.data(), data.size());

        int hadCount = m_resampledCount;
        m_resampledCount += data.size();

        int resamplerLatency = m_resampler->getLatency();

        if (hadCount < resamplerLatency) {
            int stillToDrop = resamplerLatency - hadCount;
            if (stillToDrop >= int(data.size())) {
                return fs;
            } else {
                data = vector<double>(data.begin() + stillToDrop, data.end());
            }
        }
    }

    Grid cqout = m_cq->process(data);
    transcribe(cqout, fs);
    return fs;
}

Silvet::FeatureSet
Silvet::getRemainingFeatures()
{
    Grid cqout = m_cq->getRemainingOutput();
    FeatureSet fs;

    if (m_columnCount == 0) {
        // process() was never called, but we still want these
        insertTemplateFeatures(fs);
    } else {

        // Complete the transcription

        transcribe(cqout, fs);

        // And make sure any extant playing notes are finished and returned

        m_pianoRoll.push_back({});

        auto events = noteTrack();
    
        for (const auto &f : events.notes) {
            fs[m_notesOutputNo].push_back(f);
        }
    
        for (const auto &f : events.onsets) {
            fs[m_onsetsOutputNo].push_back(f);
        }
    
        for (const auto &f : events.onOffsets) {
            fs[m_onOffsetsOutputNo].push_back(f);
        }
    }

    return fs;
}

void
Silvet::insertTemplateFeatures(FeatureSet &fs)
{
    const InstrumentPack &pack = getPack(m_instrument);
    for (int i = 0; i < int(pack.templates.size()) * pack.templateNoteCount; ++i) {
        RealTime timestamp = RealTime::fromSeconds(double(i) / m_colsPerSec);
        Feature f;
        char buffer[50];
        sprintf(buffer, "Note %d", i + 1);
        f.label = buffer;
        f.hasTimestamp = true;
        f.timestamp = timestamp;
        f.values = pack.templates[i / pack.templateNoteCount]
            .data[i % pack.templateNoteCount];
        fs[m_templateOutputNo].push_back(f);
    }
}        

int
Silvet::getShiftCount() const
{
    bool wantShifts = (m_mode == HighQualityMode) && m_fineTuning;
    int shiftCount = 1;
    if (wantShifts) {
        const InstrumentPack &pack(getPack(m_instrument));
        shiftCount = pack.templateMaxShift * 2 + 1;
    }
    return shiftCount;
}

void
Silvet::transcribe(const Grid &cqout, Silvet::FeatureSet &fs)
{
    Grid filtered = preProcess(cqout);

    if (filtered.empty()) return;
    
    const InstrumentPack &pack(getPack(m_instrument));

    int width = filtered.size();

    double silenceThreshold = 0.01;
    
    for (int i = 0; i < width; ++i) {

        RealTime timestamp = getColumnTimestamp(m_pianoRoll.size() - 1 + i);
        float inputGain = getInputGainAt(timestamp);

        Feature f;
        double rms = 0.0;

        for (int j = 0; j < pack.templateHeight; ++j) {
            double v = filtered[i][j];
            rms += v * v;
            f.values.push_back(float(v));
        }

        rms = sqrt(rms / pack.templateHeight);
        if (rms / inputGain < silenceThreshold) {
            filtered[i].clear();
        }
        
        fs[m_fcqOutputNo].push_back(f);
    }
    
    Grid localPitches(width);

    int shiftCount = getShiftCount();
    bool wantShifts = (shiftCount > 1);

    vector<vector<int> > localBestShifts;
    if (wantShifts) {
        localBestShifts = vector<vector<int> >(width);
    }

    int emThreadCount = 1;

#if (defined(MAX_EM_THREADS) && (MAX_EM_THREADS > 1))
    emThreadCount = MAX_EM_THREADS;

    if (emThreadCount > int(std::thread::hardware_concurrency())) {
        emThreadCount = std::thread::hardware_concurrency();
    }
    if (m_mode == LiveMode && pack.templates.size() == 1) {
        // The EM step is probably not slow enough to merit it
        emThreadCount = 1;
    }

    if (emThreadCount > 1) {
        for (int i = 0; i < width; ) {
            typedef future<pair<vector<double>, vector<int>>> EMFuture;
            vector<EMFuture> results;
            for (int j = 0; j < emThreadCount && i + j < width; ++j) {
                const vector<double> &column = filtered.at(i + j);
                results.push_back
                    (async(std::launch::async,
                           [&]() { return applyEM(pack, column); }));
            }
            for (int j = 0; j < emThreadCount && i + j < width; ++j) {
                auto out = results[j].get();
                localPitches[i+j] = out.first;
                if (wantShifts) localBestShifts[i+j] = out.second;
            }
            i += emThreadCount;
        }
    }
#endif

    if (emThreadCount == 1) {
        for (int i = 0; i < width; ++i) {
            auto out = applyEM(pack, filtered.at(i));
            localPitches[i] = out.first;
            if (wantShifts) localBestShifts[i] = out.second;
        }
    }
        
    for (int i = 0; i < width; ++i) {

        vector<double> filtered;

        for (int j = 0; j < pack.templateNoteCount; ++j) {
            m_postFilter[j]->push(localPitches[i][j]);
            filtered.push_back(m_postFilter[j]->get());
        }

        RealTime timestamp = getColumnTimestamp(m_pianoRoll.size() - 1);
        float inputGain = getInputGainAt(timestamp);

        Feature f;
        for (int j = 0; j < (int)filtered.size(); ++j) {
            float v = filtered[j];
            if (v < pack.levelThreshold) v = 0.f;
            f.values.push_back(v / inputGain);
        }
        fs[m_pitchOutputNo].push_back(f);

        f.values.clear();
        f.values.resize(12);
        for (int j = 0; j < (int)filtered.size(); ++j) {
            f.values[j % 12] += filtered[j] / inputGain;
        }
        fs[m_chromaOutputNo].push_back(f);

        // This pushes the up-to-max-polyphony activation column to
        // m_pianoRoll
        postProcess(filtered, localBestShifts[i]);

        auto events = noteTrack();

        for (const auto &f : events.notes) {
            fs[m_notesOutputNo].push_back(f);
        }

        for (const auto &f : events.onsets) {
            fs[m_onsetsOutputNo].push_back(f);
        }

        for (const auto &f : events.onOffsets) {
            fs[m_onOffsetsOutputNo].push_back(f);
        }
    }
}

pair<vector<double>, vector<int> >
Silvet::applyEM(const InstrumentPack &pack,
                const vector<double> &column)
{
    double columnThreshold = 1e-5;
    
    if (m_mode == LiveMode) {
        columnThreshold /= 15;
    }
    
    vector<double> pitches(pack.templateNoteCount, 0.0);
    vector<int> bestShifts;

    if (column.empty()) return { pitches, bestShifts };
    
    double sum = 0.0;
    for (int j = 0; j < pack.templateHeight; ++j) {
        sum += column.at(j);
    }
    if (sum < columnThreshold) return { pitches, bestShifts };

    EM em(&pack, m_mode == HighQualityMode);

    em.setPitchSparsity(pack.pitchSparsity);
    em.setSourceSparsity(pack.sourceSparsity);

    int iterations = (m_mode == HighQualityMode ? 20 : 10);

    for (int j = 0; j < iterations; ++j) {
        em.iterate(column.data());
    }

    const float *pitchDist = em.getPitchDistribution();
    const float *const *shiftDist = em.getShifts();

    int shiftCount = getShiftCount();
    
    for (int j = 0; j < pack.templateNoteCount; ++j) {

        pitches[j] = pitchDist[j] * sum;

        int bestShift = 0;
        float bestShiftValue = 0.0;
        if (shiftCount > 1) {
            for (int k = 0; k < shiftCount; ++k) {
                float value = shiftDist[k][j];
                if (k == 0 || value > bestShiftValue) {
                    bestShiftValue = value;
                    bestShift = k;
                }
            }
            bestShifts.push_back(bestShift);
        }                
    }

    return { pitches, bestShifts };
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

    const InstrumentPack &pack(getPack(m_instrument));

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
            vector<double> outCol(pack.templateHeight);

            // In HQ mode, the CQ returns 600 bins and we ignore the
            // lowest 55 of them (assuming binsPerSemitone == 5).
            // 
            // In live mode the CQ is an octave shorter, returning 540
            // bins or equivalent, so we instead pad them with an
            // additional 5 or equivalent zeros.
            // 
            // We also need to reverse the column as we go, since the
            // raw CQ has the high frequencies first and we need it
            // the other way around.

            int bps = (m_mode == LiveMode ?
                       binsPerSemitoneLive : binsPerSemitoneNormal);
            
            if (m_mode == HighQualityMode) {
                for (int j = 0; j < pack.templateHeight; ++j) {
                    int ix = inCol.size() - j - (11 * bps);
                    outCol[j] = inCol[ix];
                }
            } else {
                for (int j = 0; j < bps; ++j) {
                    outCol[j] = 0.0;
                }
                for (int j = bps; j < pack.templateHeight; ++j) {
                    int ix = inCol.size() - j + (bps-1);
                    outCol[j] = inCol[ix];
                }
            }

            vector<double> noiseLevel1 = 
                MedianFilter<double>::filter(8 * bps, outCol);
            for (int j = 0; j < pack.templateHeight; ++j) {
                noiseLevel1[j] = std::min(outCol[j], noiseLevel1[j]);
            }

            vector<double> noiseLevel2 = 
                MedianFilter<double>::filter(8 * bps, noiseLevel1);
            for (int j = 0; j < pack.templateHeight; ++j) {
                outCol[j] = std::max(outCol[j] - noiseLevel2[j], 0.0);
            }

            out.push_back(outCol);
        }

        ++m_columnCount;
    }

    return out;
}
    
void
Silvet::postProcess(const vector<double> &pitches,
                    const vector<int> &bestShifts)
{
    const InstrumentPack &pack(getPack(m_instrument));

    // Threshold for level and reduce number of candidate pitches

    typedef std::multimap<double, int> ValueIndexMap;

    ValueIndexMap strengths;

    for (int j = 0; j < pack.templateNoteCount; ++j) {

        double strength = pitches[j];
        if (strength < pack.levelThreshold) continue;
        
        // In live mode with only a 12-bpo CQ, we are very likely to
        // get clusters of two or three high scores at a time for
        // neighbouring semitones. Eliminate these by picking only the
        // peaks (except that we never eliminate a note that has
        // already been established as currently playing). This means
        // we can't recognise actual semitone chords if they ever
        // appear, but it's not as if live mode is good enough for
        // that to be a big deal anyway.
        if (m_mode == LiveMode) {
            if (m_current.find(j) == m_current.end() &&
                (j == 0 ||
                 j + 1 == pack.templateNoteCount ||
                 pitches[j] < pitches[j-1] ||
                 pitches[j] < pitches[j+1])) {
                // not a peak or a currently-playing note: skip it
                continue;
            }
        }

        strengths.insert(ValueIndexMap::value_type(strength, j));
    }

    ValueIndexMap::const_iterator si = strengths.end();

    map<int, double> active;
    map<int, int> activeShifts;

    int shiftCount = getShiftCount();
    
    while (int(active.size()) < pack.maxPolyphony && si != strengths.begin()) {

        --si;

        double strength = si->first;
        int j = si->second;

        active[j] = strength;

        if (shiftCount > 1) {
            activeShifts[j] = bestShifts[j];
        }
    }

    m_pianoRoll.push_back(active);

    if (shiftCount > 1) {
        m_pianoRollShifts.push_back(activeShifts);
    }

    return;
}

Silvet::FeatureChunk
Silvet::noteTrack()
{        
    // Minimum duration pruning, and conversion to notes. We can only
    // report notes that have just ended (i.e. that are absent in the
    // latest active set but present in the prior set in the piano
    // roll) -- any notes that ended earlier will have been reported
    // already, and if they haven't ended, we don't know their
    // duration.

    int width = m_pianoRoll.size() - 1;

    const map<int, double> &active = m_pianoRoll[width];

    double columnDuration = 1.0 / m_colsPerSec;

    // only keep notes >= 100ms or thereabouts
    double durationThrSec = 0.1;
    int durationThreshold = floor(durationThrSec / columnDuration); // in cols
    if (durationThreshold < 1) durationThreshold = 1;

    FeatureList noteFeatures, onsetFeatures, onOffsetFeatures;

    if (width < durationThreshold + 1) {
        return { noteFeatures, onsetFeatures, onOffsetFeatures };
    }
    
    for (map<int, double>::const_iterator ni = m_pianoRoll[width-1].begin();
         ni != m_pianoRoll[width-1].end(); ++ni) {

        int note = ni->first;

        int end = width;
        int start = end-1;

        while (m_pianoRoll[start].find(note) != m_pianoRoll[start].end()) {
            --start;
        }
        ++start;

        int duration = end - start;
            
        if (duration < durationThreshold) {
            continue;
        }

        if (duration == durationThreshold) {
            m_current.insert(note);
            emitOnset(start, note, onsetFeatures);
            emitOnset(start, note, onOffsetFeatures);
        }            
        
        if (active.find(note) == active.end()) {
            // the note was playing but just ended
            m_current.erase(note);
            emitNoteAndOffset(start, end, note, noteFeatures, onOffsetFeatures);
        } else { // still playing
            // repeated note detection: if level is greater than this
            // multiple of its previous value, then we end the note and
            // restart it with the same pitch
            double restartFactor = 1.5;
            if (duration >= durationThreshold * 2 &&
                (active.find(note)->second >
                 restartFactor * m_pianoRoll[width-1][note])) {
                m_current.erase(note);
                emitNoteAndOffset(start, end-1, note, noteFeatures, onOffsetFeatures);
                // and remove this so that we start counting the new
                // note's duration from the current position
                m_pianoRoll[width-1].erase(note);
            }
        }
    }

//    cerr << "returning " << noteFeatures.size() << " complete note(s) " << endl;

    return { noteFeatures, onsetFeatures, onOffsetFeatures };
}

void
Silvet::emitNoteAndOffset(int start, int end, int note,
                          FeatureList &noteFeatures,
                          FeatureList &onOffsetFeatures)
{
    // Emit the complete note-event feature, and its offset. We have
    // already emitted the note onset when it started -- that process
    // is separated out in order to get a faster response during live
    // tracking. However, if the note shift changes within the note
    // (which can happen only if we have fine-tuning switched on), we
    // emit an offset and then a new onset with the new shift.
    
    int partStart = start;
    int partShift = 0;
    double partStrength = 0;

    // NB this *must* be less than durationThreshold above
    int partThreshold = floor(0.05 * m_colsPerSec);

    for (int i = start; i != end; ++i) {
        
        double strength = m_pianoRoll[i][note];

        int shift = 0;

        if (getShiftCount() > 1) {

            shift = m_pianoRollShifts[i][note];

            if (i == partStart) {
                partShift = shift;
            }

            if (i > partStart + partThreshold && shift != partShift) {

                // pitch has changed, emit an intermediate note
                noteFeatures.push_back(makeNoteFeature(partStart,
                                                       i,
                                                       note,
                                                       partShift,
                                                       partStrength));

                onOffsetFeatures.push_back(makeOffsetFeature(i,
                                                             note,
                                                             partShift));
                
                partStart = i;
                partShift = shift;

                onOffsetFeatures.push_back(makeOnsetFeature(i,
                                                            note,
                                                            partShift,
                                                            partStrength));

                partStrength = 0;
            }
        }

        if (strength > partStrength) {
            partStrength = strength;
        }
    }

    if (end >= partStart + partThreshold) {

        noteFeatures.push_back(makeNoteFeature(partStart,
                                               end,
                                               note,
                                               partShift,
                                               partStrength));

        onOffsetFeatures.push_back(makeOffsetFeature(end,
                                                     note,
                                                     partShift));

    } else if (partStart > start) {

        // we have emitted an onset for this, so must add an offset
        onOffsetFeatures.push_back(makeOffsetFeature(end,
                                                     note,
                                                     partShift));
    }
}

void
Silvet::emitOnset(int start, int note, FeatureList &onOffsetFeatures)
{
    int len = int(m_pianoRoll.size());

    double onsetStrength = 0;

    int shift = 0;
    if (getShiftCount() > 1) {
        shift = m_pianoRollShifts[start][note];
    }
    
    for (int i = start; i < len; ++i) {
        double strength = m_pianoRoll[i][note];
        if (strength > onsetStrength) {
            onsetStrength = strength;
        }
    }

    if (onsetStrength == 0) return;
    
    onOffsetFeatures.push_back(makeOnsetFeature(start,
                                                note,
                                                shift,
                                                onsetStrength));
}

RealTime
Silvet::getColumnTimestamp(int column)
{
    double columnDuration = 1.0 / m_colsPerSec;
    int postFilterLatency = int(m_postFilter[0]->getSize() / 2);

    return m_startTime + RealTime::fromSeconds
        (columnDuration * (column - postFilterLatency) + 0.02);
}

Silvet::Feature
Silvet::makeNoteFeature(int start,
                        int end,
                        int note,
                        int shift,
                        double strength)
{
    Feature f;

    f.hasTimestamp = true;
    f.timestamp = getColumnTimestamp(start);

    f.hasDuration = true;
    f.duration = getColumnTimestamp(end) - f.timestamp;

    f.values.clear();
    f.values.push_back(getNoteFrequency(note, shift));
    f.values.push_back(getVelocityFor(strength, start));

    f.label = getNoteName(note, shift);

    return f;
}

Silvet::Feature
Silvet::makeOnsetFeature(int start,
                         int note,
                         int shift,
                         double strength)
{
    Feature f;

    f.hasTimestamp = true;
    f.timestamp = getColumnTimestamp(start);

    f.hasDuration = false;

    f.values.clear();
    f.values.push_back(getNoteFrequency(note, shift));
    f.values.push_back(getVelocityFor(strength, start));

    f.label = getNoteName(note, shift);

    return f;
}

Silvet::Feature
Silvet::makeOffsetFeature(int col,
                          int note,
                          int shift)
{
    Feature f;
    
    f.hasTimestamp = true;
    f.timestamp = getColumnTimestamp(col);

    f.hasDuration = false;

    f.values.clear();
    f.values.push_back(getNoteFrequency(note, shift));
    f.values.push_back(0); // velocity 0 for offset

    f.label = getNoteName(note, shift) + " off";

    return f;
}

int
Silvet::getVelocityFor(double strength, int column)
{
    RealTime rt = getColumnTimestamp(column + 1);

    float inputGain = getInputGainAt(rt);

    double scale = 2.0;
    if (m_mode == LiveMode) scale = 20.0;
        
    double velocity = round((strength * scale) / inputGain);
    
    if (velocity > 127.0) velocity = 127.0;
    if (velocity < 1.0) velocity = 1.0; // assume surpassed 0 threshold already

    return int(velocity);
}

float
Silvet::getInputGainAt(RealTime t)
{
    map<RealTime, float>::const_iterator i = m_inputGains.lower_bound(t);

    if (i == m_inputGains.end()) {
        if (i != m_inputGains.begin()) {
            --i;
        } else {
            return 1.f; // no data
        }
    }

//    cerr << "gain at time " << t << " = " << i->second << endl;

    return i->second;
}
                        
