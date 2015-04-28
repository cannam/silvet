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
#include <future>

#include <cstdio>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::future;
using std::async;
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
    desc.description = "Sets the tradeoff of processing speed against transcription quality. Draft mode is tuned in favour of overall speed; Live mode is tuned in favour of lower latency; while Intensive mode (the default) will almost always produce the best results.";
    desc.minValue = 0;
    desc.maxValue = 2;
    desc.defaultValue = int(defaultMode);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.push_back("Draft (faster)"); 
    desc.valueNames.push_back("Intensive (higher quality)");
    desc.valueNames.push_back("Live (lower latency)");
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
    d.description = "Overall note transcription. Each note has time, duration, estimated pitch, and a synthetic MIDI velocity (1-127) estimated from the strength of the pitch in the mixture.";
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
            d.binNames.push_back(noteName(i, 0, 1));
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
            d.binNames.push_back(chromaName(i));
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
Silvet::chromaName(int pitch) const
{
    static const char *names[] = {
        "A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"
    };

    return names[pitch];
}
    
std::string
Silvet::noteName(int note, int shift, int shiftCount) const
{
    string n = chromaName(note % 12);

    int oct = (note + 9) / 12; 
    
    char buf[30];

    float pshift = 0.f;
    if (shiftCount > 1) {
        // see noteFrequency below
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
Silvet::noteFrequency(int note, int shift, int shiftCount) const
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

    if (m_mode != HighQualityMode) {
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

    // For params.q, the MIREX code uses 0.8, but it seems that with
    // atomHopFactor of 0.3, using q == 0.9 or lower drops the FFT
    // size to 512 from 1024 and alters some other processing
    // parameters, making everything much, much slower. Could be a
    // flaw in the CQ parameter calculations, must check. For
    // atomHopFactor == 1, q == 0.8 is fine
    params.q = (m_mode == HighQualityMode ? 0.95 : 0.8);
    params.atomHopFactor = (m_mode == HighQualityMode ? 0.3 : 1.0);
    params.threshold = 0.0005;
    params.window = CQParameters::Hann;

    m_cq = new CQSpectrogram(params, CQSpectrogram::InterpolateLinear);

//    cerr << "CQ bins = " << m_cq->getTotalBins() << endl;
//    cerr << "CQ min freq = " << m_cq->getMinFrequency() << " (and for confirmation, freq of bin 0 = " << m_cq->getBinFrequency(0) << ")" << endl;
    
    m_colsPerSec = (m_mode == DraftMode ? 25 : 50);

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
        transcribe(cqout, fs);
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

void
Silvet::transcribe(const Grid &cqout, Silvet::FeatureSet &fs)
{
    Grid filtered = preProcess(cqout);

    if (filtered.empty()) return;
    
    const InstrumentPack &pack(getPack(m_instrument));

    for (int i = 0; i < (int)filtered.size(); ++i) {
        Feature f;
        for (int j = 0; j < pack.templateHeight; ++j) {
            f.values.push_back(float(filtered[i][j]));
        }
        fs[m_fcqOutputNo].push_back(f);
    }

    int width = filtered.size();

    Grid localPitches(width);

    bool wantShifts = (m_mode == HighQualityMode) && m_fineTuning;
    int shiftCount = 1;
    if (wantShifts) {
        shiftCount = pack.templateMaxShift * 2 + 1;
    }

    vector<vector<int> > localBestShifts;
    if (wantShifts) {
        localBestShifts = vector<vector<int> >(width);
    }

#ifndef MAX_EM_THREADS
#define MAX_EM_THREADS 8
#endif

#if (defined(MAX_EM_THREADS) && (MAX_EM_THREADS > 1))
    for (int i = 0; i < width; ) {
        typedef future<pair<vector<double>, vector<int>>> EMFuture;
        vector<EMFuture> results;
        for (int j = 0; j < MAX_EM_THREADS && i + j < width; ++j) {
            results.push_back
                (async(std::launch::async,
                       [&](int index) {
                           return applyEM(pack, filtered.at(index), wantShifts);
                       }, i + j));
        }
        for (int j = 0; j < MAX_EM_THREADS && i + j < width; ++j) {
            auto out = results[j].get();
            localPitches[i+j] = out.first;
            if (wantShifts) localBestShifts[i+j] = out.second;
        }
        i += MAX_EM_THREADS;
    }
#else
    for (int i = 0; i < width; ++i) {
        auto out = applyEM(pack, filtered.at(i), wantShifts);
        localPitches[i] = out.first;
        if (wantShifts) localBestShifts[i] = out.second;
    }
#endif
        
    for (int i = 0; i < width; ++i) {

        // This returns a filtered column, and pushes the
        // up-to-max-polyphony activation column to m_pianoRoll
        vector<double> filtered = postProcess
            (localPitches[i], localBestShifts[i], wantShifts);

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
        
        FeatureList noteFeatures = noteTrack(shiftCount);

        for (FeatureList::const_iterator fi = noteFeatures.begin();
             fi != noteFeatures.end(); ++fi) {
            fs[m_notesOutputNo].push_back(*fi);
        }
    }
}

pair<vector<double>, vector<int> >
Silvet::applyEM(const InstrumentPack &pack,
                const vector<double> &column,
                bool wantShifts)
{
    double columnThreshold = 1e-5;
    
    if (m_mode == LiveMode) {
        columnThreshold /= 20;
    }
    
    vector<double> pitches(pack.templateNoteCount, 0.0);
    vector<int> bestShifts;
    
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

    int shiftCount = 1;
    if (wantShifts) {
        shiftCount = pack.templateMaxShift * 2 + 1;
    }
    
    for (int j = 0; j < pack.templateNoteCount; ++j) {

        pitches[j] = pitchDist[j] * sum;

        int bestShift = 0;
        float bestShiftValue = 0.0;
        if (wantShifts) {
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
            // In draft and live mode the CQ is an octave shorter,
            // returning 540 bins or equivalent, so we instead pad
            // them with an additional 5 or equivalent zeros.
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
    
vector<double>
Silvet::postProcess(const vector<double> &pitches,
                    const vector<int> &bestShifts,
                    bool wantShifts)
{
    const InstrumentPack &pack(getPack(m_instrument));

    vector<double> filtered;

    for (int j = 0; j < pack.templateNoteCount; ++j) {
        m_postFilter[j]->push(pitches[j]);
        filtered.push_back(m_postFilter[j]->get());
    }

    if (m_mode == LiveMode) {
        // In live mode with only a 12-bpo CQ, we are very likely to
        // get clusters of two or three high scores at a time for
        // neighbouring semitones. Eliminate these by picking only the
        // peaks. This means we can't recognise actual semitone chords
        // if they ever appear, but it's not as if live mode is good
        // enough for that to be a big deal anyway.
        for (int j = 0; j < pack.templateNoteCount; ++j) {
            if (j > 0 && j + 1 < pack.templateNoteCount &&
                filtered[j] >= filtered[j-1] &&
                filtered[j] >= filtered[j+1]) {
            } else {
                filtered[j] = 0.0;
            }
        }
    }

    // Threshold for level and reduce number of candidate pitches

    typedef std::multimap<double, int> ValueIndexMap;

    ValueIndexMap strengths;

    for (int j = 0; j < pack.templateNoteCount; ++j) {
        double strength = filtered[j];
        if (strength < pack.levelThreshold) continue;
        strengths.insert(ValueIndexMap::value_type(strength, j));
    }

    ValueIndexMap::const_iterator si = strengths.end();

    map<int, double> active;
    map<int, int> activeShifts;

    while (int(active.size()) < pack.maxPolyphony && si != strengths.begin()) {

        --si;

        double strength = si->first;
        int j = si->second;

        active[j] = strength;

        if (wantShifts) {
            activeShifts[j] = bestShifts[j];
        }
    }

    m_pianoRoll.push_back(active);

    if (wantShifts) {
        m_pianoRollShifts.push_back(activeShifts);
    }

    return filtered;
}

Vamp::Plugin::FeatureList
Silvet::noteTrack(int shiftCount)
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
    int durationThreshold = floor(0.1 / columnDuration); // columns
    if (durationThreshold < 1) durationThreshold = 1;

    FeatureList noteFeatures;

    if (width < durationThreshold + 1) {
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

        while (m_pianoRoll[start].find(note) != m_pianoRoll[start].end()) {
            --start;
        }
        ++start;

        if ((end - start) < durationThreshold) {
            continue;
        }

        emitNote(start, end, note, shiftCount, noteFeatures);
    }

//    cerr << "returning " << noteFeatures.size() << " complete note(s) " << endl;

    return noteFeatures;
}

void
Silvet::emitNote(int start, int end, int note, int shiftCount,
                 FeatureList &noteFeatures)
{
    int partStart = start;
    int partShift = 0;
    int partVelocity = 0;

    int partThreshold = floor(0.05 * m_colsPerSec);

    for (int i = start; i != end; ++i) {
        
        double strength = m_pianoRoll[i][note];

        int shift = 0;

        if (shiftCount > 1) {

            shift = m_pianoRollShifts[i][note];

            if (i == partStart) {
                partShift = shift;
            }

            if (i > partStart + partThreshold && shift != partShift) {
                
//                cerr << "i = " << i << ", partStart = " << partStart << ", shift = " << shift << ", partShift = " << partShift << endl;

                // pitch has changed, emit an intermediate note
                noteFeatures.push_back(makeNoteFeature(partStart,
                                                       i,
                                                       note,
                                                       partShift,
                                                       shiftCount,
                                                       partVelocity));
                partStart = i;
                partShift = shift;
                partVelocity = 0;
            }
        }

        int v;
        if (m_mode == LiveMode) {
            v = round(strength * 20);
        } else {
            v = round(strength * 2);
        }
        if (v > partVelocity) {
            partVelocity = v;
        }
    }

    if (end >= partStart + partThreshold) {
        noteFeatures.push_back(makeNoteFeature(partStart,
                                               end,
                                               note,
                                               partShift,
                                               shiftCount,
                                               partVelocity));
    }
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
                        int shiftCount,
                        int velocity)
{
    Feature f;

    f.hasTimestamp = true;
    f.timestamp = getColumnTimestamp(start);

    f.hasDuration = true;
    f.duration = getColumnTimestamp(end) - f.timestamp;

    f.values.clear();

    f.values.push_back
        (noteFrequency(note, shift, shiftCount));

    float inputGain = getInputGainAt(f.timestamp);
//    cerr << "adjusting velocity from " << velocity << " to " << round(velocity/inputGain) << endl;
    velocity = round(velocity / inputGain);
    if (velocity > 127) velocity = 127;
    if (velocity < 1) velocity = 1;
    f.values.push_back(velocity);

    f.label = noteName(note, shift, shiftCount);

    return f;
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
                        
