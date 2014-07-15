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

Silvet::Silvet(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_instruments(InstrumentPack::listInstrumentPacks()),
    m_resampler(0),
    m_cq(0),
    m_hqMode(true),
    m_fineTuning(false),
    m_instrument(0),
    m_colsPerSec(50)
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
    return 1;
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
    desc.description = "Determines the tradeoff of processing speed against transcription quality";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 1;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.push_back("Draft (faster)"); 
    desc.valueNames.push_back("Intensive (higher quality)");
    list.push_back(desc);

    desc.identifier = "instrument";
    desc.name = "Instrument";
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
    } else if (identifier == "instrument") {
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

    d.identifier = "timefreq";
    d.name = "Time-frequency distribution";
    d.description = "Filtered constant-Q time-frequency distribution used as input to the expectation-maximisation algorithm";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = m_instruments[0].templateHeight;
    d.binNames.clear();
    if (m_cq) {
        char name[20];
        for (int i = 0; i < m_instruments[0].templateHeight; ++i) {
            // We have a 600-bin (10 oct 60-bin CQ) of which the
            // lowest-frequency 55 bins have been dropped, for a
            // 545-bin template. The native CQ bins go high->low
            // frequency though, so these are still the first 545 bins
            // as reported by getBinFrequency, though in reverse order
            float freq = m_cq->getBinFrequency
                (m_instruments[0].templateHeight - i - 1);
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

    return list;
}

std::string
Silvet::noteName(int note, int shift, int shiftCount) const
{
    static const char *names[] = {
        "A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"
    };

    const char *n = names[note % 12];

    int oct = (note + 9) / 12; 
    
    char buf[30];

    float pshift = 0.f;
    if (shiftCount > 1) {
        // see noteFrequency below
        pshift = 
            float((shiftCount - shift) - int(shiftCount / 2) - 1) / shiftCount;
    }

    if (pshift > 0.f) {
        sprintf(buf, "%s%d+%dc", n, oct, int(round(pshift * 100)));
    } else if (pshift < 0.f) {
        sprintf(buf, "%s%d-%dc", n, oct, int(round((-pshift) * 100)));
    } else {
        sprintf(buf, "%s%d", n, oct);
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

    return float(27.5 * pow(2.0, (note + pshift) / 12.0));
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

    double minFreq = 27.5;

    if (!m_hqMode) {
        // We don't actually return any notes from the bottom octave,
        // so we can just pad with zeros
        minFreq *= 2;
    }

    CQParameters params(processingSampleRate,
                        minFreq, 
                        processingSampleRate / 3,
                        processingBPO);

    params.q = 0.95; // MIREX code uses 0.8, but it seems 0.9 or lower
                     // drops the FFT size to 512 from 1024 and alters
                     // some other processing parameters, making
                     // everything much, much slower. Could be a flaw
                     // in the CQ parameter calculations, must check
    params.atomHopFactor = 0.3;
    params.threshold = 0.0005;
    params.window = CQParameters::Hann;

    m_cq = new CQSpectrogram(params, CQSpectrogram::InterpolateLinear);

    m_colsPerSec = m_hqMode ? 50 : 25;

    for (int i = 0; i < (int)m_postFilter.size(); ++i) {
        delete m_postFilter[i];
    }
    m_postFilter.clear();
    for (int i = 0; i < m_instruments[0].templateNoteCount; ++i) {
        m_postFilter.push_back(new MedianFilter<double>(3));
    }
    m_pianoRoll.clear();
    m_columnCount = 0;
    m_startTime = RealTime::zeroTime;
    m_signalMax = 0.0;
}

Silvet::FeatureSet
Silvet::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    if (m_columnCount == 0) {
        m_startTime = timestamp;
    }
    
    vector<double> data;
    for (int i = 0; i < m_blockSize; ++i) {
        double d = inputBuffers[0][i];
        if (fabs(d) > m_signalMax) {
            m_signalMax = fabs(d);
        } 
    }
    for (int i = 0; i < m_blockSize; ++i) {
        double d = inputBuffers[0][i];
        if (m_signalMax > 0.0) {
            data.push_back(d / m_signalMax * 0.5);
        } else {
            data.push_back(0.0);
        }
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
    
    const InstrumentPack &pack = m_instruments[m_instrument];

    for (int i = 0; i < (int)filtered.size(); ++i) {
        Feature f;
        for (int j = 0; j < pack.templateHeight; ++j) {
            f.values.push_back(float(filtered[i][j]));
        }
        fs[m_fcqOutputNo].push_back(f);
    }

    int width = filtered.size();

    int iterations = m_hqMode ? 20 : 10;

    //!!! pitches or notes? [terminology]
    Grid localPitches(width, vector<double>(pack.templateNoteCount, 0.0));

    bool wantShifts = m_hqMode && m_fineTuning;
    int shiftCount = 1;
    if (wantShifts) {
        shiftCount = pack.templateMaxShift * 2 + 1;
    }

    vector<vector<int> > localBestShifts;
    if (wantShifts) {
        localBestShifts = 
            vector<vector<int> >(width, vector<int>(pack.templateNoteCount, 0));
    }

    vector<bool> present(width, false);

#pragma omp parallel for
    for (int i = 0; i < width; ++i) {

        double sum = 0.0;
        for (int j = 0; j < pack.templateHeight; ++j) {
            sum += filtered.at(i).at(j);
        }
        if (sum < 1e-5) continue;

        present[i] = true;

        EM em(&pack, m_hqMode);

        em.setPitchSparsity(pack.pitchSparsity);
        em.setSourceSparsity(pack.sourceSparsity);

        for (int j = 0; j < iterations; ++j) {
            em.iterate(filtered.at(i).data());
        }

        const float *pitchDist = em.getPitchDistribution();
        const float *const *shiftDist = em.getShifts();

        for (int j = 0; j < pack.templateNoteCount; ++j) {

            localPitches[i][j] = pitchDist[j] * sum;

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
                localBestShifts[i][j] = bestShift;
            }                
        }
    }
        
    for (int i = 0; i < width; ++i) {

        if (!present[i]) {
            // silent column
            for (int j = 0; j < pack.templateNoteCount; ++j) {
                m_postFilter[j]->push(0.0);
            }
            m_pianoRoll.push_back(map<int, double>());
            if (wantShifts) {
                m_pianoRollShifts.push_back(map<int, int>());
            }
            continue;
        }

        postProcess(localPitches[i], localBestShifts[i], wantShifts);
        
        FeatureList noteFeatures = noteTrack(shiftCount);

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

    const InstrumentPack &pack = m_instruments[m_instrument];

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
            // lowest 55 of them.
            // 
            // In draft mode the CQ is an octave shorter, returning
            // 540 bins, so we instead pad them with an additional 5
            // zeros.
            // 
            // We also need to reverse the column as we go, since the
            // raw CQ has the high frequencies first and we need it
            // the other way around.

            if (m_hqMode) {
                for (int j = 0; j < pack.templateHeight; ++j) {
                    int ix = inCol.size() - j - 55;
                    outCol[j] = inCol[ix];
                }
            } else {
                for (int j = 0; j < 5; ++j) {
                    outCol[j] = 0.0;
                }
                for (int j = 5; j < pack.templateHeight; ++j) {
                    int ix = inCol.size() - j + 4;
                    outCol[j] = inCol[ix];
                }
            }

            vector<double> noiseLevel1 = 
                MedianFilter<double>::filter(40, outCol);
            for (int j = 0; j < pack.templateHeight; ++j) {
                noiseLevel1[j] = std::min(outCol[j], noiseLevel1[j]);
            }

            vector<double> noiseLevel2 = 
                MedianFilter<double>::filter(40, noiseLevel1);
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
                    const vector<int> &bestShifts,
                    bool wantShifts)
{
    const InstrumentPack &pack = m_instruments[m_instrument];

    vector<double> filtered;

    for (int j = 0; j < pack.templateNoteCount; ++j) {
        m_postFilter[j]->push(pitches[j]);
        filtered.push_back(m_postFilter[j]->get());
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

    Feature f;
    f.hasTimestamp = true;
    f.hasDuration = true;

    double columnDuration = 1.0 / m_colsPerSec;
    int postFilterLatency = int(m_postFilter[0]->getSize() / 2);
    int partThreshold = floor(0.05 / columnDuration);

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
                f.timestamp = RealTime::fromSeconds
                    (columnDuration * (partStart - postFilterLatency) + 0.02);
                f.duration = RealTime::fromSeconds
                    (columnDuration * (i - partStart));
                f.values.clear();
                f.values.push_back
                    (noteFrequency(note, partShift, shiftCount));
                f.values.push_back(partVelocity);
                f.label = noteName(note, partShift, shiftCount);
                noteFeatures.push_back(f);
                partStart = i;
                partShift = shift;
                partVelocity = 0;
            }
        }

        int v = strength * 2;
        if (v > 127) v = 127;

        if (v > partVelocity) {
            partVelocity = v;
        }
    }

    if (end >= partStart + partThreshold) {
        f.timestamp = RealTime::fromSeconds
            (columnDuration * (partStart - postFilterLatency) + 0.02);
        f.duration = RealTime::fromSeconds
            (columnDuration * (end - partStart));
        f.values.clear();
        f.values.push_back
            (noteFrequency(note, partShift, shiftCount));
        f.values.push_back(partVelocity);
        f.label = noteName(note, partShift, shiftCount);
        noteFeatures.push_back(f);
    }
}
