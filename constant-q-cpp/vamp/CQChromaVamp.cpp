/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/*
    Constant-Q library
    Copyright (c) 2013-2014 Queen Mary, University of London

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of the Centre for
    Digital Music; Queen Mary, University of London; and Chris Cannam
    shall not be used in advertising or otherwise to promote the sale,
    use or other dealings in this Software without prior written
    authorization.
*/

#include "CQChromaVamp.h"

#include "cq/Chromagram.h"

#include <algorithm>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

static const int defaultLowestOctave = 0;
static const int defaultOctaveCount = 7;
static const int defaultBPO = 36;
static const float defaultTuningFrequency = 440.f;

CQChromaVamp::CQChromaVamp(float inputSampleRate) :
    Vamp::Plugin(inputSampleRate),
    m_lowestOctave(defaultLowestOctave),
    m_octaveCount(defaultOctaveCount),
    m_tuningFrequency(defaultTuningFrequency),
    m_bpo(defaultBPO),
    m_chroma(0),
    m_haveStartTime(false),
    m_columnCount(0)
{
}

CQChromaVamp::~CQChromaVamp()
{
    delete m_chroma;
}

string
CQChromaVamp::getIdentifier() const
{
    return "cqchromavamp";
}

string
CQChromaVamp::getName() const
{
    return "CQ Chromagram";
}

string
CQChromaVamp::getDescription() const
{
    return "Extract a Constant-Q spectrogram with constant ratio of centre frequency to resolution from the audio, then wrap it around into a single-octave chromagram.";
}

string
CQChromaVamp::getMaker() const
{
    return "Queen Mary, University of London";
}

int
CQChromaVamp::getPluginVersion() const
{
    return 2;
}

string
CQChromaVamp::getCopyright() const
{
    return "Plugin by Chris Cannam. Method by Christian Schörkhuber and Anssi Klapuri. Copyright (c) 2015 QMUL. BSD/MIT licence.";
}

CQChromaVamp::ParameterList
CQChromaVamp::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor desc;

    desc.identifier = "lowestoct";
    desc.name = "Lowest Contributing Octave";
    desc.unit = "";
    desc.description = "Octave number of the lowest octave to include in the chromagram. Octave numbering is ASA standard, with -1 as the first octave in the MIDI range and middle-C being C4. The octave starts at C.";
    desc.minValue = -1;
    desc.maxValue = 12;
    desc.defaultValue = defaultLowestOctave;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "octaves";
    desc.name = "Contributing Octave Count";
    desc.unit = "octaves";
    desc.description = "Number of octaves to use when generating the Constant-Q transform. All octaves are wrapped around and summed to produce a single octave chromagram as output.";
    desc.minValue = 1;
    desc.maxValue = 12;
    desc.defaultValue = defaultOctaveCount;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "tuning";
    desc.name = "Tuning Frequency";
    desc.unit = "Hz";
    desc.description = "Frequency of concert A";
    desc.minValue = 360;
    desc.maxValue = 500;
    desc.defaultValue = defaultTuningFrequency;
    desc.isQuantized = false;
    list.push_back(desc);
    
    desc.identifier = "bpo";
    desc.name = "Bins per Octave";
    desc.unit = "bins";
    desc.description = "Number of constant-Q transform bins per octave";
    desc.minValue = 2;
    desc.maxValue = 480;
    desc.defaultValue = defaultBPO;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    return list;
}

float
CQChromaVamp::getParameter(std::string param) const
{
    if (param == "lowestoct") {
        return m_lowestOctave;
    }
    if (param == "octaves") {
        return m_octaveCount;
    }
    if (param == "tuning") {
        return m_tuningFrequency;
    }
    if (param == "bpo") {
        return m_bpo;
    }
    std::cerr << "WARNING: CQChromaVamp::getParameter: unknown parameter \""
              << param << "\"" << std::endl;
    return 0.0;
}

void
CQChromaVamp::setParameter(std::string param, float value)
{
    if (param == "lowestoct") {
        m_lowestOctave = int(value + 0.5f);
    } else if (param == "octaves") {
        m_octaveCount = int(value + 0.5f);
    } else if (param == "tuning") {
        m_tuningFrequency = value;
    } else if (param == "bpo") {
        m_bpo = int(value + 0.5f);
    } else {
        std::cerr << "WARNING: CQChromaVamp::setParameter: unknown parameter \""
                  << param << "\"" << std::endl;
    }
}

bool
CQChromaVamp::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (m_chroma) {
	delete m_chroma;
        m_chroma = 0;
    }

    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_stepSize = stepSize;
    m_blockSize = blockSize;

    reset();

    if (!m_chroma || !m_chroma->isValid()) {
        cerr << "CQVamp::initialise: Constant-Q parameters not valid! Not initialising" << endl;
        return false;
    }

    return true;
}

void
CQChromaVamp::reset()
{
    delete m_chroma;
    Chromagram::Parameters p(m_inputSampleRate);
    p.lowestOctave = m_lowestOctave;
    p.octaveCount = m_octaveCount;
    p.binsPerOctave = m_bpo;
    p.tuningFrequency = m_tuningFrequency;

    m_chroma = new Chromagram(p);

    m_haveStartTime = false;
    m_startTime = Vamp::RealTime::zeroTime;
    m_columnCount = 0;
}

size_t
CQChromaVamp::getPreferredStepSize() const
{
    return 0;
}

size_t
CQChromaVamp::getPreferredBlockSize() const
{
    return 0;
}

CQChromaVamp::OutputList
CQChromaVamp::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "chromagram";
    d.name = "Chromagram";
    d.unit = "";
    d.description = "Chromagram obtained from output of constant-Q transform, folding over each process block into a single-octave vector";
    d.hasFixedBinCount = true;
    d.binCount = m_bpo;

    if (m_chroma) {
        for (int i = 0; i < (int)d.binCount; ++i) {
            d.binNames.push_back(m_chroma->getBinName(i));
        }
    }

    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = m_inputSampleRate / (m_chroma ? m_chroma->getColumnHop() : 256);
    list.push_back(d);

    return list;
}

CQChromaVamp::FeatureSet
CQChromaVamp::process(const float *const *inputBuffers,
                      Vamp::RealTime timestamp)
{
    if (!m_chroma) {
	cerr << "ERROR: CQChromaVamp::process: "
	     << "Plugin has not been initialised"
	     << endl;
	return FeatureSet();
    }

    if (!m_haveStartTime) {
        m_startTime = timestamp;
        m_haveStartTime = true;
    }

    vector<double> data;
    for (int i = 0; i < m_blockSize; ++i) data.push_back(inputBuffers[0][i]);
    
    vector<vector<double> > chromaout = m_chroma->process(data);
    return convertToFeatures(chromaout);
}

CQChromaVamp::FeatureSet
CQChromaVamp::getRemainingFeatures()
{
    vector<vector<double> > chromaout = m_chroma->getRemainingOutput();
    return convertToFeatures(chromaout);
}

CQChromaVamp::FeatureSet
CQChromaVamp::convertToFeatures(const vector<vector<double> > &chromaout)
{
    FeatureSet returnFeatures;

    int width = chromaout.size();

    for (int i = 0; i < width; ++i) {

	vector<float> column(chromaout[i].begin(), chromaout[i].end());

	Feature feature;
	feature.hasTimestamp = true;
        feature.timestamp = m_startTime + Vamp::RealTime::frame2RealTime
            (m_columnCount * m_chroma->getColumnHop() - m_chroma->getLatency(),
             m_inputSampleRate);
	feature.values = column;
	feature.label = "";

        if (feature.timestamp >= m_startTime) {
            returnFeatures[0].push_back(feature);
        }

        ++m_columnCount;
    }

    return returnFeatures;
}

