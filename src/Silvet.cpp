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

#include "maths/MedianFilter.h"
#include "dsp/rateconversion/Resampler.h"

#include "constant-q-cpp/cpp-qm-dsp/CQInterpolated.h"

#include <vector>

#include <cstdio>

using std::vector;
using std::cerr;
using std::endl;

static int processingSampleRate = 44100;
static int processingBPO = 60;
static int processingHeight = 545;

Silvet::Silvet(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_resampler(0),
    m_cq(0)
{
}

Silvet::~Silvet()
{
    delete m_resampler;
    delete m_cq;
    for (int i = 0; i < (int)m_filterA.size(); ++i) {
        delete m_filterA[i];
        delete m_filterB[i];
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
    return list;
}

float
Silvet::getParameter(string identifier) const
{
    return 0;
}

void
Silvet::setParameter(string identifier, float value) 
{
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
    d.identifier = "transcription";
    d.name = "Transcription";
    d.description = ""; //!!!
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 2;
    d.binNames.push_back("Frequency");
    d.binNames.push_back("Velocity");
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = m_inputSampleRate / (m_cq ? m_cq->getColumnHop() : 256);
    d.hasDuration = true;
    m_notesOutputNo = list.size();
    list.push_back(d);

    d.identifier = "inputgrid";
    d.name = "Filtered time-frequency grid";
    d.description = "The pre-processed constant-Q time-frequency distribution used as input to the PLCA step";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = processingHeight;
    d.binNames.clear();
    if (m_cq) {
        char name[20];
        for (int i = 0; i < processingHeight; ++i) {
            float freq = m_cq->getBinFrequency(i + 55);
            sprintf(name, "%.1f Hz", freq);
            d.binNames.push_back(name);
        }
    }
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = 25;
    d.hasDuration = false;
    m_cqOutputNo = list.size();
    list.push_back(d);

    return list;
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

    m_cq = new CQInterpolated
	(processingSampleRate, 27.5, processingSampleRate / 3, processingBPO,
         CQInterpolated::Linear);

    for (int i = 0; i < (int)m_filterA.size(); ++i) {
        delete m_filterA[i];
        delete m_filterB[i];
    }
    m_filterA.clear();
    m_filterB.clear();
    for (int i = 0; i < processingHeight; ++i) {
        m_filterA.push_back(new MedianFilter<double>(40));
        m_filterB.push_back(new MedianFilter<double>(40));
    }
    m_columnCount = 0;
    m_reducedColumnCount = 0;
}

Silvet::FeatureSet
Silvet::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    vector<double> data;
    for (int i = 0; i < m_blockSize; ++i) data.push_back(inputBuffers[0][i]);

    if (m_resampler) {
	data = m_resampler->process(data.data(), data.size());
    }

    Grid cqout = m_cq->process(data);
    return transcribe(cqout);
}

Silvet::FeatureSet
Silvet::getRemainingFeatures()
{
    Grid cqout = m_cq->getRemainingBlocks();
    return transcribe(cqout);
}

Silvet::FeatureSet
Silvet::transcribe(const Grid &cqout)
{
    Grid filtered = preProcess(cqout);

    FeatureSet fs;

    for (int i = 0; i < (int)filtered.size(); ++i) {
        Feature f;
        for (int j = 0; j < processingHeight; ++j) {
            f.values.push_back(float(filtered[i][j]));
        }
        fs[m_cqOutputNo].push_back(f);
    }

    int width = filtered.size();

    int iterations = 12;

    for (int i = 0; i < width; ++i) {

        double sum = 0.0;
        for (int j = 0; j < processingHeight; ++j) {
            sum += filtered[i][j];
        }
        cerr << "sum = " << sum << endl;

        if (sum < 1e-5) continue;

        EM em;
        for (int j = 0; j < iterations; ++j) {
            em.iterate(filtered[i]);
        }

        //!!! now do something with the results from em!
        em.report();
    }

    return fs;
}

Silvet::Grid
Silvet::preProcess(const Grid &in)
{
    int width = in.size();

    // reduce to 100 columns per second, or one column every 441 samples

    int spacing = processingSampleRate / 100;

    Grid out;

    //!!! nb we count the CQ latency in terms of processing hops, but
    //!!! actually it isn't guaranteed to be an exact number (in fact
    //!!! it probably isn't) so this is imprecise -- fix
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

                double val = inCol[ix];
                m_filterA[j]->push(val);

                double a = m_filterA[j]->get();
                m_filterB[j]->push(std::min(a, val));

                double filtered = m_filterB[j]->get();
                outCol[j] = filtered;
            }

            // then we only use every fourth filtered column, for 25
            // columns per second in the eventual grid

            if (m_reducedColumnCount % 4 == 0) {
                out.push_back(outCol);
            }

            ++m_reducedColumnCount;
        }

        ++m_columnCount;
    }

    return out;
}
    
