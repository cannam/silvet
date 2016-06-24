/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "flattendynamics-ladspa.h"

#include <iostream>
#include <cmath>

using std::cerr;
using std::endl;

const float historySeconds = 4.f;
const float catchUpSeconds = 0.2f;
const float targetMaxRMS = 0.05f;
const float rmsMaxDecay = 0.999f; // per sample
const float maxGain = 10.f;

const char *const
FlattenDynamics::portNames[PortCount] =
{
    "Input",
    "Output",
    "Gain",
};

const LADSPA_PortDescriptor 
FlattenDynamics::ports[PortCount] =
{
    LADSPA_PORT_INPUT  | LADSPA_PORT_AUDIO,
    LADSPA_PORT_OUTPUT | LADSPA_PORT_AUDIO,
    LADSPA_PORT_OUTPUT | LADSPA_PORT_CONTROL,
};

const LADSPA_PortRangeHint 
FlattenDynamics::hints[PortCount] =
{
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
};

const LADSPA_Properties
FlattenDynamics::properties = LADSPA_PROPERTY_HARD_RT_CAPABLE;

const LADSPA_Descriptor 
FlattenDynamics::ladspaDescriptor =
{
    0xf0b375, // "Unique" ID
    "flattendynamics", // Label
    properties,
    "Flatten Dynamics", // Name
    "Queen Mary, University of London", // Maker
    "BSD", // Copyright
    PortCount,
    ports,
    portNames,
    hints,
    0, // Implementation data
    ladspaInstantiate,
    ladspaConnectPort,
    ladspaActivate,
    ladspaRun,
    0, // Run adding
    0, // Set run adding gain
    ladspaDeactivate,
    ladspaCleanup
};

const LADSPA_Descriptor *
FlattenDynamics::getDescriptor(unsigned long index)
{
    if (index == 0) return &ladspaDescriptor;
    return 0;
}

FlattenDynamics::FlattenDynamics(int sampleRate) :
    m_sampleRate(sampleRate),
    m_input(0),
    m_output(0),
    m_pgain(0),
    m_history(0),
    m_histlen(0),
    m_histwrite(0),
    m_histread(0),
    m_sumOfSquares(0.f),
    m_rms(0.f),
    m_maxRms(0.f),
    m_gain(1.f)
{
    reset();
}

FlattenDynamics::~FlattenDynamics()
{
    delete[] m_history;
}
    
LADSPA_Handle
FlattenDynamics::ladspaInstantiate(const LADSPA_Descriptor *, unsigned long rate)
{
    FlattenDynamics *flatten = new FlattenDynamics(rate);
    return flatten;
}

void
FlattenDynamics::ladspaConnectPort(LADSPA_Handle handle,
                                   unsigned long port, 
                                   LADSPA_Data *location)
{
    FlattenDynamics *flatten = (FlattenDynamics *)handle;
    if (ports[port] & LADSPA_PORT_INPUT) {
        flatten->connectInputPort(Port(port), location);
    } else {
        flatten->connectOutputPort(Port(port), location);
    }
}

void
FlattenDynamics::ladspaActivate(LADSPA_Handle handle)
{
    FlattenDynamics *flatten = (FlattenDynamics *)handle;
    flatten->reset();
}

void
FlattenDynamics::ladspaRun(LADSPA_Handle handle, unsigned long samples)
{
    FlattenDynamics *flatten = (FlattenDynamics *)handle;
    flatten->process(samples);
}

void
FlattenDynamics::ladspaDeactivate(LADSPA_Handle handle)
{
    ladspaActivate(handle); // both functions just reset the plugin
}

void
FlattenDynamics::ladspaCleanup(LADSPA_Handle handle)
{
    delete (FlattenDynamics *)handle;
}

void
FlattenDynamics::connectInputPort(Port p, const float *location)
{
    const float **ports[PortCount] = {
	&m_input, 0, 0,
    };

    *ports[int(p)] = location;
}

void
FlattenDynamics::connectOutputPort(Port p, float *location)
{
    float **ports[PortCount] = {
	0, &m_output, &m_pgain,
    };

    *ports[int(p)] = location;
}

void
FlattenDynamics::reset()
{
    delete[] m_history;
    m_histlen = int(round(m_sampleRate * historySeconds));
    if (m_histlen < 1) m_histlen = 1;
    m_history = new float[m_histlen];
    for (int i = 0; i < m_histlen; ++i) {
        m_history[i] = 0.f;
    }
    m_histwrite = 0;
    m_histread = 0;

    m_sumOfSquares = 0.0;
    m_rms = 0.f;
    m_maxRms = 0.f;
    m_gain = 1.f;
}

void
FlattenDynamics::updateParameters()
{
    if (m_pgain) *m_pgain = m_gain;
}

void
FlattenDynamics::process(int sampleCount)
{
    if (!m_input || !m_output) return;

    updateParameters();

    for (int i = 0; i < sampleCount; ++i) {
        m_output[i] = processSingle(m_input[i]);
    }
}

float
FlattenDynamics::processSingle(float f)
{
    updateRMS(f);

    if (m_rms == 0.f) {
        return f;
    }

    if (m_rms >= m_maxRms) {
        m_maxRms = m_rms;
    } else {
        m_maxRms = m_rms + (m_maxRms - m_rms) * rmsMaxDecay;
    }

    float targetGain = targetMaxRMS / m_maxRms;

    if (targetGain > maxGain) {
        targetGain = maxGain;
    }

    float catchUpSamples = catchUpSeconds * m_sampleRate;
    // asymptotic, could improve?
    m_gain = m_gain + (targetGain - m_gain) / catchUpSamples;

    if (fabsf(f) * m_gain > 1.f) {
        m_gain = 1.f / fabsf(f);
    }

//    cerr << "target gain = " << targetGain << ", gain = " << m_gain << endl;
    return f * m_gain;
}

void
FlattenDynamics::updateRMS(float f)
{
    // We update the RMS values by maintaining a sum-of-last-n-squares
    // total (which the RMS is the square root of 1/n of) and
    // recording the last n samples of history in a circular
    // buffer. When a sample drops off the start of the history, we
    // remove the square of it from the sum-of-squares total; then we
    // add the square of the new sample.

    int nextWrite = (m_histwrite + 1) % m_histlen;

    float lose = 0.f;

    if (nextWrite == m_histread) {
        // full
        lose = m_history[m_histread];
        m_histread = (m_histread + 1) % m_histlen;
    }

    m_history[m_histwrite] = f;
    m_histwrite = nextWrite;

    int fill = (m_histwrite - m_histread + m_histlen) % m_histlen;

    m_sumOfSquares -= lose * lose;
    m_sumOfSquares += f * f;

    m_rms = sqrt(m_sumOfSquares / fill);
//    cerr << "rms = " << m_rms << " (from " << fill << " samples of " << m_histlen << ", latest " << f << ")" << endl;
}

const LADSPA_Descriptor *
ladspa_descriptor(unsigned long ix)
{
    return FlattenDynamics::getDescriptor(ix);
}

