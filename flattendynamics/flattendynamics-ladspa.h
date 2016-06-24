/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#ifndef FLATTENDYNAMICS_LADSPA_H
#define FLATTENDYNAMICS_LADSPA_H

#include <ladspa.h>

class FlattenDynamics
{
public:
    static const LADSPA_Descriptor *getDescriptor(unsigned long index);

    enum Port {
        AudioInputPort     = 0,
	AudioOutputPort    = 1,
        GainOutputPort     = 2,
	PortCount          = 3,
    };

    // Rest of the public interface is for use when constructing the
    // class directly, rather than through LADSPA

    FlattenDynamics(int sampleRate);
    ~FlattenDynamics();

    void connectInputPort(Port p, const float *addr);
    void connectOutputPort(Port p, float *addr);
    void reset();
    void process(int nsamples);

private:
    static const char *const portNames[PortCount];
    static const LADSPA_PortDescriptor ports[PortCount];
    static const LADSPA_PortRangeHint hints[PortCount];
    static const LADSPA_Properties properties;
    static const LADSPA_Descriptor ladspaDescriptor;

    static LADSPA_Handle ladspaInstantiate(const LADSPA_Descriptor *, unsigned long);
    static void ladspaConnectPort(LADSPA_Handle, unsigned long, LADSPA_Data *);
    static void ladspaActivate(LADSPA_Handle);
    static void ladspaRun(LADSPA_Handle, unsigned long);
    static void ladspaDeactivate(LADSPA_Handle);
    static void ladspaCleanup(LADSPA_Handle);

    float processSingle(float sample);
    void updateRMS(float);
    void updateParameters();

    int m_sampleRate;
    const float *m_input;
    float *m_output;
    float *m_pgain;

    float *m_history;
    int m_histlen;
    int m_histwrite;
    int m_histread;

    double m_sumOfSquares;
    float m_rms;
    float m_maxRms;
    float m_gain;
};

#endif
