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

#include "LiveInstruments.h"

#include "data/include/templates.h"

#include <iostream>

using namespace std;

InstrumentPack
LiveAdapter::adapt(const InstrumentPack &original)
{
    vector<InstrumentPack::Templates> templates;

//            cerr << "LiveAdapter: reduced template height is " << SILVET_TEMPLATE_HEIGHT/5 << endl;

    bool merge = false;
    // The live template for piano has only one piano in it, so as
    // to process faster. We make it by averaging the originals
    if (original.name == "Piano") {
        merge = true;
    }

    InstrumentPack::Templates t;
    bool first = true;

    // The live instrument template is one octave shorter than the
    // original, as well as having only 12 bpo instead of 60
    int height = SILVET_TEMPLATE_HEIGHT/5 - 12;
    
    for (const auto &origt: original.templates) {

	t.lowestNote = origt.lowestNote;
	t.highestNote = origt.highestNote;

        t.data.resize(SILVET_TEMPLATE_NOTE_COUNT);

	for (int j = 0; j < SILVET_TEMPLATE_NOTE_COUNT; ++j) {

	    t.data[j].resize(height);

            if (j >= t.lowestNote && j <= t.highestNote) {
                for (int k = 0; k < height; ++k) {
                    // This is the index of the middle template, no. 2
                    // of 5, in each 5-bin semitone series. (We add 4
                    // because there are 2 blank shift-space slots at
                    // the start of the template data.)
                    t.data[j][k] += origt.data[j][60 + k * 5 + 4];
                }
            }
        }

        if (!merge) {
            templates.push_back(t);
            t = InstrumentPack::Templates();
        }

        first = false;
    }

    if (merge) {
        templates.push_back(t);
    }

    // re-normalise

    for (auto &t: templates) {
        for (auto &d: t.data) {
            float sum = 0.f;
            for (auto v: d) sum += v;
            for (auto &v: d) v /= sum;
        }
    }

    InstrumentPack live(templates[0].lowestNote,
			templates[0].highestNote,
			original.name,
			templates);

    live.templateHeight = height;
    live.templateMaxShift = 0;
    live.templateSize = height;
    
    live.maxPolyphony = original.maxPolyphony;
    live.pitchSparsity = original.pitchSparsity;
    live.sourceSparsity = original.sourceSparsity;
    live.levelThreshold = original.levelThreshold / 20;

    return live;
}

vector<InstrumentPack>
LiveAdapter::adaptAll(const vector<InstrumentPack> &v)
{
    vector<InstrumentPack> out;
    for (int i = 0; i < (int)v.size(); ++i) {
        InstrumentPack p(LiveAdapter::adapt(v[i]));
        out.push_back(p);
    }
    return out;
}
