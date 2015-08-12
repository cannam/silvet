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
    
    for (const auto &origt: original.templates) {

	t.lowestNote = origt.lowestNote;
	t.highestNote = origt.highestNote;
        t.data.resize(origt.data.size());

	for (int j = 0; j < int(origt.data.size()); ++j) {

	    t.data[j].resize(SILVET_TEMPLATE_HEIGHT/5);

	    for (int k = 0; k < SILVET_TEMPLATE_HEIGHT/5; ++k) {

                if (!merge || first) {
                    t.data[j][k] = 0.f;
                }

                for (int m = 0; m < 5; ++m) {
                    t.data[j][k] += origt.data[j][k * 5 + m + 2];
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
            if (sum > 0.f) {
                for (auto &v: d) v /= sum;
            }
        }
    }
    
    InstrumentPack live(original.lowestNote,
			original.highestNote,
			original.name,
			templates);

    live.templateHeight = SILVET_TEMPLATE_HEIGHT/5;
    live.templateMaxShift = 0;
    live.templateSize = live.templateHeight;
    
    live.maxPolyphony = original.maxPolyphony;
    live.pitchSparsity = original.pitchSparsity;
    live.sourceSparsity = original.sourceSparsity;
    live.levelThreshold = original.levelThreshold / 15;

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
