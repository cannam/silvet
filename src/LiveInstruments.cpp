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
            
    for (vector<InstrumentPack::Templates>::const_iterator i =
	     original.templates.begin();
	 i != original.templates.end(); ++i) {

	InstrumentPack::Templates t;
	t.lowestNote = i->lowestNote;
	t.highestNote = i->highestNote;
	t.data.resize(i->data.size());

	for (int j = 0; j < int(i->data.size()); ++j) {

	    t.data[j].resize(SILVET_TEMPLATE_HEIGHT/5);

            float sum = 0.f;

	    for (int k = 0; k < SILVET_TEMPLATE_HEIGHT/5; ++k) {

                t.data[j][k] = 0.f;

                for (int m = 0; m < 5; ++m) {
                    t.data[j][k] += i->data[j][k * 5 + m + 2];
                }
                
                sum += t.data[j][k];
	    }
            
	    // re-normalise
            if (sum > 0.f) {
                for (int k = 0; k < (int)t.data[j].size(); ++k) {
                    t.data[j][k] *= 1.f / sum;
                }
            }
	}

	templates.push_back(t);
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
