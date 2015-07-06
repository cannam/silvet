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

#include "Instruments.h"

#include "data/include/templates.h"

#include <iostream>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

InstrumentPack::InstrumentPack(int lowest, int highest,
                               std::string n, std::vector<Templates> tt) :
    templateNoteCount(SILVET_TEMPLATE_NOTE_COUNT),
    templateHeight(SILVET_TEMPLATE_HEIGHT),
    templateMaxShift(SILVET_TEMPLATE_MAX_SHIFT),
    templateSize(SILVET_TEMPLATE_SIZE),
    lowestNote(lowest),
    highestNote(highest),
    maxPolyphony(5),
    pitchSparsity(1.1),
    sourceSparsity(1.2),
    levelThreshold(5),
    name(n),
    templates(tt)
{
}

const char *simpleInstruments[] = {
    // Each instrument has two consecutive slots, one for the pack
    // name and one for the template to look up
    "Guitar", "guitar",
    "Violin", "violin",
    "Viola", "viola",
    "Cello", "cello",
    "Horn", "horn",
    "Flute", "flute",
    "Oboe", "oboe",
    "Clarinet", "clarinet",
    "Tenor Sax", "tenorsax",
    "Bassoon", "bassoon",
};

static bool
isBowedString(int i)
{
    string tname(simpleInstruments[i+1]);
    return tname == "violin"
	|| tname == "viola"
	|| tname == "cello"
	;
}

static bool
isWind(int i)
{
    string tname(simpleInstruments[i+1]);
    return tname == "horn"
	|| tname == "flute"
	|| tname == "oboe"
	|| tname == "clarinet"
	|| tname == "tenorsax"
	|| tname == "bassoon"
	;
}

static InstrumentPack::Templates
templatesFor(string name)
{
    for (int i = 0; i < SILVET_TEMPLATE_COUNT; ++i) {

	if (string(silvet_templates[i].name) == name) {

	    silvet_template_t *t = &silvet_templates[i];

	    InstrumentPack::Templates rt;
	    rt.lowestNote = t->lowest;
	    rt.highestNote = t->highest;
	    rt.data = vector<vector<float> >
		(SILVET_TEMPLATE_NOTE_COUNT,
		 vector<float>(SILVET_TEMPLATE_SIZE, 0.f));

	    for (int j = 0; j < SILVET_TEMPLATE_NOTE_COUNT; ++j) {
		for (int k = 0; k < SILVET_TEMPLATE_SIZE; ++k) {
		    rt.data[j][k] = t->data[j][k];
		}
	    }

	    return rt;
	}
    }

    return InstrumentPack::Templates();
}

static bool 
isOK(InstrumentPack &d)
{
    if (d.name == "") {
	cerr << "ERROR: Silvet::InstrumentPack: Empty name in instrument definition" << endl;
	return false;
    }
    if (d.templates.empty()) {
	cerr << "ERROR: Silvet::InstrumentPack: Instrument definition \"" << d.name << "\" contains no templates!" << endl;
	return false;
    }
    for (int i = 0; i < int(d.templates.size()); ++i) {
	if (d.templates[i].data.empty()) {
	    cerr << "ERROR: Silvet::InstrumentPack: Instrument definition \"" << d.name << "\" contains one or more empty templates!" << endl;
	    return false;
	}
    }
    return true;
}

vector<InstrumentPack> 
InstrumentPack::listInstrumentPacks()
{
    vector<InstrumentPack> ii;

    vector<Templates> allTemplates;
    allTemplates.push_back(templatesFor("piano1"));

    vector<Templates> pianoTemplates;
    pianoTemplates.push_back(templatesFor("piano1"));
    pianoTemplates.push_back(templatesFor("piano2"));
    pianoTemplates.push_back(templatesFor("piano3"));
    pianoTemplates.push_back(templatesFor("pianorwc"));
    InstrumentPack piano(silvet_templates_lowest_note,
			 silvet_templates_highest_note,
			 "Piano",
			 pianoTemplates);
    piano.maxPolyphony = 6;
    piano.levelThreshold = 5;
    piano.pitchSparsity = 1.0;
    piano.sourceSparsity = 1.0;
    if (isOK(piano)) {
	ii.push_back(piano);
    }

    vector<Templates> stringTemplates;
    vector<Templates> windTemplates;

    for (int i = 0;
	 i < int(sizeof(simpleInstruments)/sizeof(simpleInstruments[0]));
	 i += 2) {

	vector<Templates> tt;
	Templates t = templatesFor(simpleInstruments[i+1]);
	tt.push_back(t);
	allTemplates.push_back(t);
	InstrumentPack instr(t.lowestNote,
			     t.highestNote,
			     simpleInstruments[i],
			     tt);

        instr.maxPolyphony = 5;
        instr.levelThreshold = 6;

	if (isBowedString(i)) {
            instr.maxPolyphony = 2;
            instr.levelThreshold = 5;
	    stringTemplates.push_back(t);
	}
	if (isWind(i)) {
            instr.maxPolyphony = 2;
	    windTemplates.push_back(t);
	}
	if (isOK(instr)) {
	    ii.push_back(instr);
	}
    }

    InstrumentPack all(silvet_templates_lowest_note,
		       silvet_templates_highest_note,
		       "Multiple or unknown instruments",
		       allTemplates);
    all.maxPolyphony = 5;
    all.levelThreshold = 6;
    if (isOK(all)) {
	ii.insert(ii.begin(), all);
    }

    InstrumentPack strings(silvet_templates_lowest_note,  // cello
			   silvet_templates_highest_note, // violin
			   "String quartet",
			   stringTemplates);
    strings.maxPolyphony = 6;
    strings.levelThreshold = 3;
    if (isOK(strings)) {
	ii.push_back(strings);
    }

    InstrumentPack winds(silvet_templates_lowest_note,   // basson
			 silvet_templates_highest_note,  // flute
			 "Wind ensemble",
			 windTemplates);
    winds.maxPolyphony = 5;
    winds.levelThreshold = 5;
    if (isOK(winds)) {
	ii.push_back(winds);
    }

    return ii;
}

