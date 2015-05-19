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

#include <vector>
#include <string>

#ifndef SILVET_INSTRUMENTS_H
#define SILVET_INSTRUMENTS_H

/**
 * Define an instrument pack, i.e. a group of templates that are made
 * available as a single preset at the user interface level. A pack
 * might contain only a single instrument template (e.g. bassoon), or
 * it may be a compound of several templates (e.g. different piano
 * recordings forming a single piano pack), or it may be a group of
 * distinct instrument templates (e.g. a pack containing all supported
 * instruments, or potentially groupings such as string quartet or
 * rock band).
 */
class InstrumentPack
{
public:
    int templateNoteCount;
    int templateHeight;
    int templateMaxShift;
    int templateSize; // height plus space for shift at either end

    int lowestNote;
    int highestNote;

    int maxPolyphony; // realistic practical limit, not a theoretical one
    float pitchSparsity;
    float sourceSparsity;
    float levelThreshold;

    std::string name;

    struct Templates {
	int lowestNote;
	int highestNote;
	// templateNoteCount * templateSize
	std::vector<std::vector<float> > data;
    };
    
    std::vector<Templates> templates;

    static std::vector<InstrumentPack> listInstrumentPacks();

private:
    InstrumentPack(int lowest, int highest, std::string n,
		   std::vector<Templates> tt);
    
    friend class LiveAdapter;
};

#endif
