/* Do not edit: this file was automatically generated by generateTemplatesC */

#ifndef SILVET_DATA_TEMPLATES_H
#define SILVET_DATA_TEMPLATES_H

/* note: intended to parse as both C and C++ */

#define SILVET_TEMPLATE_COUNT      10   /* Number of instruments */
#define SILVET_TEMPLATE_NOTE_COUNT 88   /* Number of notes per instrument */ 
#define SILVET_TEMPLATE_HEIGHT     545  /* Frequency bins per template */
#define SILVET_TEMPLATE_MAX_SHIFT  2    /* Zeros at either end of template */ 
#define SILVET_TEMPLATE_SIZE       549  /* Height + 2 * max shift space */ 

typedef struct {
    const char *name;
    int lowest;
    int highest;
    float data[SILVET_TEMPLATE_NOTE_COUNT][SILVET_TEMPLATE_SIZE];
} silvet_template_t;

static int silvet_templates_lowest_note = 15;
static int silvet_templates_highest_note = 72;

static silvet_template_t silvet_templates[SILVET_TEMPLATE_COUNT] = {
#include "bassoon.h"
#include "cello.h"
#include "clarinet.h"
#include "flute.h"
#include "guitar.h"
#include "horn.h"
#include "oboe.h"
#include "tenorsax.h"
#include "violin.h"
#include "piano-maps-SptkBGCl.h"
};

#endif
