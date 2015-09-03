
SRC_DIR  := src

VAMPSDK_DIR  ?= ../vamp-plugin-sdk

CQ_DIR	     ?= constant-q-cpp
BQVEC_DIR    ?= bqvec
FD_DIR	     ?= flattendynamics

PLUGIN_EXT	?= .so

CXX	?= g++
CC	?= gcc

CFLAGS := $(CFLAGS)
CXXFLAGS := $(CFLAGS) -I. -I$(VAMPSDK_DIR) -I$(CQ_DIR) -I$(BQVEC_DIR) -I$(BQVEC_DIR) -I$(FD_DIR) $(CXXFLAGS) 

LDFLAGS := $(LDFLAGS) 
PLUGIN_LDFLAGS := $(LDFLAGS) $(PLUGIN_LDFLAGS)

PLUGIN	:= silvet$(PLUGIN_EXT)

PLUGIN_HEADERS := $(SRC_DIR)/Silvet.h $(SRC_DIR)/EM.h $(SRC_DIR)/Instruments.h $(SRC_DIR)/LiveInstruments.h
PLUGIN_SOURCES := $(SRC_DIR)/Silvet.cpp $(SRC_DIR)/EM.cpp $(SRC_DIR)/Instruments.cpp $(SRC_DIR)/LiveInstruments.cpp $(SRC_DIR)/libmain.cpp

BQVEC_HEADERS	:= $(BQVEC_DIR)/Allocators.h $(BQVEC_DIR)/Restrict.h $(BQVEC_DIR)/VectorOps.h
BQVEC_SOURCES	:= $(BQVEC_DIR)/src/Allocators.cpp

FD_HEADERS	:= $(FD_DIR)/flattendynamics-ladspa.h
FD_SOURCES	:= $(FD_DIR)/flattendynamics-ladspa.cpp

HEADERS	     := $(PLUGIN_HEADERS) $(BQVEC_HEADERS) $(FD_HEADERS)
SOURCES	     := $(PLUGIN_SOURCES) $(BQVEC_SOURCES) $(FD_SOURCES)
OBJECTS	     := $(SOURCES:.cpp=.o)
OBJECTS	     := $(OBJECTS:.c=.o)

LIBS	:= $(CQ_DIR)/libcq.a $(VAMPSDK_DIR)/libvamp-sdk.a

all: constant-q-cpp $(PLUGIN)

.PHONY: constant-q-cpp
constant-q-cpp: 
	$(MAKE) -C $@ -f Makefile$(MAKEFILE_EXT) libcq.a

$(PLUGIN):	$(OBJECTS) $(LIBS)
	$(CXX) -o $@ $^ $(LIBS) $(PLUGIN_LDFLAGS)

clean:		
	rm -f $(OBJECTS)
	$(MAKE) -C constant-q-cpp -f Makefile$(MAKEFILE_EXT) clean

distclean:	clean
	rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(HEADERS)

# DO NOT DELETE

src/Silvet.o: src/Silvet.h src/MedianFilter.h src/Instruments.h src/EM.h
src/Silvet.o: constant-q-cpp/src/dsp/Resampler.h src/LiveInstruments.h
src/EM.o: src/EM.h src/Instruments.h
src/Instruments.o: src/Instruments.h data/include/templates.h
src/Instruments.o: data/include/bassoon.h data/include/cello.h
src/Instruments.o: data/include/clarinet.h data/include/flute.h
src/Instruments.o: data/include/guitar.h data/include/horn.h
src/Instruments.o: data/include/oboe.h data/include/tenorsax.h
src/Instruments.o: data/include/violin.h data/include/viola.h
src/Instruments.o: data/include/piano1.h data/include/piano2.h
src/Instruments.o: data/include/piano3.h data/include/pianorwc.h
src/LiveInstruments.o: src/LiveInstruments.h src/Instruments.h
src/LiveInstruments.o: data/include/templates.h data/include/bassoon.h
src/LiveInstruments.o: data/include/cello.h data/include/clarinet.h
src/LiveInstruments.o: data/include/flute.h data/include/guitar.h
src/LiveInstruments.o: data/include/horn.h data/include/oboe.h
src/LiveInstruments.o: data/include/tenorsax.h data/include/violin.h
src/LiveInstruments.o: data/include/viola.h data/include/piano1.h
src/LiveInstruments.o: data/include/piano2.h data/include/piano3.h
src/LiveInstruments.o: data/include/pianorwc.h
src/libmain.o: src/Silvet.h src/MedianFilter.h src/Instruments.h
flattendynamics/flattendynamics-ladspa.o: flattendynamics/flattendynamics-ladspa.h
src/Silvet.o: src/MedianFilter.h src/Instruments.h
src/LiveInstruments.o: src/Instruments.h
