
SRC_DIR  := src

VAMPSDK_DIR  ?= ../vamp-plugin-sdk

CQ_DIR	     ?= constant-q-cpp
BQVEC_DIR    ?= bqvec/src

PLUGIN_EXT	?= .so

CXX	?= g++
CC	?= gcc

CFLAGS := $(CFLAGS)
CXXFLAGS := $(CFLAGS) -I. -I$(VAMPSDK_DIR) -I$(CQ_DIR) -I$(BQVEC_DIR) $(CXXFLAGS)

LDFLAGS := $(LDFLAGS) 
PLUGIN_LDFLAGS := $(LDFLAGS) $(PLUGIN_LDFLAGS)

PLUGIN	:= silvet$(PLUGIN_EXT)

PLUGIN_HEADERS := \
	$(SRC_DIR)/Silvet.h \
	$(SRC_DIR)/EM.h \
	$(SRC_DIR)/Instruments.h \
	$(SRC_DIR)/AgentFeeder.h \
	$(SRC_DIR)/AgentFeederMono.h \
	$(SRC_DIR)/AgentFeederPoly.h \
	$(SRC_DIR)/AgentHypothesis.h \
	$(SRC_DIR)/NoteHypothesis.h

PLUGIN_SOURCES := \
	$(SRC_DIR)/Silvet.cpp \
	$(SRC_DIR)/EM.cpp \
	$(SRC_DIR)/Instruments.cpp \
	$(SRC_DIR)/NoteHypothesis.cpp \
	$(SRC_DIR)/libmain.cpp

BQVEC_HEADERS	:= \
	$(BQVEC_DIR)/Allocators.h \
	$(BQVEC_DIR)/Restrict.h \
	$(BQVEC_DIR)/VectorOps.h

BQVEC_SOURCES	:= \
	$(BQVEC_DIR)/Allocators.cpp

HEADERS	     := $(PLUGIN_HEADERS) $(BQVEC_HEADERS)
SOURCES	     := $(PLUGIN_SOURCES) $(BQVEC_SOURCES)
OBJECTS	     := $(SOURCES:.cpp=.o)
OBJECTS	     := $(OBJECTS:.c=.o)

LIBS	:= $(CQ_DIR)/libcq.a $(VAMPSDK_DIR)/libvamp-sdk.a

all: constant-q-cpp $(PLUGIN)

$(PLUGIN):	$(OBJECTS) $(LIBS)
	$(CXX) -o $@ $^ $(LIBS) $(PLUGIN_LDFLAGS)

clean:		
	rm -f $(OBJECTS)

distclean:	clean
	rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(HEADERS)

# DO NOT DELETE

src/Silvet.o: src/Silvet.h src/MedianFilter.h src/Instruments.h
src/Silvet.o: src/NoteHypothesis.h src/AgentHypothesis.h src/EM.h
src/Silvet.o: src/AgentFeederPoly.h src/AgentFeeder.h
src/Silvet.o: constant-q-cpp/src/dsp/Resampler.h
src/EM.o: src/EM.h src/Instruments.h
src/Instruments.o: src/Instruments.h data/include/templates.h
src/Instruments.o: data/include/bassoon.h data/include/cello.h
src/Instruments.o: data/include/clarinet.h data/include/flute.h
src/Instruments.o: data/include/guitar.h data/include/horn.h
src/Instruments.o: data/include/oboe.h data/include/tenorsax.h
src/Instruments.o: data/include/violin.h data/include/piano1.h
src/Instruments.o: data/include/piano2.h data/include/piano3.h
src/NoteHypothesis.o: src/NoteHypothesis.h src/AgentHypothesis.h
src/NoteHypothesis.o: src/AgentFeederPoly.h src/AgentFeeder.h
src/libmain.o: src/Silvet.h src/MedianFilter.h src/Instruments.h
src/libmain.o: src/NoteHypothesis.h src/AgentHypothesis.h
bqvec/src/Allocators.o: bqvec/src/Allocators.h bqvec/src/VectorOps.h
bqvec/src/Allocators.o: bqvec/src/Restrict.h
src/Silvet.o: src/MedianFilter.h src/Instruments.h src/NoteHypothesis.h
src/Silvet.o: src/AgentHypothesis.h
src/AgentFeeder.o: src/AgentHypothesis.h
src/AgentFeederMono.o: src/AgentFeeder.h src/AgentHypothesis.h
src/AgentFeederPoly.o: src/AgentFeeder.h src/AgentHypothesis.h
src/NoteHypothesis.o: src/AgentHypothesis.h
bqvec/src/Allocators.o: bqvec/src/VectorOps.h bqvec/src/Restrict.h
bqvec/src/VectorOps.o: bqvec/src/Restrict.h
