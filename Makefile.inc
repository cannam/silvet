
SRC_DIR  := src

QMDSP_DIR ?= ../qm-dsp
VAMPSDK_DIR  ?= ../vamp-plugin-sdk

CQ_DIR	     ?= constant-q-cpp/cpp-qm-dsp

PLUGIN_EXT	?= .so

CXX	?= g++
CC	?= gcc

CFLAGS := $(CFLAGS) 
CXXFLAGS := -I. -I$(VAMPSDK_DIR) -I$(QMDSP_DIR) $(CXXFLAGS)

LDFLAGS := $(LDFLAGS) 
PLUGIN_LDFLAGS := $(LDFLAGS) $(PLUGIN_LDFLAGS)

PLUGIN	:= silvet$(PLUGIN_EXT)

VAMP_HEADERS := $(SRC_DIR)/Silvet.h
VAMP_SOURCES := $(SRC_DIR)/Silvet.cpp $(SRC_DIR)/libmain.cpp

CQ_HEADERS   := $(CQ_DIR)/CQKernel.h $(CQ_DIR)/ConstantQ.h $(CQ_DIR)/CQInterpolated.h
CQ_SOURCES   := $(CQ_DIR)/CQKernel.cpp $(CQ_DIR)/ConstantQ.cpp $(CQ_DIR)/CQInterpolated.cpp

HEADERS	     := $(VAMP_HEADERS) $(CQ_HEADERS)
SOURCES	     := $(VAMP_SOURCES) $(CQ_SOURCES)
OBJECTS	     := $(SOURCES:.cpp=.o)

LIBS	:= $(QMDSP_DIR)/libqm-dsp.a $(VAMPSDK_DIR)/libvamp-sdk.a -lpthread

all: $(PLUGIN)

$(PLUGIN):	$(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(PLUGIN_LDFLAGS)

clean:		
	rm -f $(OBJECTS)

distclean:	clean
	rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(HEADERS)

# DO NOT DELETE

src/Silvet.o: src/Silvet.h data/include/templates.h data/include/bassoon.h
src/Silvet.o: data/include/cello.h data/include/clarinet.h
src/Silvet.o: data/include/flute.h data/include/guitar.h data/include/horn.h
src/Silvet.o: data/include/oboe.h data/include/tenorsax.h
src/Silvet.o: data/include/violin.h data/include/piano-maps-SptkBGCl.h
src/Silvet.o: data/include/piano1.h data/include/piano2.h
src/Silvet.o: data/include/piano3.h
src/Silvet.o: constant-q-cpp/cpp-qm-dsp/CQInterpolated.h
src/Silvet.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
src/Silvet.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
src/libmain.o: src/Silvet.h
constant-q-cpp/cpp-qm-dsp/CQKernel.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
constant-q-cpp/cpp-qm-dsp/ConstantQ.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
constant-q-cpp/cpp-qm-dsp/ConstantQ.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/CQInterpolated.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
constant-q-cpp/cpp-qm-dsp/ConstantQ.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
