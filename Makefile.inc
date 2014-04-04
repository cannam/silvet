
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

VAMP_HEADERS := $(SRC_DIR)/Silvet.h $(SRC_DIR)/EM.h
VAMP_SOURCES := $(SRC_DIR)/Silvet.cpp $(SRC_DIR)/EM.cpp $(SRC_DIR)/libmain.cpp

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

src/Silvet.o: src/Silvet.h src/EM.h
src/Silvet.o: constant-q-cpp/cpp-qm-dsp/CQInterpolated.h
src/Silvet.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
src/Silvet.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
src/EM.o: src/EM.h data/include/templates.h data/include/bassoon.h
src/EM.o: data/include/cello.h data/include/clarinet.h data/include/flute.h
src/EM.o: data/include/guitar.h data/include/horn.h data/include/oboe.h
src/EM.o: data/include/tenorsax.h data/include/violin.h
src/EM.o: data/include/piano-maps-SptkBGCl.h data/include/piano1.h
src/EM.o: data/include/piano2.h data/include/piano3.h
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
