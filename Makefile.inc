
SRC_DIR  := src

VAMPSDK_DIR  ?= ../vamp-plugin-sdk

CQ_DIR	     ?= constant-q-cpp/cpp-qm-dsp
BQVEC_DIR    ?= bqvec/src
QMDSP_DIR    ?= qm-dsp

PLUGIN_EXT	?= .so

CXX	?= g++
CC	?= gcc

CFLAGS := $(CFLAGS) -I$(QMDSP_DIR)/ext/kissfft -I$(QMDSP_DIR)/ext/kissfft/tools -Dkiss_fft_scalar=double
CXXFLAGS := $(CFLAGS) -I. -I$(VAMPSDK_DIR) -I$(QMDSP_DIR) -I$(BQVEC_DIR) $(CXXFLAGS)

LDFLAGS := $(LDFLAGS) 
PLUGIN_LDFLAGS := $(LDFLAGS) $(PLUGIN_LDFLAGS)

PLUGIN	:= silvet$(PLUGIN_EXT)

PLUGIN_HEADERS := $(SRC_DIR)/Silvet.h $(SRC_DIR)/EM.h
PLUGIN_SOURCES := $(SRC_DIR)/Silvet.cpp $(SRC_DIR)/EM.cpp $(SRC_DIR)/libmain.cpp

CQ_HEADERS   := $(CQ_DIR)/CQKernel.h $(CQ_DIR)/ConstantQ.h $(CQ_DIR)/CQInterpolated.h
CQ_SOURCES   := $(CQ_DIR)/CQKernel.cpp $(CQ_DIR)/ConstantQ.cpp $(CQ_DIR)/CQInterpolated.cpp

BQVEC_HEADERS	:= $(BQVEC_DIR)/Allocators.h $(BQVEC_DIR)/Restrict.h $(BQVEC_DIR)/VectorOps.h
BQVEC_SOURCES	:= $(BQVEC_DIR)/Allocators.cpp

QMDSP_HEADERS	:= 
QMDSP_SOURCES	:= $(QMDSP_DIR)/dsp/rateconversion/Resampler.cpp $(QMDSP_DIR)/base/KaiserWindow.cpp $(QMDSP_DIR)/base/SincWindow.cpp $(QMDSP_DIR)/thread/Thread.cpp $(QMDSP_DIR)/dsp/transforms/FFT.cpp $(QMDSP_DIR)/maths/MathUtilities.cpp $(QMDSP_DIR)/ext/kissfft/kiss_fft.c $(QMDSP_DIR)/ext/kissfft/tools/kiss_fftr.c

HEADERS	     := $(PLUGIN_HEADERS) $(CQ_HEADERS) $(BQVEC_HEADERS)
SOURCES	     := $(PLUGIN_SOURCES) $(CQ_SOURCES) $(BQVEC_SOURCES) $(QMDSP_SOURCES)
OBJECTS	     := $(SOURCES:.cpp=.o)
OBJECTS	     := $(OBJECTS:.c=.o)

LIBS	:= $(VAMPSDK_DIR)/libvamp-sdk.a -lpthread

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
src/EM.o: data/include/piano-maps-SptkBGCl.h
src/libmain.o: src/Silvet.h
constant-q-cpp/cpp-qm-dsp/CQKernel.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
constant-q-cpp/cpp-qm-dsp/ConstantQ.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
constant-q-cpp/cpp-qm-dsp/ConstantQ.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/CQInterpolated.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
bqvec/src/Allocators.o: bqvec/src/Allocators.h bqvec/src/VectorOps.h
bqvec/src/Allocators.o: bqvec/src/Restrict.h
constant-q-cpp/cpp-qm-dsp/ConstantQ.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/ConstantQ.h
constant-q-cpp/cpp-qm-dsp/CQInterpolated.o: constant-q-cpp/cpp-qm-dsp/CQKernel.h
bqvec/src/Allocators.o: bqvec/src/VectorOps.h bqvec/src/Restrict.h
bqvec/src/VectorOps.o: bqvec/src/Restrict.h
