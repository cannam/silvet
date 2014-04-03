
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

CQ_HEADERS   := $(CQ_DIR)/CQKernel.h $(CQ_DIR)/ConstantQ.h
CQ_SOURCES   := $(CQ_DIR)/CQKernel.cpp $(CQ_DIR)/ConstantQ.cpp

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

