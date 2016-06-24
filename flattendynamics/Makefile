
#CXXFLAGS	:= -DLADSPA -fpic -O3 -ffast-math -msse -msse2 -mfpmath=sse -ftree-vectorize

CXXFLAGS	:= -DLADSPA -fpic -g

SOURCES		:= flattendynamics-ladspa.cpp 
OBJECTS		:= flattendynamics-ladspa.o 
HEADERS		:= flattendynamics-ladspa.h ladspa.h

flattendynamics.so:	$(OBJECTS)
	$(CXX) $^ $(CXXFLAGS) -shared -o $@

clean:		
	rm -f $(OBJECTS)

flattendynamics.o:	$(HEADERS)
flattendynamics-ladspa.o:	$(HEADERS)

