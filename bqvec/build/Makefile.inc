
SRC_DIR	:= src
TEST_DIR := test
HEADER_DIR := bqvec

SOURCES	:= $(wildcard $(SRC_DIR)/*.cpp)
HEADERS	:= $(wildcard $(HEADER_DIR)/*.h) $(wildcard $(SRC_DIR)/*.h)

OBJECTS	:= $(SOURCES:.cpp=.o)
OBJECTS	:= $(OBJECTS:.c=.o)

TIMINGS_SOURCES	:= $(TEST_DIR)/Timings.cpp
TIMINGS_OBJECTS	:= $(TIMINGS_SOURCES:.cpp=.o)

TEST_SOURCES	:= $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS	:= $(TEST_SOURCES:.cpp=.o)

OPTFLAGS := -O3 -ffast-math

CXXFLAGS := -std=c++98 -fpic -Wall -Wextra -Werror $(VECTOR_DEFINES) $(ALLOCATOR_DEFINES) -I. $(THIRD_PARTY_INCLUDES) -I$(HEADER_DIR) $(OPTFLAGS)

LIBRARY	:= libbqvec.a

all:	$(LIBRARY) timings

test:	$(LIBRARY) timings test-allocators test-vectorops test-vectorops-complex
	./test-allocators && ./test-vectorops && ./test-vectorops-complex

valgrind:	$(LIBRARY) timings test-allocators test-vectorops test-vectorops-complex
	valgrind ./test-allocators && valgrind ./test-vectorops && valgrind ./test-vectorops-complex

$(LIBRARY):	$(OBJECTS)
	$(AR) rc $@ $^

timings: $(TIMINGS_OBJECTS) $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(THIRD_PARTY_LIBS)

test-allocators:	test/TestAllocators.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework -L. -lbqvec $(THIRD_PARTY_LIBS)

test-vectorops:	test/TestVectorOps.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework -L. -lbqvec $(THIRD_PARTY_LIBS)

test-vectorops-complex:	test/TestVectorOpsComplex.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework -L. -lbqvec $(THIRD_PARTY_LIBS)

clean:		
	rm -f $(OBJECTS) $(TEST_OBJECTS) $(TIMINGS_OBJECTS)

distclean:	clean
	rm -f $(LIBRARY) test-allocators test-vectorops test-vectorops-complex

depend:
	makedepend -Y -fbuild/Makefile.inc $(TIMINGS_SOURCES) $(TEST_SOURCES) $(HEADERS) $(SOURCES)


# DO NOT DELETE

test/Timings.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h bqvec/Restrict.h
test/Timings.o: bqvec/ComplexTypes.h
test/Timings.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h bqvec/Restrict.h
test/Timings.o: bqvec/ComplexTypes.h
test/TestVectorOpsComplex.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h
test/TestVectorOpsComplex.o: bqvec/Restrict.h bqvec/ComplexTypes.h
test/TestVectorOpsComplex.o: bqvec/VectorOps.h
test/TestVectorOps.o: bqvec/VectorOps.h
test/TestAllocators.o: bqvec/VectorOps.h bqvec/Allocators.h
bqvec/RingBuffer.o: bqvec/Barrier.h bqvec/Allocators.h bqvec/Restrict.h
bqvec/RingBuffer.o: bqvec/VectorOps.h
bqvec/VectorOpsComplex.o: bqvec/VectorOps.h bqvec/Restrict.h
bqvec/VectorOpsComplex.o: bqvec/ComplexTypes.h
bqvec/VectorOps.o: bqvec/Restrict.h
src/Allocators.o: bqvec/Allocators.h
src/Barrier.o: bqvec/Barrier.h
