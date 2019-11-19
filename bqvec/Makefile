
# Add to VECTOR_DEFINES the relevant options for your desired
# third-party library support.
#
# Available options are
#
#  -DHAVE_IPP    Intel's Integrated Performance Primitives are available
#  -DHAVE_VDSP   Apple's Accelerate framework is available
#
# The above are optional (they affect performance, not function) and
# you may define more than one of them.
#
# The following two options trade off speed against precision for single-
# precision paths in cases where IPP and VDSP are not available:
#
#  -DUSE_POMMIER_MATHFUN Use Julien Pommier's SSE/NEON implementation
#   of sincos in 32-bit polar-to-cartesian conversion
#  -DUSE_APPROXIMATE_ATAN2 Use a quick but *very* approximate atan2
#   function in 32-bit cartesian-to-polar conversion
#
# And a handful of miscellaneous flags:
#
#  -DLACK_SINCOS  Math library lacks sincos() function
#  -DNO_COMPLEX_TYPES  Don't bother defining bq_complex_t functions
#  -DUSE_SINGLE_PRECISION_COMPLEX  Use float, not double, for bq_complex_t
#  -DNO_EXCEPTIONS  Don't throw exceptions (abort instead)
#
# Add any relevant -I flags for include paths as well.
#
# Note that you must supply the same flags when including bqvec
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqvec source files to your
# application's build system and not build a bqvec library at all.)

VECTOR_DEFINES		:=


# Add to ALLOCATOR_DEFINES options relating to aligned malloc.
# These are not usually necessary.
#
# Available options are
#
#  -DHAVE_POSIX_MEMALIGN       The posix_memalign call is available in sys/mman.h
#  -DLACK_POSIX_MEMALIGN       The posix_memalign call is not available
#
#  -DMALLOC_IS_ALIGNED         The malloc call already returns aligned memory
#  -DMALLOC_IS_NOT_ALIGNED     The malloc call does not return aligned memory
#
#  -DUSE_OWN_ALIGNED_MALLOC    If no aligned malloc is available, roll your own
#  -DAVOID_OWN_ALIGNED_MALLOC  If no aligned malloc is available, refuse to build
#
#  -DLACK_BAD_ALLOC            The C++ library lacks the std::bad_alloc exception
#
# Here "aligned" is assumed to mean "aligned enough for whatever
# vector stuff the space will be used for" which likely means at least
# 16-byte alignment.
#
# If no options are provided, we will use IPP functions if HAVE_IPP is
# defined, or else use _aligned_malloc when building with Visual C++
# on Windows, roll our own when building with some other compiler on
# Windows, use system malloc when building on OS/X, and use
# posix_memalign elsewhere.
#
# Note that you must supply the same flags when including bqvec
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqvec source files to your
# application's build system and not build a bqvec library at all.)

ALLOCATOR_DEFINES 	:= 


# Add any related includes and libraries here
#
THIRD_PARTY_INCLUDES	:=
THIRD_PARTY_LIBS	:=


# If you are including a set of bq libraries into a project, you can
# override variables for all of them (including all of the above) in
# the following file, which all bq* Makefiles will include if found

-include ../Makefile.inc-bq


# This project-local Makefile describes the source files and contains
# no routinely user-modifiable parts

include build/Makefile.inc
