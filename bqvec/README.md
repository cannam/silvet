
bqvec
=====

A small library for vector management and arithmetic in C++ using raw
C pointer arrays, designed for simple audio buffer-shuffling. Also
includes aligned malloc wrappers and a lock-free ring buffer.

The code can call out to vector arithmetic helpers (IPP, vDSP) in some
places, and has loops written with an eye to auto-vectorising
compilers, but mostly this is a convenience library rather than for
performance -- it initially exists to give a fairly consistent API to
useful functions over audio buffer arrays.

This code originated as part of the Rubber Band Library written by the
same authors (see https://hg.sr.ht/~breakfastquay/rubberband/).
It has been pulled out into a separate library and relicensed under a
more permissive licence.

Generally expected to be vendored in to local project builds rather
than being installed as a system library.

C++ standard required: C++98 (does not use C++11 or newer features)

 * To compile on Linux: make
   
 * To compile on macOS: make -f build/Makefile.osx

 * To build and run tests: as above, but add the "test" target -
   requires Boost test headers installed

[![Build Status](https://travis-ci.org/breakfastquay/bqvec.svg?branch=master)](https://travis-ci.org/breakfastquay/bqvec)

Copyright 2007-2018 Particular Programs Ltd. See the file COPYING for
(BSD/MIT-style) licence terms.

