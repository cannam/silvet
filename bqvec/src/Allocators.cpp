/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    bqvec

    A small library for vector arithmetic and allocation in C++ using
    raw C pointer arrays.

    Copyright 2007-2016 Particular Programs Ltd.

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
    ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of Chris Cannam and
    Particular Programs Ltd shall not be used in advertising or
    otherwise to promote the sale, use or other dealings in this
    Software without prior written authorization.
*/

#include "Allocators.h"

#ifdef HAVE_IPP
#include <ipps.h>
#endif

#include <iostream>
#include <climits>

using std::cerr;
using std::endl;

namespace breakfastquay {

#ifdef HAVE_IPP

template <>
float *allocate(size_t count)
{
    if (count > INT_MAX) {
#ifndef NO_EXCEPTIONS
        throw std::length_error("Size overflow in allocate");
#else
        abort();
#endif
    }
    
    float *ptr = ippsMalloc_32f(int(count));
    if (!ptr) {
#ifndef NO_EXCEPTIONS
        throw (std::bad_alloc());
#else
        abort();
#endif
    }

    for (size_t i = 0; i < count; ++i) {
        new (ptr + i) float;
    }
    return ptr;
}

template <>
double *allocate(size_t count)
{
    if (count > INT_MAX) {
#ifndef NO_EXCEPTIONS
        throw std::length_error("Size overflow in allocate");
#else
        abort();
#endif
    }
    
    double *ptr = ippsMalloc_64f(int(count));
    if (!ptr) {
#ifndef NO_EXCEPTIONS
        throw (std::bad_alloc());
#else
        abort();
#endif
    }

    for (size_t i = 0; i < count; ++i) {
        new (ptr + i) double;
    }
    return ptr;
}

template <>
void deallocate(float *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

template <>
void deallocate(double *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

#endif

}

