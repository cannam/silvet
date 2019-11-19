/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    bqvec

    A small library for vector arithmetic and allocation in C++ using
    raw C pointer arrays.

    Copyright 2007-2018 Particular Programs Ltd.

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

#include "Barrier.h"

#if defined __APPLE__
#if !defined __MAC_10_12
#include <libkern/OSAtomic.h>
#endif
#endif
#if defined _WIN32 && defined _MSC_VER
#include <Windows.h>
#endif

namespace breakfastquay {

void system_memorybarrier()
{
#if defined __APPLE__

#if defined __MAC_10_12
    __sync_synchronize();
#else
    OSMemoryBarrier();
#endif
    
#elif (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1)

    __sync_synchronize();

#elif defined _WIN32 

#if defined _MSC_VER
    MemoryBarrier();
#else /* (mingw) */
    LONG Barrier = 0;
    __asm__ __volatile__("xchgl %%eax,%0 "
                         : "=r" (Barrier));
#endif

#else
#warning "No memory barrier defined"
#endif
    
}

}

