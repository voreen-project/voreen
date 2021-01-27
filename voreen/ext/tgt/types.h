/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_TYPES_H
#define TGT_TYPES_H

/*
 * This file is a derived work of the following file:
 */

/*
 * stdint.h - integer types for Win32
 * Copyright (C) 2003 ETC s.r.o.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the ETC s.r.o. nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Written by Marcel Telka <marcel@telka.sk>, 2003.
 *
 */

/*
    This file defines the very practical C99 types et al.
    The most useful ones are:

    int8_t
    int16_t
    int32_t
    int64_t
    intptr_t

    uint8_t
    uint16_t
    uint32_t
    uint64_t
    uintptr_t
    uint

    size_t
    ulong

    TGT_UNUSED
*/

/**
 * This is needed for .dll or .so support respectively
 */
#ifdef TGT_BUILD_DLL
    // building library -> export symbols
    #ifdef WIN32
        #define TGT_API __declspec(dllexport)
    #else
        #define TGT_API
    #endif
#else
    // including library -> import symbols
    #ifdef WIN32
        #define TGT_API __declspec(dllimport)
    #else
        #define TGT_API
    #endif
#endif


/**
 * With this macro you can get rid of annoying "unused parameter" warnings
 */
#ifdef __GNUC__
    #define TGT_UNUSED __attribute__ ((unused))
#else
    /// empty macro for unsupported compilers
    #define TGT_UNUSED
#endif

//#ifdef _MSC_VER
//    #pragma warning( disable : 4100 )
//#endif


// For size_t
#include <cstddef>

// Have not found a good place for this
#include <cstdlib>

/// ulong for 32 or 64 bit sizes. This depends on your arch and/or OS
typedef unsigned long ulong;

#if defined(__APPLE__) || (defined(WIN32) && defined(__GNUC__))
    /* unsigned int type */
    typedef unsigned int uint;
#endif

#if defined(_MSC_VER) && (_MSC_VER < 1600)
    // MSVC++ prior to Visual Studio 2010 does not provide stdint.h
    #define _USE_MATH_DEFINES

    #include <windows.h>
    #include <limits.h>

    /*
     * Integer Types
     */

    /* Exact-width integer types */

    #ifndef TGT_USE_SDL_WIN32
    typedef CHAR     int8_t;
    #endif

    typedef SHORT    int16_t;
    typedef INT32    int32_t;
    typedef INT64    int64_t;
    typedef UCHAR   uint8_t;
    typedef WORD    uint16_t;
    typedef UINT32  uint32_t;
    typedef UINT64  uint64_t;

    /* unsigned int type */

    typedef unsigned int uint;

    /* Minimum-width integer types */

    typedef CHAR     int_least8_t;
    typedef SHORT    int_least16_t;
    typedef INT32    int_least32_t;
    typedef INT64    int_least64_t;
    typedef UCHAR   uint_least8_t;
    typedef WORD    uint_least16_t;
    typedef UINT32  uint_least32_t;
    typedef UINT64  uint_least64_t;

    /* Fastest minimum-width integer types */

    typedef CHAR     int_fast8_t;
    typedef SHORT    int_fast16_t;
    typedef INT32    int_fast32_t;
    typedef INT64    int_fast64_t;
    typedef UCHAR   uint_fast8_t;
    typedef WORD    uint_fast16_t;
    typedef UINT32  uint_fast32_t;
    typedef UINT64  uint_fast64_t;

    /* Integer types capable of holding object pointers */

    typedef INT_PTR   intptr_t;
    typedef UINT_PTR uintptr_t;

    /* Greatest-width integer types */

    typedef INT64    intmax_t;
    typedef UINT64  uintmax_t;

    #if !defined(__cplusplus) || defined(__STDC_LIMIT_MACROS)

        /*
         * Limits of Specified-Width Integer Types
         */

        /* Limits of exact-width integer types */

        #define INT8_MIN        (-127 - 1)
        #define INT16_MIN       (-32767 - 1)
        #define INT32_MIN       (-2147483647 - 1)
        #define INT64_MIN       (-9223372036854775807i64 - 1)

        #define INT8_MAX        127
        #define INT16_MAX       32767
        #define INT32_MAX       2147483647
        #define INT64_MAX       9223372036854775807i64

        #define UINT8_MAX       0xFF
        #define UINT16_MAX      0xFFFF
        #define UINT32_MAX      0xFFFFFFFF
        #define UINT64_MAX      0xFFFFFFFFFFFFFFFFui64

        /* Limits of minimum-width integer types */

        #define INT_LEAST8_MIN       INT8_MIN
        #define INT_LEAST16_MIN      INT16_MIN
        #define INT_LEAST32_MIN      INT32_MIN
        #define INT_LEAST64_MIN      INT64_MIN

        #define INT_LEAST8_MAX       INT8_MAX
        #define INT_LEAST16_MAX      INT16_MAX
        #define INT_LEAST32_MAX      INT32_MAX
        #define INT_LEAST64_MAX      INT64_MAX

        #define UINT_LEAST8_MAX     UINT8_MAX
        #define UINT_LEAST16_MAX    UINT16_MAX
        #define UINT_LEAST32_MAX    UINT32_MAX
        #define UINT_LEAST64_MAX    UINT64_MAX

        /* Limits of fastest minimum-width integer types */

        #define INT_FAST8_MIN        INT8_MIN
        #define INT_FAST16_MIN       INT16_MIN
        #define INT_FAST32_MIN       INT32_MIN
        #define INT_FAST64_MIN       INT64_MIN

        #define INT_FAST8_MAX        INT8_MAX
        #define INT_FAST16_MAX       INT16_MAX
        #define INT_FAST32_MAX       INT32_MAX
        #define INT_FAST64_MAX       INT64_MAX

        #define UINT_FAST8_MAX      UINT8_MAX
        #define UINT_FAST16_MAX     UINT16_MAX
        #define UINT_FAST32_MAX     UINT32_MAX
        #define UINT_FAST64_MAX     UINT64_MAX

        /* Limits of integer types capable of holding object pointers */

        #define INTPTR_MIN       INT32_MIN
        #define INTPTR_MAX       INT32_MAX
        #define UINTPTR_MAX     UINT32_MAX

        /* Limits of greatest-width integer types */

        #define INTMAX_MIN       INT64_MIN
        #define INTMAX_MAX       INT64_MAX
        #define UINTMAX_MAX     UINT64_MAX

        /*
         * Limits of Other Integer Types
         */

        /* Limits of ptrdiff_t */

        #define PTRDIFF_MIN      INT_MIN
        #define PTRDIFF_MAX      INT_MAX

        /* Limits of sig_atomic_t */

        /* N/A for Windows */

        /* Limit of size_t */

        #define SIZE_MAX        UINT_MAX

        /* Limits of wchar_t */

        #define WCHAR_MIN       0
        #define WCHAR_MAX       USHRT_MAX

        /* Limits of wint_t */

        #define WINT_MIN        WCHAR_MIN
        #define WINT_MAX        WCHAR_MAX

    #endif /* __STDC_LIMIT_MACROS */

    #if !defined(__cplusplus) || defined(__STDC_CONSTANT_MACROS)

        /*
         * Macros for Integer Constant Expressions
         */

        /* Macros for minimum-width integer constant expressions */

        #define INT8_C(value)       (value)
        #define INT16_C(value)      (value)
        #define INT32_C(value)      (value)
        #define INT64_C(value)      (value##i64)

        #define UINT8_C(value)      (value)
        #define UINT16_C(value)     (value)
        #define UINT32_C(value)     (value)
        #define UINT64_C(value)     (value##ui64)

        /* Macros for greatest-width integer constant expressions */

        #define INTMAX_C(value)      INT64_C(value)
        #define UINTMAX_C(value)    UINT64_C(value)

    #endif /* __STDC_CONSTANT_MACROS */

#else // So it must be UNIX or VS2010 (or later)

    /*
        The ISO C99 standard specifies that in C++ implementations limit
        macros should only be defined if explicitly requested.
    */
    //FIXME: so this always defines the macros. the application should set
    // this define instead only if needed.
    #define __STDC_LIMIT_MACROS
    #include <stdint.h>
    // For uint - hope this works on all *NIX variants
    #include <sys/types.h>

#endif //WIN32

#ifdef WIN32
    #define TGT_PATH_SEPARATOR "\\"
#else
    #define TGT_PATH_SEPARATOR "/"
#endif


#endif //TGT_TYPES_H
