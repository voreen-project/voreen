/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
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

#ifndef TGT_ASSERT_H
#define TGT_ASSERT_H

#include "tgt/types.h"

#ifdef TGT_DEBUG

    /**
     * This is the assert function. It is used by the tgtAssert macro.
     * @param line          line number
     * @param filename      file name
     * @param functionName  function name
     * @param description   description of the assertion
     * @param always        should this assertion be ignored always? This is a kind of a return value
     * @return              true if a breakpoint should be thrown
    */
    TGT_API bool myCustomAssert(long line, const char* filename, const char* functionName,
                        const char* description, bool& always);

    /*
     * THROW_BREAKPOINT activates a Breakpoint in debug mode.
     * Does it work with PowerPC? I have no idea. I have googled the asm instruction...
    */

    //__GNUC__ is defined by all GCC versions
    #ifdef __GNUC__
        #ifdef __APPLE__
             //#define TGT_THROW_BREAKPOINT asm("trap")
             //FIXME: need to check architecture here!
             #define TGT_THROW_BREAKPOINT {}
        #else //Then it is a Linux system or Mingw
            #define TGT_THROW_BREAKPOINT asm("int $3")
        #endif //__APPLE__
    #elif defined(WIN32)
        #define TGT_THROW_BREAKPOINT __debugbreak();
    #else //Ok, so it is an unsupported arch -> no breakpoints :(
        #define TGT_THROW_BREAKPOINT {}
    #endif

    /**
     * the assertion macro
     * @param e boolean expression. This one should be true
     * @param description a c-string. The description for the assertion.
    */
    #define tgtAssert(e, description) \
        if ( !(bool(e)) ) { \
            static bool always = false; \
            if (!always) \
                if (myCustomAssert(__LINE__, __FILE__, __FUNCTION__, (description), always)) \
                    TGT_THROW_BREAKPOINT; \
        }

#else //TGT_DEBUG

    /**
     * No functionality in the release version
    */
    #define tgtAssert(b, description)

#endif //TGT_DEBUG


#endif //TGT_ASSERT_H
