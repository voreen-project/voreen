/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#include "tgt/assert.h"

#include <sstream>
#include <iostream>
#include <cstdlib>

#ifdef TGT_DEBUG

#ifndef TGT_NON_INTERACTIVE_ASSERT
bool myCustomAssert(long line, const char* filename, const char* functionName, const char* description, bool& always) {
#else
bool myCustomAssert(long line, const char* filename, const char* functionName, const char* description, bool& /*always*/) {
#endif
    std::cout << "Assertion: " << filename << " " << line << ":" << std::endl << "     " << functionName << std::endl;
    std::cout << "   " << description << std::endl;
#ifndef TGT_NON_INTERACTIVE_ASSERT
    do {
        std::cout << "   (b)reakpoint, (i)gnore, ignore (a)lways, (q)uit" << std::endl;
        std::string inStr;
        std::cin >> inStr;
        switch (inStr[0]) {
            case 'b':
                return true;
            case 'i':
                return false;
            case 'a':
                always = true;
                return false;
            case 'q':
                exit(EXIT_FAILURE);
            default:
                std::cout << "An assertion in a assertion... :(" << std::endl;
        }
    } while (true);
#else
    exit(EXIT_FAILURE);
#endif
}

#endif // TGT_DEBUG
