/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#ifndef VRN_PPTERMINALS_H
#define VRN_PPTERMINALS_H

#include "voreen/core/utils/GLSLparser/terminalsmap.h"

namespace voreen {

namespace glslparser {

class PreprocessorTerminals : public TerminalsMap {
public:
    enum SymbolIDs {
        ID_START = 0, ID_FENCE, ID_DEFINE, ID_UNDEF,
        ID_IF, ID_ELIF, ID_ELSE, ID_ENDIF, ID_IFDEF, ID_IFNDEF,
        ID_ERROR, ID_PRAGMA, ID_VERSION, ID_LINE, ID_EXTENSION,
        ID_INCLUDE,

        ID_OP_DEFINED, ID_OP_FFENCE, // operator ##

        ID_LPAREN, ID_RPAREN, ID_COMMA, ID_QUOTE,
        ID_PLUS, ID_DASH, ID_OP_COMPLEMENT, ID_OP_NOT,
        ID_OP_MUL, ID_OP_DIV, ID_OP_MOD,
        ID_OP_LSHIFT, ID_OP_RSHIFT,
        ID_LESS, ID_LEQ, ID_GEQ, ID_GREATER,
        ID_EQ, ID_NEQ,
        ID_BIT_AND, ID_BIT_XOR, ID_BIT_OR,
        ID_LOGICAL_AND, ID_LOGICAL_OR,

        ID_TOKEN, ID_TEXT, ID_IDENTIFIER, ID_CONSTANT, ID_STRING,
        ID_NEWLINE, ID_MACRO, ID_UNKNOWN, ID_PIRATE, ID_END, ID_EMPTY
    };

public:
    PreprocessorTerminals();
};

}   // namespace glslparser

}   // namespace

#endif
