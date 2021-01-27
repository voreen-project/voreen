/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
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

#ifndef VRN_GLSLTERMINALS_H
#define VRN_GLSLTERMINALS_H

#include "voreen/core/utils/GLSLparser/terminalsmap.h"

#include <set>

namespace voreen {

namespace glslparser {

class GLSLTerminals : public TerminalsMap {
public:
     enum SymbolIDs {
        ID_START = 0,

        ID_LSHIFT_OP, ID_RSHIFT_OP, ID_INC_OP, ID_DEC_OP,
        ID_LE_OP, ID_GE_OP, ID_EQ_OP, ID_NE_OP,
        ID_AND_OP, ID_OR_OP, ID_XOR_OP,

        ID_LSHIFT_ASSIGN, ID_RSHIFT_ASSIGN,
        ID_ADD_ASSIGN, ID_DIV_ASSIGN, ID_MUL_ASSIGN, ID_SUB_ASSIGN,
        ID_MOD_ASSIGN, ID_AND_ASSIGN, ID_OR_ASSIGN, ID_XOR_ASSIGN,

        ID_LPAREN, ID_RPAREN, ID_LBRACKET, ID_RBRACKET,
        ID_LBRACE, ID_RBRACE, ID_LANGLE, ID_RANGLE,
        ID_DOT, ID_COMMA, ID_COLON, ID_SEMICOLON, ID_EQUAL, ID_BANG, ID_DASH, ID_TILDE,
        ID_PLUS, ID_STAR, ID_SLASH, ID_PERCENT, ID_VERTICAL_BAR, ID_CARET,
        ID_AMPERSAND, ID_QUESTION,

        ID_ATTRIBUTE, ID_CONST, ID_IN, ID_INOUT,
        ID_OUT, ID_UNIFORM, ID_VARYING, ID_CENTROID,
        ID_INVARIANT,
        ID_LAYOUT, ID_FLAT, ID_SMOOTH, ID_NONPERSPECTIVE,
        ID_LOWP, ID_HIGHP, ID_MEDIUMP, ID_PRECISION,

        ID_DO, ID_ELSE, ID_FOR, ID_IF, ID_WHILE,
        ID_CASE, ID_DEFAULT, ID_SWITCH,
        ID_BREAK, ID_CONTINUE, ID_DISCARD, ID_RETURN,
        ID_STRUCT, ID_FALSE, ID_TRUE, ID_VOID,

        // DO NOT INSERT ANY SYMBOLS BELOW THIS COMMENT !!!
        // THE ORDER IS USED FOR SIMPLIFYING THE PARSER GENERATOR SETUP!
        //
        ID_BOOL, ID_FLOAT, ID_INT, ID_UINT,

        ID_BVEC2, ID_BVEC3, ID_BVEC4,
        ID_IVEC2, ID_IVEC3, ID_IVEC4,
        ID_UVEC2, ID_UVEC3, ID_UVEC4,
        ID_VEC2, ID_VEC3, ID_VEC4,

        ID_MAT2, ID_MAT3, ID_MAT4,
        ID_MAT2X2, ID_MAT2X3, ID_MAT2X4,
        ID_MAT3X2, ID_MAT3X3, ID_MAT3X4,
        ID_MAT4X2, ID_MAT4X3, ID_MAT4X4,

        ID_SAMPLER1D, ID_SAMPLER2D, ID_SAMPLER3D,
        ID_SAMPLERCUBE, ID_SAMPLER1DSHADOW, ID_SAMPLER2DSHADOW,
        ID_SAMPLER1DARRAY, ID_SAMPLER2DARRAY,
        ID_SAMPLER1DARRAYSHADOW, ID_SAMPLER2DARRAYSHADOW,
        ID_ISAMPLER1D, ID_ISAMPLER2D, ID_ISAMPLER3D, ID_ISAMPLERCUBE,
        ID_ISAMPLER1DARRAY, ID_ISAMPLER2DARRAY,
        ID_USAMPLER1D, ID_USAMPLER2D, ID_USAMPLER3D, ID_USAMPLERCUBE,
        ID_USAMPLER1DARRAY, ID_USAMPLER2DARRAY,
        ID_SAMPLER2DRECT, ID_SAMPLER2DRECTSHADOW,
        ID_ISAMPLER2DRECT, ID_USAMPLER2DRECT,
        ID_SAMPLERBUFFER, ID_ISAMPLERBUFFER, ID_USAMPLEBUFFER,
        ID_SAMPLER2DMS, ID_ISAMPLER2DMS, ID_USAMPLER2DMS,
        ID_SAMPLER2DMSARRY, ID_ISAMPLER2DMSARRAY, ID_USAMPLER2DMSARRAY,
        //
        // DO NOT INSERT ANY SYMBOLS ABOVE THIS COMMENT !!!
        // THE ORDER IS USED FOR SIMPLIFYING THE PARSER GENERATOR SETUP!

        ID_IDENTIFIER, ID_TYPE_NAME, ID_FLOATCONST, ID_INTCONST, ID_UINTCONST,
        ID_FIELDSELECTION, ID_NATIVE_TYPE, ID_LEADING_ANNOTATION, ID_TRAILING_ANNOTATION,
        ID_RESERVED_KEYWORD, ID_UNKNOWN,

        ID_END, ID_EMPTY
    };

public:
    GLSLTerminals();

    static std::set<std::string>& getPredefinedStructs();
};

}   // namespace glslparser

}   // namespace voreen

#endif
