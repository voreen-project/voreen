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

#include "voreen/core/utils/GLSLparser/glsl/glslterminals.h"

namespace voreen {

namespace glslparser {

GLSLTerminals::GLSLTerminals() {
    insertSymbol("<<", ID_LSHIFT_OP);
    insertSymbol(">>", ID_RSHIFT_OP);
    insertSymbol("++", ID_INC_OP);
    insertSymbol("--", ID_DEC_OP);
    insertSymbol("<=", ID_LE_OP);
    insertSymbol(">=", ID_GE_OP);
    insertSymbol("==", ID_EQ_OP);
    insertSymbol("!=", ID_NE_OP);
    insertSymbol("&&", ID_AND_OP);
    insertSymbol("||", ID_OR_OP);
    insertSymbol("^^", ID_XOR_OP);

    insertSymbol("<<=", ID_LSHIFT_ASSIGN);
    insertSymbol(">>=", ID_RSHIFT_ASSIGN);
    insertSymbol("+=", ID_ADD_ASSIGN);
    insertSymbol("/=", ID_DIV_ASSIGN);
    insertSymbol("*=", ID_MUL_ASSIGN);
    insertSymbol("-=", ID_SUB_ASSIGN);
    insertSymbol("%=", ID_MOD_ASSIGN);
    insertSymbol("&=", ID_AND_ASSIGN);
    insertSymbol("|=", ID_OR_ASSIGN);
    insertSymbol("^=", ID_XOR_ASSIGN);

    insertSymbol("(", ID_LPAREN);
    insertSymbol(")", ID_RPAREN);
    insertSymbol("[", ID_LBRACKET);
    insertSymbol("]", ID_RBRACKET);
    insertSymbol("{", ID_LBRACE);
    insertSymbol("}", ID_RBRACE);
    insertSymbol("<", ID_LANGLE);
    insertSymbol(">", ID_RANGLE);
    insertSymbol(".", ID_DOT);
    insertSymbol(",", ID_COMMA);
    insertSymbol(":", ID_COLON);
    insertSymbol(";", ID_SEMICOLON);
    insertSymbol("=", ID_EQUAL);
    insertSymbol("!", ID_BANG);
    insertSymbol("-", ID_DASH);
    insertSymbol("~", ID_TILDE);
    insertSymbol("+", ID_PLUS);
    insertSymbol("*", ID_STAR);
    insertSymbol("/", ID_SLASH);
    insertSymbol("%", ID_PERCENT);
    insertSymbol("|", ID_VERTICAL_BAR);
    insertSymbol("^", ID_CARET);
    insertSymbol("&", ID_AMPERSAND);
    insertSymbol("?", ID_QUESTION);

    insertSymbol("attribute", ID_ATTRIBUTE, true);
    insertSymbol("const", ID_CONST, true);
    insertSymbol("in", ID_IN, true);
    insertSymbol("inout", ID_INOUT, true);
    insertSymbol("out", ID_OUT, true);
    insertSymbol("uniform", ID_UNIFORM, true);
    insertSymbol("varying", ID_VARYING, true);
    insertSymbol("centroid", ID_CENTROID, true);
    insertSymbol("invariant", ID_INVARIANT, true);

    insertSymbol("flat", ID_FLAT, true);
    insertSymbol("layout", ID_LAYOUT, true);
    insertSymbol("noperspective", ID_NONPERSPECTIVE, true);
    insertSymbol("smooth", ID_SMOOTH, true);

    // precision qualifiers
    //
    insertSymbol("highp", ID_HIGHP, true);
    insertSymbol("lowp", ID_LOWP, true);
    insertSymbol("mediump", ID_MEDIUMP, true);
    insertSymbol("precision", ID_PRECISION, true);

    // program control / statements
    //
    insertSymbol("case", ID_CASE, true);
    insertSymbol("default", ID_DEFAULT, true);
    insertSymbol("do", ID_DO, true);
    insertSymbol("else", ID_ELSE, true);
    insertSymbol("for", ID_FOR, true);
    insertSymbol("if", ID_IF, true);
    insertSymbol("switch", ID_SWITCH, true);
    insertSymbol("while", ID_WHILE, true);

    insertSymbol("break", ID_BREAK, true);
    insertSymbol("continue", ID_CONTINUE, true);
    insertSymbol("discard", ID_DISCARD, true);
    insertSymbol("return", ID_RETURN, true);
    insertSymbol("false", ID_FALSE, true);
    insertSymbol("true", ID_TRUE, true);
    insertSymbol("struct", ID_STRUCT, true);

    // "primitive" types
    //
    insertSymbol("bool", ID_BOOL, true);
    insertSymbol("float", ID_FLOAT, true);
    insertSymbol("int", ID_INT, true);
    insertSymbol("uint", ID_UINT, true);
    insertSymbol("void", ID_VOID, true);

    // "complexe" types
    //
    insertSymbol("bvec2", ID_BVEC2, true);
    insertSymbol("bvec3", ID_BVEC3, true);
    insertSymbol("bvec4", ID_BVEC4, true);

    insertSymbol("ivec2", ID_IVEC2, true);
    insertSymbol("ivec3", ID_IVEC3, true);
    insertSymbol("ivec4", ID_IVEC4, true);

    insertSymbol("uvec2", ID_UVEC2, true);
    insertSymbol("uvec3", ID_UVEC3, true);
    insertSymbol("uvec4", ID_UVEC4, true);

    insertSymbol("vec2", ID_VEC2, true);
    insertSymbol("vec3", ID_VEC3, true);
    insertSymbol("vec4", ID_VEC4, true);

    insertSymbol("mat2", ID_MAT2, true);
    insertSymbol("mat3", ID_MAT3, true);
    insertSymbol("mat4", ID_MAT4, true);

    insertSymbol("mat2x2", ID_MAT2X2, true);
    insertSymbol("mat2x3", ID_MAT2X3, true);
    insertSymbol("mat2x4", ID_MAT2X4, true);

    insertSymbol("mat3x2", ID_MAT3X2, true);
    insertSymbol("mat3x3", ID_MAT3X3, true);
    insertSymbol("mat3x4", ID_MAT3X4, true);

    insertSymbol("mat4x2", ID_MAT4X2, true);
    insertSymbol("mat4x3", ID_MAT4X3, true);
    insertSymbol("mat4x4", ID_MAT4X4, true);

    insertSymbol("sampler1D", ID_SAMPLER1D, true);
    insertSymbol("sampler2D", ID_SAMPLER2D, true);
    insertSymbol("sampler3D", ID_SAMPLER3D, true);
    insertSymbol("samplerCube", ID_SAMPLERCUBE, true);
    insertSymbol("sampler1DShadow", ID_SAMPLER1DSHADOW, true);
    insertSymbol("sampler2DShadow", ID_SAMPLER2DSHADOW, true);

    insertSymbol("sampler1DArray", ID_SAMPLER1DARRAY, true);
    insertSymbol("sampler2DArray", ID_SAMPLER2DARRAY, true);
    insertSymbol("sampler1DArrayShadow", ID_SAMPLER1DARRAYSHADOW, true);
    insertSymbol("sampler2DArrayShadow", ID_SAMPLER2DARRAYSHADOW, true);
    insertSymbol("isampler1D", ID_ISAMPLER1D, true);
    insertSymbol("isampler2D", ID_ISAMPLER2D, true);
    insertSymbol("isampler3D", ID_ISAMPLER3D, true);
    insertSymbol("isamplerCube", ID_ISAMPLERCUBE, true);
    insertSymbol("isampler1DArray", ID_ISAMPLER1DARRAY, true);
    insertSymbol("isampler2DArray", ID_ISAMPLER2DARRAY, true);
    insertSymbol("usampler1D", ID_USAMPLER1D, true);
    insertSymbol("usampler2D", ID_USAMPLER2D, true);
    insertSymbol("usampler3D", ID_USAMPLER3D, true);
    insertSymbol("usamplerCube", ID_USAMPLERCUBE, true);
    insertSymbol("usampler1DArray", ID_USAMPLER1DARRAY, true);
    insertSymbol("usampler2DArray", ID_USAMPLER2DARRAY, true);
    insertSymbol("sampler2DRect", ID_SAMPLER2DRECT, true);
    insertSymbol("sampler2DRectShadow", ID_SAMPLER2DRECTSHADOW, true);
    insertSymbol("isampler2DRect", ID_ISAMPLER2DRECT, true);
    insertSymbol("usampler2DRect", ID_USAMPLER2DRECT, true);
    insertSymbol("samplerBuffer", ID_SAMPLERBUFFER, true);
    insertSymbol("isamplerBuffer", ID_ISAMPLERBUFFER, true);
    insertSymbol("usamplerBuffer", ID_USAMPLEBUFFER, true);
    insertSymbol("sampler2DMS", ID_SAMPLER2DMS, true);
    insertSymbol("isampler2DMS", ID_ISAMPLER2DMS, true);
    insertSymbol("usampler2DMS", ID_USAMPLER2DMS, true);
    insertSymbol("sampler2DMSArray", ID_SAMPLER2DMSARRY, true);
    insertSymbol("isampler2DMSArray", ID_ISAMPLER2DMSARRAY, true);
    insertSymbol("usampler2DMSArray", ID_USAMPLER2DMSARRAY, true);

    insertSymbol("IDENTIFIER", ID_IDENTIFIER);
    insertSymbol("INTCONSTANT", ID_INTCONST);
    insertSymbol("UINTCONSTANT", ID_UINTCONST);
    insertSymbol("FLOATCONSTANT", ID_FLOATCONST);
    insertSymbol("TYPE-NAME", ID_TYPE_NAME);
    insertSymbol("FIELD-SELECTION", ID_FIELDSELECTION);
    insertSymbol("NATIVE-TYPE", ID_NATIVE_TYPE);
    insertSymbol("LEADING-ANNOTATION", ID_LEADING_ANNOTATION);
    insertSymbol("TRAILING-ANNOTATION", ID_TRAILING_ANNOTATION);

    // other keywords reserved by GLSL specification 1.50
    //
    const size_t NUM_OTHERS = 70;
    std::string others[NUM_OTHERS] = {
        "common", "partition", "active", "asm", "class",
        "union", "enum", "typedef", "template", "this",
        "packed", "goto", "inline", "noinline", "volatile",
        "public", "static", "extern", "external", "interface",
        "long", "short", "double", "half", "fixed",
        "unsigned", "superp", "input", "output", "hvec2",
        "hvec3", "hvec4", "dvec2", "dvec3", "dvec4",
        "fvec2", "fvec3", "fvec4", "sampler3DRect", "filter",
        "image1D", "image2D", "image3D", "imageCube", "iimage1D",
        "iimage2D", "iimage3D", "iimageCube", "uimage1D", "uimage2D",
        "uimage3D", "uimageCube", "image1DArray", "image2DArray", "iimage1DArray",
        "iimage2DArray", "uimage1DArray", "uimage2DArray", "image1DShadow", "image2DShadow",
        "image1DArrayShadow", "image2DArrayShadow", "imageBuffer", "iimageBuffer", "uimageBuffer",
        "sizeof", "cast", "namespace", "using", "row_major"
    };

    for (size_t i = 0; i < NUM_OTHERS; ++i)
        insertSymbol(others[i], ID_RESERVED_KEYWORD, true);

}

std::set<std::string>& GLSLTerminals::getPredefinedStructs() {
    static std::set<std::string> predefinedStructs_;
    if (predefinedStructs_.empty()) {
        // Add the names of predefined structs into a set in order to correctly identify them as
        // TYPE-NAME tokens during lexing.
        //
        predefinedStructs_.insert("gl_FogParameters");
        predefinedStructs_.insert("gl_LightModelParameters");
        predefinedStructs_.insert("gl_LightModelProducts");
        predefinedStructs_.insert("gl_LightProducts");
        predefinedStructs_.insert("gl_LightSourceParameters");
        predefinedStructs_.insert("gl_MaterialParameters");
        predefinedStructs_.insert("gl_PointParameters");
    }

    return predefinedStructs_;
}

}   // namespace glslparser

}   // namespace voreen
