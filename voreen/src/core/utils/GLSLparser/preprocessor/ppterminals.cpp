/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/GLSLparser/preprocessor/ppterminals.h"

namespace voreen {

namespace glslparser {

PreprocessorTerminals::PreprocessorTerminals() {
    insertSymbol("#", ID_FENCE);
    insertSymbol("#define", ID_DEFINE, true);
    insertSymbol("#undef", ID_UNDEF, true);
    insertSymbol("#if", ID_IF, true);
    insertSymbol("#elif", ID_ELIF, true);
    insertSymbol("#else", ID_ELSE, true);
    insertSymbol("#endif", ID_ENDIF, true);
    insertSymbol("#ifdef", ID_IFDEF, true);
    insertSymbol("#ifndef", ID_IFNDEF, true);
    insertSymbol("#error", ID_ERROR, true);
    insertSymbol("#pragma", ID_PRAGMA, true);
    insertSymbol("#version", ID_VERSION, true);
    insertSymbol("#line", ID_LINE, true);
    insertSymbol("#extension", ID_EXTENSION, true);
    insertSymbol("#include", ID_INCLUDE, true);

    insertSymbol("defined", ID_OP_DEFINED, true);
    insertSymbol("##", ID_OP_FFENCE);

    insertSymbol("(", ID_LPAREN);
    insertSymbol(")", ID_RPAREN);
    insertSymbol(",", ID_COMMA);
    insertSymbol("\"", ID_QUOTE);
    insertSymbol("+", ID_PLUS);
    insertSymbol("-", ID_DASH);
    insertSymbol("~", ID_OP_COMPLEMENT);
    insertSymbol("!", ID_OP_NOT);

    insertSymbol("*", ID_OP_MUL);
    insertSymbol("/", ID_OP_DIV);
    insertSymbol("%", ID_OP_MOD);
    insertSymbol("<<", ID_OP_LSHIFT);
    insertSymbol(">>", ID_OP_RSHIFT);
    insertSymbol("<", ID_LESS);
    insertSymbol("<=", ID_LEQ);
    insertSymbol(">=", ID_GEQ);
    insertSymbol(">", ID_GREATER);
    insertSymbol("==", ID_EQ);
    insertSymbol("!=", ID_NEQ);
    insertSymbol("&", ID_BIT_AND);
    insertSymbol("^", ID_BIT_XOR);
    insertSymbol("|", ID_BIT_OR);
    insertSymbol("&&", ID_LOGICAL_AND);
    insertSymbol("||", ID_LOGICAL_OR);

    insertSymbol("identifier", ID_IDENTIFIER);
    insertSymbol("token", ID_TOKEN);
    insertSymbol("text", ID_TEXT);
    insertSymbol("int-constant", ID_CONSTANT);
    insertSymbol("string", ID_STRING);
    insertSymbol("newline", ID_NEWLINE);
    insertSymbol("macro-evaluation", ID_MACRO);
    insertSymbol("pirate", ID_PIRATE);
}

}   // namespace glslparser

}   // namespace voreen
