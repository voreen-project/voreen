/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "plotfunctionterminals.h"

namespace voreen {

namespace glslparser {

const std::vector<std::string> PlotFunctionTerminals::functionList_ = PlotFunctionTerminals::setFunctionList();

PlotFunctionTerminals::PlotFunctionTerminals() {
    insertSymbol("(", ID_LPAREN);
    insertSymbol(")", ID_RPAREN);
    insertSymbol("[", ID_LBRACKET);
    insertSymbol("]", ID_RBRACKET);
    insertSymbol(".", ID_DOT);
    insertSymbol(",", ID_COMMA);
    insertSymbol(":", ID_COLON);
    insertSymbol(";", ID_SEMICOLON);
    insertSymbol("-", ID_DASH);
    insertSymbol("+", ID_PLUS);
    insertSymbol("*", ID_STAR);
    insertSymbol("/", ID_SLASH);
    insertSymbol("^", ID_CARET);
    insertSymbol("|", ID_VERTICAL_BAR);

    insertSymbol("RANGE_TERM", ID_RANGE_TERM);
    insertSymbol("FUNCTION-TERM",ID_FUNCTION_TERM);
    insertSymbol("FUNCTION", ID_FUNCTION);
    insertSymbol("INTCONST", ID_INTCONST);
    insertSymbol("FLOATCONST", ID_FLOATCONST);
    insertSymbol("VARIABLE", ID_VARIABLE);
    insertSymbol("EMPTY",ID_EMPTY);
}

const std::vector<std::string>& PlotFunctionTerminals::getAllowedFunctionList() {
    return functionList_;
}

std::vector<std::string> PlotFunctionTerminals::setFunctionList() {
    std::vector<std::string> functionList;
    functionList.resize(21);
    functionList[0] = "abs";
    functionList[1] = "sqrt";
    functionList[2] = "sin";
    functionList[3] = "cos";
    functionList[4] = "tan";
    functionList[5] = "arcsin";
    functionList[6] = "arccos";
    functionList[7] = "arctan";
    functionList[8] = "sinh";
    functionList[9] = "cosh";
    functionList[10] = "tanh";
    functionList[11] = "ln";
    functionList[12] = "log";
    functionList[13] = "exp";
    functionList[14] = "fac";
    functionList[15] = "int";
    functionList[16] = "floor";
    functionList[17] = "ceil";
    functionList[18] = "rnd";
    functionList[19] = "sgn";
    functionList[20] = "sgx";

    return functionList;
}


}   // namespace glslparser

}   // namespace voreen
