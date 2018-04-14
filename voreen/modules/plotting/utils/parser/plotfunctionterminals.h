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

#ifndef VRN_ANNOTATION_TERMINALS_H
#define VRN_ANNOTATION_TERMINALS_H

#include "voreen/core/utils/GLSLparser/terminalsmap.h"

#include <vector>
#include <string>

namespace voreen {

namespace glslparser {

class PlotFunctionTerminals : public TerminalsMap {
public:
    /// possible symbolids in this grammar
    enum SymbolIDs {
        ID_START = 0, ID_LPAREN, ID_RPAREN, ID_LBRACKET, ID_RBRACKET, ID_DOT, ID_COMMA, ID_COLON, ID_SEMICOLON,
        ID_DASH, ID_PLUS, ID_STAR, ID_SLASH, ID_CARET,
        ID_VERTICAL_BAR,
        ID_RANGE_TERM, ID_FUNCTION_TERM, ID_FUNCTION,
        ID_INTCONST, ID_FLOATCONST, ID_VARIABLE,
        ID_UNKNOWN,
        ID_EMPTY,
        ID_END,
    };

    PlotFunctionTerminals();

    virtual ~PlotFunctionTerminals() {};
    // gives back the vector of allowed functions, the can interpret
    static const std::vector<std::string>& getAllowedFunctionList();

private:

    static std::vector<std::string> setFunctionList();

    static const std::vector<std::string> functionList_;

};

}   // namespace glslparser

}   // namespace voreen

#endif
