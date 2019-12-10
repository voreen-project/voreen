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

#ifndef VRN_GLSLPARSER_H
#define VRN_GLSLPARSER_H

#include "voreen/core/utils/GLSLparser/parser.h"

namespace voreen {

namespace glslparser {

class GLSLParser : public Parser {
public:
    GLSLParser(std::istream* const is);
    GLSLParser(const std::string& fileName);

    virtual ~GLSLParser();

protected:
    /**
     * Nodes-stack operations to perform when reducing the production related
     * to the given ID. This actually builds the ParseTreeNodes which hold the
     * semantics of the parsed input.
     */
    virtual void expandParseTree(const int productionID,
        const std::vector<Parser::ParserSymbol*>& reductionBody);

    virtual bool isProxySymbol(const int symbolID, const int originalID) const;

    // The following methods are from generated code.

    virtual ParserAction* action(const int state, const int realSymbolID) const;

    /**
     * Returns a ProductionStub
     */
    virtual ProductionStub* findProduction(const int productionID) const;

    /**
     * Goto-Table for this parser.
     */
    virtual int gotoState(const int state, const int symbol) const;

    /**
     * Returns the corresponding lexeme for the symbol with the given ID
     * if such a symbol exists. Otherwise the returned string is empty.
     */
    virtual std::string symbolID2String(const int symbolID) const;
};

}   // namespace glslparser

}   // namespace voreen

#endif // VRN_GLSLPARSER_H
