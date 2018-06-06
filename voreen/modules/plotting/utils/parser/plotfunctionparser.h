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

#ifndef VRN_PLOTFUNCTIONPARSER_H
#define VRN_PLOTFUNCTIONPARSER_H

#include "voreen/core/utils/GLSLparser/parser.h"

#include <list>

namespace voreen {

namespace glslparser {

class PlotFunctionLexer;

class PlotFunctionParser : public Parser {
public:
    PlotFunctionParser(PlotFunctionLexer* const lexer);

    PlotFunctionParser(std::istream* const is);

    PlotFunctionParser(const std::string& fileName);

    PlotFunctionParser(const std::vector<Token*>& tokens);

    virtual ~PlotFunctionParser();

protected:
    virtual Token* nextToken();

    virtual void expandParseTree(const int productionID,
        const std::vector<Parser::ParserSymbol*>& reductionBody);

    // The following methods are from generated code.

    virtual ParserAction* action(const int stateID, const int SymbolID) const;

    virtual ProductionStub* findProduction(const int productionID) const;

    virtual int gotoState(const int stateID, const int symbolID) const;

    virtual std::string symbolID2String(const int symbolID) const;

private:

    std::list<Token*> tokens_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
