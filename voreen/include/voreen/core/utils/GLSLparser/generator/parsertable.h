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

#ifndef VRN_PARSERTABLE_H
#define VRN_PARSERTABLE_H

#include "voreen/core/utils/GLSLparser/parseraction.h"

#include <list>

namespace voreen {

namespace glslparser {

class ParserAction;

class ParserTableRow {
public:
    ParserTableRow(const int rowNumber);
    ~ParserTableRow();

    bool addAction(GrammarSymbol* const symbol, ParserAction* const action);
    ParserAction* findAction(GrammarSymbol* const symbol) const;
    int getNumColumns() const { return static_cast<int>(columns_.size()); }
    int getRowNumber() const { return rowNumber_; }

    std::string generateCode(const std::string& symbolVar, const bool gotoTable) const;

    std::string toString() const;

private:
    typedef std::pair<GrammarSymbol*, ParserAction*> RowEntry;
    std::vector<RowEntry> columns_;
    const int rowNumber_;
};

// ============================================================================

class ParserTable {
public:
    ParserTable(const int numStates, const std::vector<GrammarSymbol*>& symbols,
        const std::vector<Production>& productions);
    ~ParserTable();

    bool addAction(const int state, GrammarSymbol* const symbol, ParserAction* const action);
    ParserAction* findAction(const int state, GrammarSymbol* const symbol) const;

    void generateActionCode(std::ostream& oss, const std::string& stateVar,
        const std::string& symbolVar, const std::string& className = "") const;
    void generateGotoCode(std::ostream& oss, const std::string& stateVar,
        const std::string& symbolVar, const std::string& className = "") const;
    void generateProductionsCode(std::ostream& oss, const std::string& prodVar,
        const std::string& className = "") const;
    void generateSymbolsCode(std::ostream& oss, const std::string& symbolVar,
        const std::string& className = "") const;

    void htmlOutput(const std::string& fileName) const;

private:
    static bool ProductionsCompare(const Production& p1, const Production& p2) {
        return (p1.getProductionID() < p2.getProductionID());
    }

    std::string addTabs(const std::string& input, const int numTabs) const;

private:
    const int numStates_;
    ParserTableRow** rows_;
    std::vector<GrammarSymbol*> symbols_;
    std::vector<Production> productions_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
