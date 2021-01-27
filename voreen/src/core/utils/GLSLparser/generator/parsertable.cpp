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

#include "voreen/core/utils/GLSLparser/generator/parsertable.h"

#include "voreen/core/utils/GLSLparser/grammarsymbol.h"
#include "voreen/core/utils/GLSLparser/parseraction.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <new>
#include <sstream>

namespace voreen {

namespace glslparser {

ParserTableRow::ParserTableRow(const int rowNumber)
    : rowNumber_(rowNumber)
{
}

ParserTableRow::~ParserTableRow() {
    for (size_t i = 0; i < columns_.size(); ++i)
        delete columns_[i].second;
}

bool ParserTableRow::addAction(GrammarSymbol* const symbol, ParserAction* const action) {
    if ((symbol == 0) || (action == 0))
        return false;

    if (findAction(symbol) != 0)
        return false;

    columns_.push_back(RowEntry(symbol, action));
    return true;
}

ParserAction* ParserTableRow::findAction(GrammarSymbol* const symbol) const {
    if (symbol == 0)
        return 0;

    ParserAction* action = 0;
    for (size_t i = 0; ((action == 0) && (i < columns_.size())); ++i)
    {
        if (*(columns_[i].first) == *symbol)
            action = columns_[i].second;
    }

    return action;
}

std::string ParserTableRow::generateCode(const std::string& symbolVar, const bool gotoTable) const {
    if (columns_.empty())
        return "";

    // Group all actions which are equal for different symbols using their
    // string representation and a map
    //
    typedef std::map<std::string, std::list<GrammarSymbol*> > CasesMap;
    CasesMap cases;

    for (size_t i = 0; i < columns_.size(); ++i)
    {
        if (columns_[i].first->isTerminal() == gotoTable)
            continue;

        std::list<GrammarSymbol*> symbols;
        symbols.push_back(columns_[i].first);

        std::ostringstream oss;
        if (! gotoTable)
            oss << "new " << columns_[i].second->generateCode();
        else {
            ParserActionTransition* trans = dynamic_cast<ParserActionTransition*>(columns_[i].second);
            if (trans != 0)
                oss << trans->getNextState() << ";";
        }

        // try to insert the code for the action into the map or add the symbol to the list,
        // if that action is already present.
        //
        std::pair<CasesMap::iterator, bool> res = cases.insert(std::make_pair(oss.str(), symbols));
        if (res.second == false)
            res.first->second.push_back(symbols.back());
    }

    if (cases.empty() == true)
        return "";

    std::ostringstream oss;
    oss << "switch (" << symbolVar << ") {\n";

    // Generate a case label for each symbol occuring in this row and an return-statement
    // for the code constructing the new corresponding action object. This object will be
    // used by the parser.
    //
    for (CasesMap::iterator it = cases.begin(); it != cases.end(); ++it) {
        const std::list<GrammarSymbol*>& symbols = it->second;
        for (std::list<GrammarSymbol*>::const_iterator itSym = symbols.begin();
            itSym != symbols.end(); ++itSym)
        {
            oss << "    case " << (*itSym)->getSymbolID() << ":  ";
            oss << "// " << (*itSym)->toString() << "\n";
        }

        oss << "        return " << it->first << "\n\n";
    }
    oss << "    default:\n";
    if (! gotoTable)
        oss << "        return new ParserActionError();\n";
    else
        oss << "        return -1;\n";
    oss << "}   // switch (" << symbolVar << ")\n";

    return oss.str();
}

std::string ParserTableRow::toString() const {
    std::ostringstream oss;
    for (size_t i = 0; i < columns_.size(); ++i)
        oss << columns_[i].first->toString() << ": " << columns_[i].second->toString() << " | ";

    return oss.str();
}

// ============================================================================

ParserTable::ParserTable(const int numStates, const std::vector<GrammarSymbol*>& symbols,
                         const std::vector<Production>& productions)
    : numStates_(numStates),
    rows_(0),
    symbols_(symbols),
    productions_(productions)
{
    if (numStates_ <= 0)
        throw std::bad_alloc();

    if (productions_.empty() == true)
        throw std::bad_alloc();

    rows_ = new ParserTableRow*[numStates_];
    for (int i = 0; i < numStates_; ++i)
        rows_[i] = new ParserTableRow(i);

    std::sort(productions_.begin(), productions_.end(), ProductionsCompare);
}

ParserTable::~ParserTable() {
    for (int i = 0; i < numStates_; ++i)
        delete rows_[i];
    delete [] rows_;
}

bool ParserTable::addAction(const int state, GrammarSymbol* const symbol, ParserAction* const action) {
    if ((state >= 0) && (state < numStates_))
        return rows_[state]->addAction(symbol, action);

    return false;
}

ParserAction* ParserTable::findAction(const int state, GrammarSymbol* const symbol) const {
    if ((state >= 0) && (state < numStates_))
        return rows_[state]->findAction(symbol);

    return 0;
}

void ParserTable::generateActionCode(std::ostream& oss, const std::string& stateVar,
                                     const std::string& symbolVar, const std::string& className) const
{
    oss << "ParserAction* ";
    if (className.size() > 0)
        oss << className << "::";
    oss << "action(const int " << stateVar << ", const int " << symbolVar <<") const {\n";
    oss << "    switch (" << stateVar << ") {\n";

    for (int i = 0; i < numStates_; ++i) {
        std::string caseString = rows_[i]->generateCode(symbolVar, false);
        if (! caseString.empty()) {
            oss << "        case " << i << ":\n";
            oss << addTabs(caseString, 3);
            oss << "break;\n\n";
        }

        /*
        oss << "        case " << i << ":\n";
        if (caseString.empty() == true)
            oss << "            return new ParserActionError();\n\n";
        else {
            oss << addTabs(caseString, 3);
            oss << "break;\n\n";
        }
        */
    }

    oss << "        default:\n";
    oss << "            return new ParserActionError();\n\n";
    //oss << "            break;\n";
    oss << "    } // switch (" << stateVar << ")\n\n";
    oss << "    return 0;\n";
    oss << "}\n";
}

void ParserTable::generateGotoCode(std::ostream& oss, const std::string& stateVar,
                                   const std::string& symbolVar, const std::string& className) const
{
    oss << "int ";
    if (className.size() > 0)
        oss << className << "::";
    oss << "gotoState(const int " << stateVar << ", const int " << symbolVar <<") const {\n";
    oss << "    switch (" << stateVar << ") {\n";

    for (int i = 0; i < numStates_; ++i) {
        std::string caseString = rows_[i]->generateCode(symbolVar, true);
        if (! caseString.empty()) {
            oss << "        case " << i << ":\n";
            oss << addTabs(caseString, 3);
            oss << "break;\n\n";
        }
        /*
        oss << "        case " << i << ":\n";
        if (caseString.empty() == true)
            oss << "            return -1;\n\n";
        else {
            oss << addTabs(caseString, 3);
            oss << "break;\n\n";
        }
        */
    }

    //oss << "        default:\n";
    //oss << "            break;\n";
    oss << "    } // switch (" << stateVar << ")\n\n";
    //oss << "    return 0;\n";
    oss << "    return -1;\n";
    oss << "}\n";
}

void ParserTable::generateProductionsCode(std::ostream& oss, const std::string& prodVar,
                                          const std::string& className) const {
    oss << "ProductionStub* ";
    if (className.size() > 0)
        oss << className << "::";
    oss << "findProduction(const int " << prodVar << ") const {\n";
    oss << "    switch (" << prodVar << ") {\n";

    for (size_t i = 0; i < productions_.size(); ++i) {
        std::vector<int> bodyIDs = productions_[i].getBodyIDs();

        oss << "        case " << productions_[i].getProductionID() << ":  ";
        oss << "// " << productions_[i].toString() << "\n";
        oss << "            {\n";
        oss << "                int bodyIDs[] = {";

        for (size_t j = 0; j < bodyIDs.size(); ++j) {
            if (j > 0)
                oss << ", ";
            oss << bodyIDs[j];
        }

        oss << "};\n";
        oss << "                return new ProductionStub(" << productions_[i].getHeadID() << ",";
        oss << " bodyIDs, " << productions_[i].getBody().size() << ");\n";
        oss << "            }\n";
    }

    oss << "    } // switch (productionID)\n\n";
    oss << "    return 0;\n";
    oss << "}\n";
}

void ParserTable::generateSymbolsCode(std::ostream& oss, const std::string& symbolVar,
                                      const std::string& className) const {
    oss << "std::string ";
    if (className.size() > 0)
        oss << className << "::";
    oss << "symbolID2String(const int " << symbolVar << ") const {\n";
    oss << "    switch (" << symbolVar << ") {\n";

    for (size_t i = 0; i < symbols_.size(); ++i) {
        if (symbols_[i] == 0)
            continue;

        oss << "        case " << i << ":";
        oss << " return \"" << symbols_[i]->toString() << "\";\n";
    }

    //oss << "        default:\n";
    //oss << "            break;\n";
    oss << "    }  // switch (symbolID)\n";
    oss << "    return \"\";\n";
    oss << "}\n";
}

void ParserTable::htmlOutput(const std::string& fileName) const {
    std::ofstream ofs(fileName.c_str());
    if (ofs.fail())
        return;

    ofs << "<html>\n<head></head>\n<body>\n";
    ofs << "<h1>Action-/Goto-Table</h1>\n";
    ofs << "<table border=\"1\">\n";
    ofs << "  <tr>\n";
    ofs << "    <th>&nbsp;</th>\n";
    for (size_t j = 0; j < symbols_.size(); ++j) {
        if (symbols_[j] == 0)
            continue;

        if ((j % 5) == 0)
            ofs << "    ";

        ofs << "<th>" << symbols_[j]->toString() << "</th>";

        if ((j % 5) == 4)
            ofs << "\n";
    }
    ofs << "  </tr>\n";

    for (int i = 0; i < numStates_; ++i) {
        ofs << "  <tr>\n";
        ofs << "    <td>" << i << "</td>\n";
        const ParserTableRow& row = *(rows_[i]);
        for (size_t j = 0; j < symbols_.size(); ++j) {
            if (symbols_[j] == 0)
                continue;

            if ((j % 5) == 0)
                ofs << "    ";

            ofs << "<td>";
            ParserAction* a = row.findAction(symbols_[j]);
            if (a != 0)
                ofs << a->toString();
            else
                ofs << "&nbsp;";

            ofs << "</td>";

            if ((j % 5) == 4)
                ofs << "\n";

        }   // for (j
        ofs << "  </tr>\n";
    }   // for (i
    ofs << "</table>\n";

    ofs << "<h1>Productions</h1>\n";
    ofs << "<table border=\"1\">\n";
    ofs << "  <tr>\n";
    ofs << "    <th>PID</th>\n";
    ofs << "    <th>Production</th>\n";
    ofs << "  </tr>\n";

    for (size_t i = 0; i < productions_.size(); ++i) {
        ofs << "  <tr>\n";
        ofs << "    <td>" << productions_[i].getProductionID() << "</td>\n";
        ofs << "    <td>" << productions_[i].toString() << "</td>\n";
        ofs << "  </tr>\n";
    }

    ofs << "</table>\n</body>\n</html>\n";
}

// private methods
//

std::string ParserTable::addTabs(const std::string& input, const int numTabs) const {
    if ((numTabs <= 0) || (input.empty() == true))
        return input;

    std::string tabString(numTabs * 4, ' ');
    std::ostringstream oss;
    for (size_t lastPos = 0, pos = 0; pos != std::string::npos; ) {
        pos = input.find('\n', lastPos);
        if (pos != std::string::npos) {
            oss << tabString << input.substr(lastPos, (pos - lastPos) + 1);
            lastPos = pos + 1;
        } else
            oss << tabString << input.substr(lastPos);
    }
    return oss.str();
}

}   // namespace glslparser

}   // namespace voreen
