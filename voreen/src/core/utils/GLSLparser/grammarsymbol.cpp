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

#include "voreen/core/utils/GLSLparser/grammarsymbol.h"

#include "voreen/core/utils/GLSLparser/generator/grammar.h"
#include <sstream>

namespace voreen {

namespace glslparser {

GrammarSymbol::~GrammarSymbol() {
}

bool GrammarSymbol::operator==(const GrammarSymbol& rhs) const {
    return ((name_ == rhs.name_) && (productions_.size() == rhs.productions_.size()));
}

bool GrammarSymbol::operator!=(const GrammarSymbol& rhs) const {
    return ((name_ != rhs.name_) || (productions_.size() != rhs.productions_.size()));
}

const std::set<GrammarSymbol*>& GrammarSymbol::getFirstSet(const Grammar& myGrammar) {
    if (! firstSetValid_) {
        firstSetValid_ = true;

        if (isEmptySymbol() == true)
            firstSet_.insert(myGrammar.getEmptySymbol());
        else if (isTerminal() == true)
            firstSet_.insert(this);
        else {
            for (size_t i = 0; i < productions_.size(); ++i) {
                const std::vector<GrammarSymbol*>& body = productions_[i].getBody();

                bool allEmpty = true;
                for (size_t k = 0; k < body.size(); ++k) {
                    std::set<GrammarSymbol*> firstK = body[k]->getFirstSet(myGrammar);
                    std::set<GrammarSymbol*>::iterator it = firstK.find(myGrammar.getEmptySymbol());
                    if (it == firstK.end())
                        allEmpty = false;
                    else
                        firstK.erase(it);

                    if ((allEmpty == false) && (firstK.empty() == false)) {
                        firstSet_.insert(firstK.begin(), firstK.end());
                        break;
                    }
                }   // for (k

                if (allEmpty == true)
                    firstSet_.insert(myGrammar.getEmptySymbol());
            }   // for (i
        }   // else
    }   // if (! firstSetValid_

    return firstSet_;
}

std::string GrammarSymbol::toString() const {
    std::ostringstream oss;

    if (isTerminal() == true)
        oss << name_;
        //oss << "'" << name_ << "'";
    else
        oss << "[" << name_ << "]";

    return oss.str();
}

// private methods
//

GrammarSymbol::GrammarSymbol(const int symbolID, const std::string& name)//, const SymbolFlag flag)
    : symbolID_(symbolID),
    name_(name),
    productions_(),
    firstSetValid_(false),
    flag_(SYMBOL_COMMON)//flag_(flag)
{
}

bool GrammarSymbol::addProduction(const std::vector<GrammarSymbol*>& body,
                                  int& productionIDCounter_)
{
    Production newProduction(productionIDCounter_++, this, body);

    bool exists = false;
    for (size_t i = 0; ((exists == false) && (i < productions_.size())); ++i)
        exists = (productions_[i] == newProduction);

    if (exists == false) {
        productions_.push_back(newProduction);
        return true;
    } else
        --productionIDCounter_;

    return false;
}

}   // namespace glslparser

}   // namespace voreen
