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

#ifndef VRN_GRAMMARSYMBOL_H
#define VRN_GRAMMARSYMBOL_H

#include "voreen/core/utils/GLSLparser/production.h"

#include <set>
#include <string>

namespace voreen {

namespace glslparser {

class Grammar;

class GrammarSymbol {
friend class Grammar;

public:
    enum SymbolFlag {
        SYMBOL_COMMON = 0,
        SYMBOL_START = 1,
        SYMBOL_EMPTY = 2,
        SYMBOL_END = 3
    };

public:
    ~GrammarSymbol();
    bool operator==(const GrammarSymbol& rhs) const;
    bool operator!=(const GrammarSymbol& rhs) const;

    const std::set<GrammarSymbol*>& getFirstSet(const Grammar& myGrammar);
    SymbolFlag getFlag() const { return flag_; }
    const std::string& getName() const { return name_; }
    const std::vector<Production>& getProductions() const { return productions_; }
    int getSymbolID() const { return symbolID_; }

    bool isEmptySymbol() const { return (flag_ == SYMBOL_EMPTY); }
    bool isEndSymbol() const { return (flag_ == SYMBOL_END); }
    bool isStartSymbol() const { return (flag_ == SYMBOL_START); }
    bool isTerminal() const { return productions_.empty(); }

    std::string toString() const;

private:
    GrammarSymbol(const int symbolID, const std::string& name);//, const SymbolFlag flag = SYMBOL_COMMON);
    GrammarSymbol(const GrammarSymbol&);
    GrammarSymbol& operator=(const GrammarSymbol&);

    bool addProduction(const std::vector<GrammarSymbol*>& body, int& productionIDCounter);
    void setFlag(const SymbolFlag flag) { flag_ = flag; }

private:
    int symbolID_;
    std::string name_;
    std::vector<Production> productions_;
    std::set<GrammarSymbol*> firstSet_;
    bool firstSetValid_;
    SymbolFlag flag_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
