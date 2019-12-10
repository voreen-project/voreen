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

#ifndef VRN_GRAMMAR_H
#define VRN_GRAMMAR_H

#include "voreen/core/utils/GLSLparser/grammarsymbol.h"
#include "voreen/core/utils/GLSLparser/generator/lritem.h"
#include "voreen/core/utils/GLSLparser/terminalsmap.h"

#include <list>

namespace voreen {

namespace glslparser {

class ParserTable;

class Grammar {
public:
    Grammar();
    virtual ~Grammar();

    /**
     * Implemented though pure virtual to initialize terminals! Should ALWAYS
     * be called by implementations in derived classes.
     */
    virtual void setup(const TerminalsMap& terminals) = 0;

    void addTerminal(const std::string& name, const int symbolID,
        const GrammarSymbol::SymbolFlag flag = GrammarSymbol::SYMBOL_COMMON);
    void addProduction(const std::string& head, const std::string& body);
    void addProduction(const std::string& head, const std::vector<std::string>& body);

    GrammarSymbol* getEmptySymbol() const { return emptySymbol_; }
    GrammarSymbol* getEndSymbol() const { return endSymbol_; }
    GrammarSymbol* getStartSymbol() const { return startSymbol_; }

    const GrammarSymbol* findSymbol(const std::string& name) const;
    GrammarSymbol* findSymbol(const std::string& name);
    std::set<GrammarSymbol*> getNonTerminals() const;
    std::set<GrammarSymbol*> getTerminals() const;
    std::vector<Production> getProductions() const;
    std::vector<GrammarSymbol*> getSymbols() const;

    /**
     * Closure for the given set of items.
     * The caller MUST delete the returned pointer using C++ operator delete.
     */
    LR1ItemSet* closure(LR1ItemSet* const items) const;

    /**
     * Set of symbols which may preceed the given production.
     */
    std::set<GrammarSymbol*> first(const Production& production) const;

    /**
     * Set of Items representing the next state to which the parse can go into
     * when reading the given symbol.
     * The caller MUST delete the returned pointer using C++ operator delete.
     */
    LR1ItemSet* gotoNext(const LR1ItemSet& items, GrammarSymbol* const symbol) const;


    /**
     * Creates the set of canonical LR1 item sets. Theses sets of LR1 items may only
     * differ in their look-ahead symbols and share the same kernel.
     * The sets have to be unioned to create a LALR(1) parser.
     *
     * They could be used to generate a simple shift-reduce parser, but this is not
     * recommmended for larger grammars. The memory consumption of this method
     * is also very extensive.
     *
     * Rather use <code>createdUinionedItemSets()</code> instead.
     *
     * NOTE: after all, the elements of the returned list have to be deleted by the
     * caller using C++ delete operator! Otherwise you will create a HUGE memory
     * leak.
     */
    std::list<LR1ItemSet*> createCanonicalItemSets(LR1ItemSet* const init) const;

    /**
     * Creates a set of LR1 item sets in which two sets each have different kernels.
     * These set of item sets can be used directly to create a LALR(1) parser and
     * it consumes much less memory and it is also faster due to some optimizations.
     *
     * NOTE: after all, the elements of the returned list have to be deleted by the
     * caller using C++ delete operator!
     */
    std::list<LR1ItemSet*> createUnionedItemSets(LR1ItemSet* const init) const;

    ParserTable* createParserTable(const bool unionAtOnce) const;

private:
    bool addTransitionAction(ParserTable* const table, const int state,
        const int gotoState, GrammarSymbol* const symbol) const;
    LR0ItemSet getKernel(const LR1ItemSet& itemSet) const;
    void unionItemSets(std::list<LR1ItemSet*>& canonical) const;

private:
    GrammarSymbol* emptySymbol_;
    GrammarSymbol* endSymbol_;
    GrammarSymbol* startSymbol_;

    typedef std::map<std::string, GrammarSymbol*> Symbols;
    Symbols symbols_;
    int productionIDCounter_;
    int nextSymbolID_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
