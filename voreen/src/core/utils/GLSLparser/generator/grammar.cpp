/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/GLSLparser/generator/grammar.h"
#include "voreen/core/utils/GLSLparser/generator/parsertable.h"

#include <iostream>
#include <sstream>

namespace voreen {

namespace glslparser {

Grammar::Grammar()
    : emptySymbol_(new GrammarSymbol(0x80000000, "$EMPTY$")),
    endSymbol_(new GrammarSymbol(-1, "$END$")),
    startSymbol_(new GrammarSymbol(0, "$START$")),
    productionIDCounter_(0),
    nextSymbolID_(1)
{
    // Setup some special symbols required from operations on this
    // grammar
    //
    emptySymbol_->setFlag(GrammarSymbol::SYMBOL_EMPTY);
    endSymbol_->setFlag(GrammarSymbol::SYMBOL_END);
    startSymbol_->setFlag(GrammarSymbol::SYMBOL_START);
}

Grammar::~Grammar() {
    for (Symbols::iterator it = symbols_.begin(); it != symbols_.end(); ++it) {
        // Do not delete the special symbols within this loop, rather delete
        // them afterwards for they might still be the default symbols
        //
        if (it->second->getFlag() == GrammarSymbol::SYMBOL_COMMON)
            delete it->second;
    }
    symbols_.clear();

    // Delete special symbols
    //
    delete emptySymbol_;
    delete endSymbol_;
    delete startSymbol_;
}

void Grammar::addTerminal(const std::string& name, const int symbolID,
                          const GrammarSymbol::SymbolFlag flag)
{
    GrammarSymbol* s = new GrammarSymbol(symbolID, name);

    std::pair<Symbols::iterator, bool> res = symbols_.insert(std::make_pair(name, s));
    if (! res.second)
        delete s;
    else {
        if (symbolID >= nextSymbolID_)
            nextSymbolID_ = symbolID + 1;

        // Identify special symbols among the terminals and use them instead
        // of the default values for this grammar.
        //
        switch (flag) {
            case GrammarSymbol::SYMBOL_EMPTY:
                s->setFlag(flag);
                emptySymbol_ = s;
                break;

            case GrammarSymbol::SYMBOL_END:
                s->setFlag(flag);
                endSymbol_ = s;
                break;

            case GrammarSymbol::SYMBOL_START:
                s->setFlag(flag);
                startSymbol_ = s;
                break;

            default:
                break;
        }   // switch
    }
}

void Grammar::addProduction(const std::string& head, const std::string& body) {
    std::vector<std::string> bodyVec;

    // Split the string at positions containing space characters and pass
    // it on to the overload method.
    //
    if (! body.empty()) {
        size_t pos = 0, prevPos = 0;
        do {
            prevPos = pos;
            pos = body.find(' ', pos);

            std::string part;
            if (pos == std::string::npos)
                part = body.substr(prevPos);
            else
                part = body.substr(prevPos, (pos++ - prevPos));

            if (! part.empty())
                bodyVec.push_back(part);

        } while (pos != std::string::npos);
    }

    addProduction(head, bodyVec);
}

void Grammar::addProduction(const std::string& head, const std::vector<std::string>& body) {
    // Try to find already existing parts of the body of the new production
    // among the productions or add them otherwise.
    //
    std::vector<GrammarSymbol*> productionBody;
    for (size_t i = 0; i < body.size(); ++i) {
        Symbols::const_iterator it = symbols_.find(body[i]);

        if (it != symbols_.end())
            productionBody.push_back(it->second);
        else {
            GrammarSymbol* p = new GrammarSymbol(nextSymbolID_, body[i]);
            productionBody.push_back(p);
            symbols_.insert(std::make_pair(body[i], p));
            ++nextSymbolID_;
        }
    }

    // Try to insert the new production into the symbols map.
    //
    GrammarSymbol* symbol = new GrammarSymbol(nextSymbolID_, head);
    std::pair<Symbols::iterator, bool> res = symbols_.insert(std::make_pair(head, symbol));

    // If a production for that name already exists and the new production is not identical to
    // the existing one, insert the body as a new alternative
    //
    if (! productionBody.empty()) {
        if (res.second == false) {
            GrammarSymbol* const existing = res.first->second;
            existing->addProduction(productionBody, productionIDCounter_);
            delete symbol;
        } else {
            symbol->addProduction(productionBody, productionIDCounter_);
            ++nextSymbolID_;
        }
    }
}

const GrammarSymbol* Grammar::findSymbol(const std::string& name) const {
    Symbols::const_iterator it = symbols_.find(name);
    if (it == symbols_.end())
        return 0;
    return it->second;
}

GrammarSymbol* Grammar::findSymbol(const std::string& name) {
    Symbols::iterator it = symbols_.find(name);
    if (it == symbols_.end())
        return 0;
    return it->second;
}

std::set<GrammarSymbol*> Grammar::getTerminals() const {
    std::set<GrammarSymbol*> terminals;
    for (Symbols::const_iterator it = symbols_.begin(); it != symbols_.end(); ++it) {
        if (it->second->isTerminal())
            terminals.insert(it->second);
    }
    return terminals;
}

std::set<GrammarSymbol*> Grammar::getNonTerminals() const {
    std::set<GrammarSymbol*> nonTerminals;
    for (Symbols::const_iterator it = symbols_.begin(); it != symbols_.end(); ++it) {
        if (! it->second->isTerminal())
            nonTerminals.insert(it->second);
    }
    return nonTerminals;
}

std::vector<Production> Grammar::getProductions() const {
    std::vector<Production> productions;

    for (Symbols::const_iterator itSym = symbols_.begin();
        itSym != symbols_.end(); ++itSym)
    {
        if (itSym->second->isTerminal())
            continue;

        const std::vector<Production>& p = itSym->second->getProductions();
        productions.insert(productions.end(), p.begin(), p.end());
    }
    return productions;
}

std::vector<GrammarSymbol*> Grammar::getSymbols() const {
    // Create a vector large enough to index the symbols by their IDs
    //
    std::vector<GrammarSymbol*> symVec(nextSymbolID_);
    for (Symbols::const_iterator it = symbols_.begin(); it != symbols_.end(); ++it) {
        const int symID = it->second->getSymbolID();
        if (symID >= static_cast<int>(symVec.size()))
            continue;

        symVec[symID] = it->second;
    }

    symVec.push_back(endSymbol_);
    return symVec;
}

void Grammar::setup(const TerminalsMap& terminals) {
    const std::map<std::string, int>& terminalSymbols = terminals.getLexemeMap();
    for (std::map<std::string, int>::const_iterator it = terminalSymbols.begin();
        it != terminalSymbols.end(); ++it)
    {
        addTerminal(it->first, it->second);
    }
}

LR1ItemSet* Grammar::closure(LR1ItemSet* const items) const {
    //LR1ItemSet* closure = items;
//std::cout << "closure (1) = " << closure.toString() << "\n";

    for (size_t i = 0; i < items->getNumItems(); ++i) {
        // Get the symbol right beneath the dot's position
        //
        GrammarSymbol* const s = (*items)[i].getRightSymbol();
        if (s == 0)
            continue;

        // Terminals do not have any productions, so skip them if encountered
        // to avoid the calculations of the postfix and the coresponding first-set
        //
        Symbols::const_iterator prodIt = symbols_.find(s->getName());
        if ((prodIt == symbols_.end()) || (prodIt->second->isTerminal()))
            continue;

        // Retrieve all remaining symbols right from the position of the dot and
        // calculate its first-set.
        //
        Production postfix = (*items)[i].getRightString();
//std::cout << "    postfix (1) = " << postfix.toString() << "\n";
        postfix.removeSymbol(0);    // remove the current symbol itself
        postfix.appendSymbol((*items)[i].getLookAhead());    // append the lookahead symbol
        std::set<GrammarSymbol*> firstPostfix = first(postfix);

/*std::cout << "    postfix (2) = " << postfix.toString() << "\n";
std::cout << "first Set for production " << postfix.toString() << ":\n";
std::cout << "  {";
for (std::set<GrammarSymbol*>::const_iterator it = firstPostfix.begin(); it != firstPostfix.end(); ++it)
    std::cout << (*it)->toString() << ", ";
std::cout << "}\n\n";*/

        const std::vector<Production>& productions = prodIt->second->getProductions();
        for (size_t p = 0; p < productions.size(); ++p)
        {
            for (std::set<GrammarSymbol*>::const_iterator itFirst = firstPostfix.begin();
                itFirst != firstPostfix.end(); ++itFirst)
            {
                if ((*itFirst)->isTerminal())
                    items->addItem(LR1Item(s->getName(), productions[p], 0, *itFirst));
            }   // for (itFirst
        }   // for (p
    }   // for (i

    return items;
}

std::set<GrammarSymbol*> Grammar::first(const Production& production) const {
    std::set<GrammarSymbol*> firstSet;
    const std::vector<GrammarSymbol*>& body = production.getBody();
    bool allContainEmptySymbol = true;

    //for (size_t i = 0; ((foundEmptySymbol == false) && (i < body.size())); ++i) {
    for (size_t i = 0; i < body.size(); ++i) {
        bool foundEmptySymbol = (i > 0) ? false : true;
        std::set<GrammarSymbol*> firstI = body[i]->getFirstSet(*this);
        std::set<GrammarSymbol*>::iterator itEmpty = firstI.find(emptySymbol_);

        if (itEmpty != firstI.end()) {
            foundEmptySymbol = true;
            firstI.erase(itEmpty);
        } else
            allContainEmptySymbol = false;

        // Add all non-empty symbols, if the empty symbol was found in firstI
        //
        if ((! firstI.empty()) && (foundEmptySymbol))
            firstSet.insert(firstI.begin(), firstI.end());
    }   // for (i

    if (allContainEmptySymbol)
        firstSet.insert(emptySymbol_);

    return firstSet;
}

LR1ItemSet* Grammar::gotoNext(const LR1ItemSet& items, GrammarSymbol* const symbol) const
{
    //std::vector<LR1Item> gotoSet;
    LR1ItemSet* gotoSet = new LR1ItemSet();

    for (size_t i = 0; i < items.getNumItems(); ++i) {
        const LR1Item& item = items[i];
        GrammarSymbol* const right = item.getRightSymbol();
        if ((right != 0) && (*right == *symbol)) {
            gotoSet->addItem( LR1Item(item.getName(), item.getProduction(),
                item.getPosition() + 1, item.getLookAhead()) );
        }
    }

    return closure(gotoSet);
}

std::list<LR1ItemSet*> Grammar::createCanonicalItemSets(LR1ItemSet* const init) const
{
    std::list<LR1ItemSet*> collection;
    collection.push_back(closure(init));

    int numItemSets = 0;
    for (std::list<LR1ItemSet*>::iterator it = collection.begin();
        it != collection.end(); ++it)
    {
std::cout << "creating ItemSet #" << numItemSets++ << "...\n";
//std::cout << "  I" << (counter - 1) << " = " << it->toString() << "\n";
        if (! (*it)->hasNextSet())
            continue;

        for (Symbols::const_iterator itSym = symbols_.begin();
            itSym != symbols_.end(); ++itSym)
        {
            LR1ItemSet* gotoSet = gotoNext(*(*it), itSym->second);
//std::cout << "  Symbol X = " << itSym->second->getName().c_str() << "\n";
//std::cout << "    => gotoSet = " << gotoSet.toString() << "\n";

            if ((! gotoSet->isEmpty()) && (! gotoSet->contained(collection)))
                collection.push_back(gotoSet);
            else
                delete gotoSet;
        }   // for (itSym
    } // for (it

std::cout << "\n";
    return collection;
}

std::list<LR1ItemSet*> Grammar::createUnionedItemSets(LR1ItemSet* const init) const
{
    std::list<LR1ItemSet*> collection;
    collection.push_back(closure(init));

    std::list<LR1ItemSet> candidates;
    candidates.push_back(*(collection.back()));

    int numItemSets = 0;
    for (std::list<LR1ItemSet>::iterator it = candidates.begin(); it != candidates.end(); ++it)
    {
std::cout << "investigating ItemSet #" << numItemSets++ << "...\n";

        if (! it->hasNextSet())
            continue;

        for (Symbols::const_iterator itSym = symbols_.begin(); itSym != symbols_.end(); ++itSym)
        {
            LR1ItemSet* gotoSet = gotoNext(*it, itSym->second);
            if (! gotoSet->isEmpty()) {
                int numMerged = -1;
                std::list<LR1ItemSet*>::iterator existing = collection.begin();
                for (; ((numMerged < 0) && (existing != collection.end())); ++existing) {
                    if ( *((*existing)->getKernel()) == *(gotoSet->getKernel()))
                        numMerged = (*existing)->merge(*gotoSet);
                }

                if (numMerged < 0) {
                    collection.push_back(gotoSet);
                    candidates.push_back(*gotoSet); // insert copy to candidates list
                } else if (numMerged > 0) {
                    candidates.push_back(*gotoSet); // insert copy to candidates list
                    delete gotoSet;
                } else
                    delete gotoSet;
            } else
                delete gotoSet;
        }   // for (itSym
    } // for (it

std::cout << "\n";
    return collection;
}

ParserTable* Grammar::createParserTable(const bool unionAtOnce) const {
    if (startSymbol_ == 0) {
        std::cout << "Grammar::createParserTable(): start symbol is NULL!\n";
        return 0;
    }

    if (startSymbol_->isTerminal()) {
        std::cout << "Grammar::createParserTable(): start symbol is lacks any productions!\n";
        return 0;
    }

    LR1ItemSet* init = new LR1ItemSet();
    init->addItem(LR1Item(startSymbol_->getName(), startSymbol_->getProductions()[0], 0, endSymbol_));

    typedef std::list<LR1ItemSet*> SetList;
    SetList collection;
    if (! unionAtOnce) {
        collection = createCanonicalItemSets(init);

        // Retrieve the union of all item-sets with the same kernel and the number of
        // distinct sets.
        //
        std::cout << "number of ItemSets before union: " << collection.size() << "\n";
        unionItemSets(collection);
    } else
        collection = createUnionedItemSets(init);

    std::cout << "number of ItemSets after union: " << collection.size() << "\n";

    ParserTable* table = new ParserTable(static_cast<int>(collection.size()), getSymbols(),
        getProductions());

    // map holding all previously calculated goto-sets for each on the right-hand side
    // of the productions and their positions within the collection.
    //
    typedef std::map<int, int> GotoSets;

    bool noConflict = true;
    SetList::const_iterator it = collection.begin();
    size_t totalNumItemSets = 0;
    int totalNumSaved = 0, totalNumGotos = 0;
    for (int stateNo = 0; it != collection.end(); ++it, ++stateNo) {

        LR1ItemSet* const itemSet = *it;
        GotoSets gotos;
        int saved = 0;
        int numGotos = 0;

        for (size_t i = 0; i < itemSet->getNumItems(); ++i) {
            GrammarSymbol* const right = (*itemSet)[i].getRightSymbol();

            if (right != 0) {
                int gotoState = -1;
                GotoSets::const_iterator existingGoto = gotos.find(right->getSymbolID());

                if (existingGoto == gotos.end()) {
                    LR1ItemSet* gotoSet = gotoNext(*itemSet, right);
                    if (! gotoSet->isEmpty()) {
                        gotoState = gotoSet->findKernelsPositionIn(collection);
                        if (gotoState != -1) {
                            noConflict = addTransitionAction(table, stateNo, gotoState, right);

                            // Add this new goto set to the map so that it can be reused.
                            //
                            //std::pair<LR1ItemSet*, int> second(gotoSet, gotoState);
                            gotos.insert( std::make_pair(right->getSymbolID(), gotoState) );
                        }
                    }
                    delete gotoSet;
                    ++numGotos;
                } else {
                    gotoState = existingGoto->second;
                    noConflict = addTransitionAction(table, stateNo, gotoState, right);
                    ++saved;
                }

                if (! noConflict) {
                    const Production& p = (*itemSet)[i].getProduction();
                    GrammarSymbol* const lookahead = (*itemSet)[i].getLookAhead();

                    std::cout << "  failed to add transistion-action for state " << stateNo << "!\n";
                    std::cout << "    production = " << p.toString() << "\n";
                    std::cout << "    lookahead = " << lookahead->toString() << "\n";
                    ParserAction* existing = table->findAction(stateNo, lookahead);
                    if (existing != 0) {
                        std::cout << "    exsiting action = ";
                        std::cout << existing->toString() << "\n";
                    }
                    std::cout << "    new state = " << gotoState << "\n";
                }
            } else {    // right symbol is not present
                const GrammarSymbol* const head = (*itemSet)[i].getProduction().getHead();
                GrammarSymbol* const lookahead = (*itemSet)[i].getLookAhead();

                if ((head->isStartSymbol()) &&  (*lookahead == *endSymbol_))
                    noConflict = table->addAction(stateNo, endSymbol_, new ParserActionAccept());
                else {
                    const Production& p = (*itemSet)[i].getProduction();
                    ParserAction* action = new ParserActionReduce(p.getProductionID());

                    noConflict = table->addAction(stateNo, lookahead, action);
                    if (! noConflict) {
                        std::cout << "  failed to add reduce-action for state " << stateNo << "!\n";
                        std::cout << "    production = " << p.toString() << "\n";
                        std::cout << "    lookahead = " << lookahead->toString() << "\n";
                        std::cout << "    exsiting action = ";
                        std::cout << table->findAction(stateNo, lookahead)->toString() << "\n";
                        std::cout << "    new action = " << action->toString() << "\n";
                        delete action;
                    }
                }   // else (head != start symbol
            }   // else (right symbol present
        }   // for (i

        totalNumItemSets += itemSet->getNumItems();
        totalNumGotos += numGotos;
        totalNumSaved += saved;
    }   // for (it

    for (SetList::iterator c = collection.begin(); c != collection.end(); ++c)
        delete *c;
    collection.clear();

/*
    if (noConflict == false) {
        std::cout << "A conflict among the parser actions occured while creating ";
        std::cout << "the parser table! No table has been created.\n";
        delete table;
        table = 0;
    }
*/
    std::cout << "\ntotal number of ItemSets: " << totalNumItemSets << "\n";
    std::cout << "total number of goto set calculations: " << totalNumGotos << "\n";
    std::cout << "total number of saved goto set calculations: " << totalNumSaved << "\n";
    std::cout << "\nfinished creating parser table.\n";

    return table;
}

// private methods
//

bool Grammar::addTransitionAction(ParserTable* const table, const int state,
                                  const int gotoState, GrammarSymbol* const symbol) const
{
    bool noConflict = true;
    ParserAction* const existing = table->findAction(state, symbol);
    if (existing == 0) {
        ParserAction* action = new ParserActionTransition(gotoState);
        table->addAction(state, symbol, action);
    } else {
        // If there is also an action present, it is only a conflict,
        // if the actions are of different types or if they do not
        // transition to the same next state.
        //
        ParserActionTransition* trans = dynamic_cast<ParserActionTransition*>(existing);
        if ((trans == 0) || (trans->getNextState() != gotoState))
            noConflict = false;
    }
    return noConflict;
}

void Grammar::unionItemSets(std::list<LR1ItemSet*>& canonical) const {
    typedef std::list<LR1ItemSet*> ItemCollection;

    for (ItemCollection::iterator left = canonical.begin(); left != canonical.end(); ++left) {
        ItemCollection::iterator inner = left;

        for (++inner; inner != canonical.end(); ) {
            // If there is another ItemSet with the same kernel, move all Items with
            // different lookahead symbols to the current ItemSet and delete that
            // redundant ItemSet
            //
            if (*((*left)->getKernel()) == *((*inner)->getKernel())) {
                (*left)->merge(*(*inner));
                canonical.erase(inner++);
            } else
                ++inner;
        }   // for (inner
    }   // for (left
}

}   // namespace glslparser

}   // namespace voreen
