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

#include "voreen/core/utils/GLSLparser/parser.h"

namespace voreen {

namespace glslparser {

bool operator==(const Parser::ParserSymbol& lhs, const Parser::ParserSymbol& rhs) {
    return (lhs.getSymbolID() == rhs.getSymbolID());
}

bool operator!=(const Parser::ParserSymbol& lhs, const Parser::ParserSymbol& rhs) {
    return (lhs.getSymbolID() != rhs.getSymbolID());
}

// ============================================================================

Parser::Parser(Lexer* const lexer, const bool deleteLexer)
    : debugging_(false),
    lexer_(lexer),
    deleteLexer_(deleteLexer),
    log_(std::ios_base::out | std::ios_base::binary)
{
}

Parser::~Parser() {
    if (deleteLexer_ == true)
        delete lexer_;

    cleanup();
}


ParseTreeNode* Parser::parse() {
    log_.str("");
    states_.push(0);

    Token* token = nextToken();
    bool error = false, accepted = false, fatal = false;

    while ((! accepted) && (! fatal)) {
        if (debugging_) {
            log_ << "current token: " << ((token != 0) ? token->toString() : "NULL") << "\n";
            logState();
        }

        int state = states_.top();
        const int symbolID = ((token != 0) ? token->getTokenID() : -1);
        ParserAction* a = action(state, symbolID);

        if (a == 0) {
            log_ << "  parser error: no action was found for topmost state '" << state;
            log_ << "' and symbol '" << symbolID2String(symbolID) << "'!\n";
            delete a;
            continue;
        }

        switch (a->getType()) {
            case ParserAction::ACTION_TRANSITION:
                if (ParserActionTransition* t = dynamic_cast<ParserActionTransition*>(a)) {
                    states_.push(t->getNextState());

                    if (token != 0)
                        symbols_.push(new ParserSymbol(token));
                    token = nextToken();
                }
                break;

            case ParserAction::ACTION_REDUCE:
                if (! reduce(dynamic_cast<ParserActionReduce*>(a)))
                    fatal = true;
                break;

            case ParserAction::ACTION_ACCEPT:
                if (token != 0)
                    symbols_.push(new ParserSymbol(token));
                token = 0;
                accepted = true;
                break;

            case ParserAction::ACTION_ERROR:
                logError(dynamic_cast<ParserActionError*>(a), token);
                //logState();
                fatal = true;
                // no break here
            default:
                delete token;
                token = nextToken();
                error = true;
                break;
        }   // switch

        delete a;
    }   // while (! accepted

    ParseTreeNode* root = 0;
    if ((accepted == true) && (error == false)) {
        log_ << "  accepted." << std::endl;

        if (nodes_.empty())
            log_ << " !!! Error: parse tree has no  root !!!" << std::endl;

        while (! nodes_.empty()) {
            if (root != 0)
                nodes_.top()->addChild(root);
            root = nodes_.top();
            nodes_.pop();
        }
    }

    cleanup();
    return root;
}

ParseTreeNode* Parser::popNode() {
    ParseTreeNode* top = 0;
    if (! nodes_.empty()) {
        top = nodes_.top();
        nodes_.pop();
    }
    return top;
}

void Parser::pushNode(ParseTreeNode* const node) {
    if (node != 0)
        nodes_.push(node);
}

Token* Parser::nextToken() {
    if (lexer_ != 0)
        return lexer_->scan();

    return 0;
}

// private methods
//

void Parser::cleanup() {
    for ( ; (! nodes_.empty()); nodes_.pop())
        delete nodes_.top();

    for ( ; (! symbols_.empty()); symbols_.pop())
        delete symbols_.top();
}

void Parser::logError(ParserActionError* const error, Token* const token) {
    if (error == 0)
        return;

    std::ostringstream oss;
    if (token != 0) {
        oss << "unexpected symbol '" << symbolID2String(token->getTokenID()) << "'";
        if (GenericToken<std::string>* const gen =
            dynamic_cast<GenericToken<std::string>* const>(token))
        {
            oss << "(= '" << gen->getValue() << "')";
        }
    } // else oss << "unexpected end of input";

    if (! symbols_.empty())
        oss << " following symbol '" << symbolID2String(symbols_.top()->getSymbolID()) << "'";

    error->setErrorMessage(oss.str());
    if (token != 0) {
        error->setLineNumber(token->getLineNumber());
        error->setCharPosition(token->getCharNumber());
    }

    if ((error->getError()).length() > 0)
        log_ << error->getError() << "\n";
}

void Parser::logState() {
    log_ << "states: [";
    std::stack<int> states = states_;
    for (; (! states.empty()); states.pop()) {
        log_ << states.top() << " ";
    }
    log_ << "]\n";

    log_ << "symbols: [";
    std::stack<ParserSymbol*> symbols = symbols_;
    for (; (! symbols.empty()); symbols.pop()) {
        log_ << symbolID2String(symbols.top()->getSymbolID()) << " ";
    }
    log_ << "]\n";

    log_ << "nodes: [";
    std::stack<ParseTreeNode*> nodes = nodes_;
    for (; (! nodes.empty()); nodes.pop()) {
        log_ << symbolID2String(nodes.top()->getSymbolID()) << " ";
    }
    log_ << "]\n\n";
}

bool Parser::reduce(ParserActionReduce* const r) {
    if (r == 0)
        return false;

    const int productionID = r->getProductionID();
    ProductionStub* prod = findProduction(productionID);
    if (prod == 0) {
        log_ << "fatal error: there is no production with ID = " << productionID << "!\n";
        return false;
    }

    bool result = true;
    const std::list<int>& body = prod->getBodyIDs();
    std::vector<ParserSymbol*> removedSymbols;
    for (std::list<int>::const_reverse_iterator it = body.rbegin();
        it != body.rend(); ++it)
    {
        if (! states_.empty())
            states_.pop();
        else {
            log_ << "fatal error: stack is empty!\n";
            result = false;
        }

        if (! symbols_.empty())
        {
            if ((*(symbols_.top()) != *it) && (! isProxySymbol(*it, symbols_.top()->getSymbolID())) )
            {
                log_ << "fatal error: top symbol ('" << symbolID2String(symbols_.top()->getSymbolID());
                log_ << "') does not match the symbol in the production's body ('";
                log_ << symbolID2String(*it) << "')\n";
                log_ << "ProdcutionID = " << productionID << "\n";
                result = false;
            }

            // Attention: this reverses the order of the symbols from the stack!
            //
            removedSymbols.push_back(symbols_.top());
            symbols_.pop();
        }
    }

    symbols_.push(new ParserSymbol(prod->getHeadID()));
    const int nextState = gotoState(states_.top(), prod->getHeadID());
    if (nextState >= 0)
        states_.push(nextState);
    else {
        log_ << "fatal: no following state found for state ";
        log_ << states_.top() << " and symbol " << symbolID2String(prod->getHeadID()) << "\n";

        result = false;
    }
    delete prod;

    if (result == true)
        expandParseTree(productionID, removedSymbols);

    for (size_t i = 0; i < removedSymbols.size(); ++i)
        delete removedSymbols[i];

    return result;
}

}   // namespace glslparser

}   // namespace voreen
