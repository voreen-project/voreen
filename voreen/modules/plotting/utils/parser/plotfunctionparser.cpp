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

#include "plotfunctionparser.h"

#include "plotfunctionlexer.h"
#include "plotfunctionnode.h"
#include "plotfunctionterminals.h"

namespace voreen {

namespace glslparser {

PlotFunctionParser::PlotFunctionParser(PlotFunctionLexer *const lexer)
    : Parser(lexer,false)
    , tokens_()
{}

PlotFunctionParser::PlotFunctionParser(std::istream* const is) 
    : Parser(new PlotFunctionLexer(is), true)
    , tokens_()
{}

PlotFunctionParser::PlotFunctionParser(const std::string& fileName) 
    : Parser(new PlotFunctionLexer(fileName,  PlotFunctionTerminals().getKeywords()), true)
    , tokens_()
{}

PlotFunctionParser::PlotFunctionParser(const std::vector<Token*>& tokens)
    : Parser(0)
{
    for (size_t i = 0; i < tokens.size(); ++i) {
        tokens_.push_back(tokens[i]->getCopy());
    }
}

PlotFunctionParser::~PlotFunctionParser() {
    while (! tokens_.empty()) {
        delete tokens_.front();
        tokens_.pop_front();
    }
}

// protected

Token* PlotFunctionParser::nextToken() {
    if (! tokens_.empty()) {
        Token* const t = tokens_.front();
        tokens_.pop_front();
        return t;
    }
    else
        return Parser::nextToken();

    return 0;
}

void PlotFunctionParser::expandParseTree(const int productionID,
                                       const std::vector<Parser::ParserSymbol*>& reductionBody)
{
    switch (productionID) {
        case 0:  // [$START$] ::= [expression]
            break;
        case 1:  // [expression] ::= [term]
            break;
        case 2:  // [expression] ::= [term-list]
            break;
        case 3:  // [expression] ::= [term-list] ; [term]
            if (reductionBody.size() == 3) {
                IsolatorToken* const isolator =
                    dynamic_cast<IsolatorToken* const>(reductionBody[1]->getToken());
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((isolator != 0) && (preterm != 0) && (postterm != 0)) {
                    pushNode(new PlotFunctionIsolatorTerm(isolator,preterm,postterm));
                }
            }
            break;
        case 4:  // [expression] ::= EMPTY
            {
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (name != 0)
                    pushNode(new PlotFunctionTagName(name));
            }
            break;
        case 5:  // [term-list] ::= [termwithinterval]
            break;
        case 6:  // [term-list] ::= [term-list] ; [termwithinterval]
            if (reductionBody.size() == 3) {
                IsolatorToken* const isolator =
                    dynamic_cast<IsolatorToken* const>(reductionBody[1]->getToken());
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((isolator != 0) && (preterm != 0) && (postterm != 0)) {
                    pushNode(new PlotFunctionIsolatorTerm(isolator,preterm,postterm));
                }
            }
            break;
        case 7:  // [termwithinterval] ::= [term] : [interval]
            if (reductionBody.size() == 3) {
                IsolatorToken* const isolator =
                    dynamic_cast<IsolatorToken* const>(reductionBody[1]->getToken());
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((isolator != 0) && (preterm != 0) && (postterm != 0)) {
                    pushNode(new PlotFunctionIsolatorTerm(isolator,preterm,postterm));
                }
            }
            break;
        case 8:  // [interval] ::= [range]
            break;
        case 9:  // [interval] ::= [interval] | [range]
            if (reductionBody.size() == 3) {
                IsolatorToken* const isolator =
                    dynamic_cast<IsolatorToken* const>(reductionBody[1]->getToken());
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((isolator != 0) && (preterm != 0) && (postterm != 0)) {
                    pushNode(new PlotFunctionIsolatorTerm(isolator,preterm,postterm));
                }
            }
            break;
        case 10:  // [range] ::= - [range]
            if (reductionBody.size() == 2) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(1);
                PlotFunctionTerm* const term = popNode<PlotFunctionTerm>();
                if ((op != 0) && (term != 0))
                    pushNode(new PlotFunctionPreOperation(op,term));
            }
            break;
        case 11:  // [range] ::= [leftrange] , [rightrange]
            if (reductionBody.size() == 3) {
                IsolatorToken* const isolator =
                    dynamic_cast<IsolatorToken* const>(reductionBody[1]->getToken());
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((isolator != 0) && (preterm != 0) && (postterm != 0)) {
                    pushNode(new PlotFunctionIsolatorTerm(isolator,preterm,postterm));
                }
            }
            break;
        case 12:  // [leftrange] ::= ( [interval-term]
            if (reductionBody.size() == 2) {
                BracketToken* const op = dynamic_cast<BracketToken* const>(reductionBody[1]->getToken());
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (postterm != 0))
                    pushNode(new PlotFunctionRangeTerm(op,postterm));
            }
            break;
        case 13:  // [leftrange] ::= [ [interval-term]
            if (reductionBody.size() == 2) {
                BracketToken* const op = dynamic_cast<BracketToken* const>(reductionBody[1]->getToken());
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (postterm != 0))
                    pushNode(new PlotFunctionRangeTerm(op,postterm));
            }
            break;
        case 14:  // [rightrange] ::= [interval-term] )
            if (reductionBody.size() == 2) {
                BracketToken* const op = dynamic_cast<BracketToken* const>(reductionBody[0]->getToken());
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0))
                    pushNode(new PlotFunctionRangeTerm(op,preterm));
            }
            break;
        case 15:  // [rightrange] ::= [interval-term] ]
            if (reductionBody.size() == 2) {
                BracketToken* const op = dynamic_cast<BracketToken* const>(reductionBody[0]->getToken());
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0))
                    pushNode(new PlotFunctionRangeTerm(op,preterm));
            }
            break;
        case 16:  // [term] ::= [term] + [addend]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 17:  // [term] ::= [term] - [addend]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 18:  // [term] ::= + [addend]
            if (reductionBody.size() == 2) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(1);
                PlotFunctionTerm* const term = popNode<PlotFunctionTerm>();
                if ((op != 0) && (term != 0))
                    pushNode(new PlotFunctionPreOperation(op,term));
            }
            break;
        case 19:  // [term] ::= - [addend]
            if (reductionBody.size() == 2) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(1);
                PlotFunctionTerm* const term = popNode<PlotFunctionTerm>();
                if ((op != 0) && (term != 0))
                    pushNode(new PlotFunctionPreOperation(op,term));
            }
            break;
        case 20:  // [term] ::= - - [addend]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = new OperatorToken(PlotFunctionTerminals::ID_PLUS,'+');
                op->setParameter(1);
                PlotFunctionTerm* const term = popNode<PlotFunctionTerm>();
                if ((op != 0) && (term != 0))
                    pushNode(new PlotFunctionPreOperation(op,term,true));
            }
            break;
        case 21:  // [term] ::= [term] + - [addend]
            if (reductionBody.size() == 4) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 22:  // [term] ::= [term] - + [addend]
            if (reductionBody.size() == 4) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[2]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 23:  // [term] ::= [term] - - [addend]
            if (reductionBody.size() == 4) {
                OperatorToken* const op = new OperatorToken(PlotFunctionTerminals::ID_PLUS,'+');
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm,true));
            }
            break;
        case 24:  // [term] ::= [addend]
            break;
        case 25:  // [addend] ::= [factor]
            break;
        case 26:  // [addend] ::= [addend] * [factor]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 27:  // [addend] ::= [addend] / [factor]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 28:  // [factor] ::= [variable]
            break;
        case 29:  // [factor] ::= [variable] ^ [factor]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 30:  // [variable] ::= [atom]
            break;
        case 31:  // [variable] ::= [bracket-term]
            break;
        case 32:  // [bracket-term] ::= ( [term] )
            break;
        case 33:  // [atom] ::= [function] [bracket-term]
            if (reductionBody.size() == 2) {
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionFunctionTerm(preterm,postterm));
            }
            break;
        case 34:  // [atom] ::= [single-value]
            break;
        case 35:  // [single-value] ::= INTCONST
            if (reductionBody.size() == 1) {
                IdentifierToken* const token = dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (token != 0)
                    pushNode(new PlotFunctionTagValues(token));
            }
            break;
        case 36:  // [single-value] ::= FLOATCONST
            if (reductionBody.size() == 1) {
                IdentifierToken* const token = dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (token != 0)
                    pushNode(new PlotFunctionTagValues(token));
            }
            break;
        case 37:  // [single-value] ::= VARIABLE
            if (reductionBody.size() == 1) {
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (name != 0)
                    pushNode(new PlotFunctionTagName(name));
            }
            break;
        case 38:  // [function] ::= FUNCTION
            if (reductionBody.size() == 1) {
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (name != 0)
                    pushNode(new PlotFunctionTagName(name));
            }
            break;
        case 39:  // [interval-term] ::= [interval-term] + [interval-addend]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 40:  // [interval-term] ::= [interval-term] - [interval-addend]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 41:  // [interval-term] ::= + [interval-addend]
            if (reductionBody.size() == 2) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(1);
                PlotFunctionTerm* const term = popNode<PlotFunctionTerm>();
                if ((op != 0) && (term != 0))
                    pushNode(new PlotFunctionPreOperation(op,term));
            }
            break;
        case 42:  // [interval-term] ::= - [interval-addend]
            if (reductionBody.size() == 2) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(1);
                PlotFunctionTerm* const term = popNode<PlotFunctionTerm>();
                if ((op != 0) && (term != 0))
                    pushNode(new PlotFunctionPreOperation(op,term));
            }
            break;
        case 43:  // [interval-term] ::= - - [interval-addend]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = new OperatorToken(PlotFunctionTerminals::ID_PLUS,'+');
                op->setParameter(1);
                PlotFunctionTerm* const term = popNode<PlotFunctionTerm>();
                if ((op != 0) && (term != 0))
                    pushNode(new PlotFunctionPreOperation(op,term,true));
            }
            break;
        case 44:  // [interval-term] ::= [interval-term] + - [interval-addend]
            if (reductionBody.size() == 4) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 45:  // [interval-term] ::= [interval-term] - + [interval-addend]
            if (reductionBody.size() == 4) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[2]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 46:  // [interval-term] ::= [interval-term] - - [interval-addend]
            if (reductionBody.size() == 4) {
                OperatorToken* const op = new OperatorToken(PlotFunctionTerminals::ID_PLUS,'+');
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm,true));
            }
            break;
        case 47:  // [interval-term] ::= [interval-addend]
            break;
        case 48:  // [interval-addend] ::= [interval-factor]
            break;
        case 49:  // [interval-addend] ::= [interval-addend] * [interval-factor]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 50:  // [interval-addend] ::= [interval-addend] / [interval-factor]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 51:  // [interval-factor] ::= [interval-variable]
            break;
        case 52:  // [interval-factor] ::= [interval-variable] ^ [interval-factor]
            if (reductionBody.size() == 3) {
                OperatorToken* const op = dynamic_cast<OperatorToken* const>(reductionBody[1]->getToken());
                op->setParameter(2);
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((op != 0) && (preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionCenterOperation(op,preterm,postterm));
            }
            break;
        case 53:  // [interval-variable] ::= [interval-atom]
            break;
        case 54:  // [interval-variable] ::= [interval-bracket-term]
            break;
        case 55:  // [interval-bracket-term] ::= ( [interval-term] )
            break;
        case 56:  // [interval-atom] ::= [function] [interval-bracket-term]
            if (reductionBody.size() == 2) {
                PlotFunctionTerm* const postterm = popNode<PlotFunctionTerm>();
                PlotFunctionTerm* const preterm = popNode<PlotFunctionTerm>();
                if ((preterm != 0) && (postterm != 0))
                    pushNode(new PlotFunctionFunctionTerm(preterm,postterm));
            }
            break;
        case 57:  // [interval-atom] ::= [interval-single-value]
            break;
        case 58:  // [interval-single-value] ::= INTCONST
            if (reductionBody.size() == 1) {
                IdentifierToken* const token = dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (token != 0)
                    pushNode(new PlotFunctionTagValues(token));
            }
            break;
        case 59:  // [interval-single-value] ::= FLOATCONST
            if (reductionBody.size() == 1) {
                IdentifierToken* const token = dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (token != 0)
                    pushNode(new PlotFunctionTagValues(token));
            }
            break;
    } // switch (productionID)
}

ParserAction* PlotFunctionParser::action(const int stateID, const int symbolID) const {
    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 10:  // +
                    return new ParserActionTransition(2);

                case 9:  // -
                    return new ParserActionTransition(3);

                case 22:  // EMPTY
                    return new ParserActionTransition(4);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 1:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 10:  // +
                    return new ParserActionTransition(2);

                case 9:  // -
                    return new ParserActionTransition(3);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 2:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 3:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 9:  // -
                    return new ParserActionTransition(22);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 4:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(4);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 5:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(36);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 6:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionReduce(38);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 7:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(35);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 8:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(37);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 9:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(24);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 10:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(30);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 11:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(31);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 12:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionAccept();

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 13:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 14:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 15:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(34);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 16:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(1);

                case 10:  // +
                    return new ParserActionTransition(27);

                case 9:  // -
                    return new ParserActionTransition(28);

                case 7:  // :
                    return new ParserActionTransition(29);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 17:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(2);

                case 8:  // ;
                    return new ParserActionTransition(30);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 18:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                    return new ParserActionReduce(5);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 19:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(28);

                case 13:  // ^
                    return new ParserActionTransition(31);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 20:
            switch (symbolID) {
                case 10:  // +
                    return new ParserActionTransition(27);

                case 9:  // -
                    return new ParserActionTransition(28);

                case 2:  // )
                    return new ParserActionTransition(32);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 21:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(18);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 22:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 23:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(19);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 24:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 25:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 26:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(33);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 27:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 9:  // -
                    return new ParserActionTransition(36);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 28:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 10:  // +
                    return new ParserActionTransition(38);

                case 9:  // -
                    return new ParserActionTransition(39);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 29:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(41);

                case 9:  // -
                    return new ParserActionTransition(42);

                case 3:  // [
                    return new ParserActionTransition(43);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 30:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 10:  // +
                    return new ParserActionTransition(2);

                case 9:  // -
                    return new ParserActionTransition(3);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 31:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 32:
            switch (symbolID) {
                case -1:  // $END$
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(32);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 33:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(20);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 34:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(26);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 35:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(27);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 36:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 37:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(16);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 38:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 39:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(1);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(5);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 18:  // INTCONST
                    return new ParserActionTransition(7);

                case 20:  // VARIABLE
                    return new ParserActionTransition(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 40:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(17);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 41:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 10:  // +
                    return new ParserActionTransition(54);

                case 9:  // -
                    return new ParserActionTransition(55);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 42:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(41);

                case 9:  // -
                    return new ParserActionTransition(42);

                case 3:  // [
                    return new ParserActionTransition(43);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 43:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 10:  // +
                    return new ParserActionTransition(54);

                case 9:  // -
                    return new ParserActionTransition(55);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 44:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                    return new ParserActionReduce(7);

                case 14:  // |
                    return new ParserActionTransition(68);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 45:
            switch (symbolID) {
                case 6:  // ,
                    return new ParserActionTransition(69);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 46:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                case 14:  // |
                    return new ParserActionReduce(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 47:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(3);

                case 10:  // +
                    return new ParserActionTransition(27);

                case 9:  // -
                    return new ParserActionTransition(28);

                case 7:  // :
                    return new ParserActionTransition(29);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 48:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                    return new ParserActionReduce(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 49:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(29);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 50:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(21);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 51:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(22);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 52:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // +
                case 9:  // -
                case 7:  // :
                case 2:  // )
                    return new ParserActionReduce(23);

                case 11:  // *
                    return new ParserActionTransition(24);

                case 12:  // /
                    return new ParserActionTransition(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 53:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 10:  // +
                    return new ParserActionTransition(54);

                case 9:  // -
                    return new ParserActionTransition(55);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 54:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 55:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 9:  // -
                    return new ParserActionTransition(72);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 56:
            switch (symbolID) {
                case 6:  // ,
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(59);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 57:
            switch (symbolID) {
                case 6:  // ,
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(58);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 58:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 59:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(47);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 60:
            switch (symbolID) {
                case 6:  // ,
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(53);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 61:
            switch (symbolID) {
                case 6:  // ,
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(54);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 62:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(48);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 63:
            switch (symbolID) {
                case 6:  // ,
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(57);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 64:
            switch (symbolID) {
                case 6:  // ,
                    return new ParserActionReduce(12);

                case 10:  // +
                    return new ParserActionTransition(77);

                case 9:  // -
                    return new ParserActionTransition(78);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 65:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(51);

                case 13:  // ^
                    return new ParserActionTransition(79);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 66:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                case 14:  // |
                    return new ParserActionReduce(10);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 67:
            switch (symbolID) {
                case 6:  // ,
                    return new ParserActionReduce(13);

                case 10:  // +
                    return new ParserActionTransition(77);

                case 9:  // -
                    return new ParserActionTransition(78);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 68:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(41);

                case 9:  // -
                    return new ParserActionTransition(42);

                case 3:  // [
                    return new ParserActionTransition(43);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 69:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 10:  // +
                    return new ParserActionTransition(54);

                case 9:  // -
                    return new ParserActionTransition(55);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 70:
            switch (symbolID) {
                case 10:  // +
                    return new ParserActionTransition(77);

                case 9:  // -
                    return new ParserActionTransition(78);

                case 2:  // )
                    return new ParserActionTransition(83);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 71:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(41);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 72:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 73:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(42);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 74:
            switch (symbolID) {
                case 6:  // ,
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(56);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 75:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 76:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 77:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 9:  // -
                    return new ParserActionTransition(87);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 78:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                case 10:  // +
                    return new ParserActionTransition(89);

                case 9:  // -
                    return new ParserActionTransition(90);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 79:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 80:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                case 14:  // |
                    return new ParserActionReduce(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 81:
            switch (symbolID) {
                case 10:  // +
                    return new ParserActionTransition(77);

                case 9:  // -
                    return new ParserActionTransition(78);

                case 2:  // )
                    return new ParserActionTransition(93);

                case 4:  // ]
                    return new ParserActionTransition(94);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 82:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                case 14:  // |
                    return new ParserActionReduce(11);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 83:
            switch (symbolID) {
                case 6:  // ,
                case 13:  // ^
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(55);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 84:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(43);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 85:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(49);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 86:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(50);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 87:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 88:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(39);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 89:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 90:
            switch (symbolID) {
                case 1:  // (
                    return new ParserActionTransition(53);

                case 19:  // FLOATCONST
                    return new ParserActionTransition(56);

                case 18:  // INTCONST
                    return new ParserActionTransition(57);

                case 17:  // FUNCTION
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 91:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(40);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 92:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 11:  // *
                case 12:  // /
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(52);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 93:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                case 14:  // |
                    return new ParserActionReduce(14);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 94:
            switch (symbolID) {
                case -1:  // $END$
                case 8:  // ;
                case 14:  // |
                    return new ParserActionReduce(15);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 95:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(44);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 96:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(45);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        case 97:
            switch (symbolID) {
                case 6:  // ,
                case 10:  // +
                case 9:  // -
                case 2:  // )
                case 4:  // ]
                    return new ParserActionReduce(46);

                case 11:  // *
                    return new ParserActionTransition(75);

                case 12:  // /
                    return new ParserActionTransition(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID)
            break;

        default:
            return new ParserActionError();

    } // switch (stateID)

    return 0;
}

int PlotFunctionParser::gotoState(const int stateID, const int symbolID) const {
    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 23:  // [expression]
                    return 12;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 24:  // [term]
                    return 16;

                case 25:  // [term-list]
                    return 17;

                case 26:  // [termwithinterval]
                    return 18;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 9;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 1:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 24:  // [term]
                    return 20;

                case 32:  // [addend]
                    return 9;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 2:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 21;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 3:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 23;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 14:
            switch (symbolID) {
                case 36:  // [bracket-term]
                    return 26;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 22:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 33;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 24:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 33:  // [factor]
                    return 34;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 25:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 33:  // [factor]
                    return 35;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 27:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 37;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 28:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 40;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 29:
            switch (symbolID) {
                case 27:  // [interval]
                    return 44;

                case 29:  // [leftrange]
                    return 45;

                case 28:  // [range]
                    return 46;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 30:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 24:  // [term]
                    return 47;

                case 26:  // [termwithinterval]
                    return 48;

                case 32:  // [addend]
                    return 9;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 31:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 33:  // [factor]
                    return 49;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 36:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 50;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 38:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 51;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 39:
            switch (symbolID) {
                case 35:  // [atom]
                    return 10;

                case 36:  // [bracket-term]
                    return 11;

                case 33:  // [factor]
                    return 13;

                case 37:  // [function]
                    return 14;

                case 38:  // [single-value]
                    return 15;

                case 34:  // [variable]
                    return 19;

                case 32:  // [addend]
                    return 52;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 41:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 39:  // [interval-addend]
                    return 59;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 31:  // [interval-term]
                    return 64;

                case 41:  // [interval-variable]
                    return 65;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 42:
            switch (symbolID) {
                case 29:  // [leftrange]
                    return 45;

                case 28:  // [range]
                    return 66;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 43:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 39:  // [interval-addend]
                    return 59;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 31:  // [interval-term]
                    return 67;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 53:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 39:  // [interval-addend]
                    return 59;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 31:  // [interval-term]
                    return 70;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 54:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 71;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 55:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 73;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 58:
            switch (symbolID) {
                case 43:  // [interval-bracket-term]
                    return 74;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 68:
            switch (symbolID) {
                case 29:  // [leftrange]
                    return 45;

                case 28:  // [range]
                    return 80;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 69:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 39:  // [interval-addend]
                    return 59;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 31:  // [interval-term]
                    return 81;

                case 30:  // [rightrange]
                    return 82;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 72:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 75:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 40:  // [interval-factor]
                    return 85;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 76:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 40:  // [interval-factor]
                    return 86;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 77:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 88;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 78:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 91;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 79:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 40:  // [interval-factor]
                    return 92;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 87:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 95;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 89:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 96;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

        case 90:
            switch (symbolID) {
                case 37:  // [function]
                    return 58;

                case 42:  // [interval-atom]
                    return 60;

                case 43:  // [interval-bracket-term]
                    return 61;

                case 40:  // [interval-factor]
                    return 62;

                case 44:  // [interval-single-value]
                    return 63;

                case 41:  // [interval-variable]
                    return 65;

                case 39:  // [interval-addend]
                    return 97;

                default:
                    return -1;
            }   // switch (symbolID)
            break;

    } // switch (stateID)

    return -1;
}

ProductionStub* PlotFunctionParser::findProduction(const int productionID) const {
    switch (productionID) {
        case 0:  // [$START$] ::= [expression]
            {
                int bodyIDs[] = {23};
                return new ProductionStub(0, bodyIDs, 1);
            }
        case 1:  // [expression] ::= [term]
            {
                int bodyIDs[] = {24};
                return new ProductionStub(23, bodyIDs, 1);
            }
        case 2:  // [expression] ::= [term-list]
            {
                int bodyIDs[] = {25};
                return new ProductionStub(23, bodyIDs, 1);
            }
        case 3:  // [expression] ::= [term-list] ; [term]
            {
                int bodyIDs[] = {25, 8, 24};
                return new ProductionStub(23, bodyIDs, 3);
            }
        case 4:  // [expression] ::= EMPTY
            {
                int bodyIDs[] = {22};
                return new ProductionStub(23, bodyIDs, 1);
            }
        case 5:  // [term-list] ::= [termwithinterval]
            {
                int bodyIDs[] = {26};
                return new ProductionStub(25, bodyIDs, 1);
            }
        case 6:  // [term-list] ::= [term-list] ; [termwithinterval]
            {
                int bodyIDs[] = {25, 8, 26};
                return new ProductionStub(25, bodyIDs, 3);
            }
        case 7:  // [termwithinterval] ::= [term] : [interval]
            {
                int bodyIDs[] = {24, 7, 27};
                return new ProductionStub(26, bodyIDs, 3);
            }
        case 8:  // [interval] ::= [range]
            {
                int bodyIDs[] = {28};
                return new ProductionStub(27, bodyIDs, 1);
            }
        case 9:  // [interval] ::= [interval] | [range]
            {
                int bodyIDs[] = {27, 14, 28};
                return new ProductionStub(27, bodyIDs, 3);
            }
        case 10:  // [range] ::= - [range]
            {
                int bodyIDs[] = {9, 28};
                return new ProductionStub(28, bodyIDs, 2);
            }
        case 11:  // [range] ::= [leftrange] , [rightrange]
            {
                int bodyIDs[] = {29, 6, 30};
                return new ProductionStub(28, bodyIDs, 3);
            }
        case 12:  // [leftrange] ::= ( [interval-term]
            {
                int bodyIDs[] = {1, 31};
                return new ProductionStub(29, bodyIDs, 2);
            }
        case 13:  // [leftrange] ::= [ [interval-term]
            {
                int bodyIDs[] = {3, 31};
                return new ProductionStub(29, bodyIDs, 2);
            }
        case 14:  // [rightrange] ::= [interval-term] )
            {
                int bodyIDs[] = {31, 2};
                return new ProductionStub(30, bodyIDs, 2);
            }
        case 15:  // [rightrange] ::= [interval-term] ]
            {
                int bodyIDs[] = {31, 4};
                return new ProductionStub(30, bodyIDs, 2);
            }
        case 16:  // [term] ::= [term] + [addend]
            {
                int bodyIDs[] = {24, 10, 32};
                return new ProductionStub(24, bodyIDs, 3);
            }
        case 17:  // [term] ::= [term] - [addend]
            {
                int bodyIDs[] = {24, 9, 32};
                return new ProductionStub(24, bodyIDs, 3);
            }
        case 18:  // [term] ::= + [addend]
            {
                int bodyIDs[] = {10, 32};
                return new ProductionStub(24, bodyIDs, 2);
            }
        case 19:  // [term] ::= - [addend]
            {
                int bodyIDs[] = {9, 32};
                return new ProductionStub(24, bodyIDs, 2);
            }
        case 20:  // [term] ::= - - [addend]
            {
                int bodyIDs[] = {9, 9, 32};
                return new ProductionStub(24, bodyIDs, 3);
            }
        case 21:  // [term] ::= [term] + - [addend]
            {
                int bodyIDs[] = {24, 10, 9, 32};
                return new ProductionStub(24, bodyIDs, 4);
            }
        case 22:  // [term] ::= [term] - + [addend]
            {
                int bodyIDs[] = {24, 9, 10, 32};
                return new ProductionStub(24, bodyIDs, 4);
            }
        case 23:  // [term] ::= [term] - - [addend]
            {
                int bodyIDs[] = {24, 9, 9, 32};
                return new ProductionStub(24, bodyIDs, 4);
            }
        case 24:  // [term] ::= [addend]
            {
                int bodyIDs[] = {32};
                return new ProductionStub(24, bodyIDs, 1);
            }
        case 25:  // [addend] ::= [factor]
            {
                int bodyIDs[] = {33};
                return new ProductionStub(32, bodyIDs, 1);
            }
        case 26:  // [addend] ::= [addend] * [factor]
            {
                int bodyIDs[] = {32, 11, 33};
                return new ProductionStub(32, bodyIDs, 3);
            }
        case 27:  // [addend] ::= [addend] / [factor]
            {
                int bodyIDs[] = {32, 12, 33};
                return new ProductionStub(32, bodyIDs, 3);
            }
        case 28:  // [factor] ::= [variable]
            {
                int bodyIDs[] = {34};
                return new ProductionStub(33, bodyIDs, 1);
            }
        case 29:  // [factor] ::= [variable] ^ [factor]
            {
                int bodyIDs[] = {34, 13, 33};
                return new ProductionStub(33, bodyIDs, 3);
            }
        case 30:  // [variable] ::= [atom]
            {
                int bodyIDs[] = {35};
                return new ProductionStub(34, bodyIDs, 1);
            }
        case 31:  // [variable] ::= [bracket-term]
            {
                int bodyIDs[] = {36};
                return new ProductionStub(34, bodyIDs, 1);
            }
        case 32:  // [bracket-term] ::= ( [term] )
            {
                int bodyIDs[] = {1, 24, 2};
                return new ProductionStub(36, bodyIDs, 3);
            }
        case 33:  // [atom] ::= [function] [bracket-term]
            {
                int bodyIDs[] = {37, 36};
                return new ProductionStub(35, bodyIDs, 2);
            }
        case 34:  // [atom] ::= [single-value]
            {
                int bodyIDs[] = {38};
                return new ProductionStub(35, bodyIDs, 1);
            }
        case 35:  // [single-value] ::= INTCONST
            {
                int bodyIDs[] = {18};
                return new ProductionStub(38, bodyIDs, 1);
            }
        case 36:  // [single-value] ::= FLOATCONST
            {
                int bodyIDs[] = {19};
                return new ProductionStub(38, bodyIDs, 1);
            }
        case 37:  // [single-value] ::= VARIABLE
            {
                int bodyIDs[] = {20};
                return new ProductionStub(38, bodyIDs, 1);
            }
        case 38:  // [function] ::= FUNCTION
            {
                int bodyIDs[] = {17};
                return new ProductionStub(37, bodyIDs, 1);
            }
        case 39:  // [interval-term] ::= [interval-term] + [interval-addend]
            {
                int bodyIDs[] = {31, 10, 39};
                return new ProductionStub(31, bodyIDs, 3);
            }
        case 40:  // [interval-term] ::= [interval-term] - [interval-addend]
            {
                int bodyIDs[] = {31, 9, 39};
                return new ProductionStub(31, bodyIDs, 3);
            }
        case 41:  // [interval-term] ::= + [interval-addend]
            {
                int bodyIDs[] = {10, 39};
                return new ProductionStub(31, bodyIDs, 2);
            }
        case 42:  // [interval-term] ::= - [interval-addend]
            {
                int bodyIDs[] = {9, 39};
                return new ProductionStub(31, bodyIDs, 2);
            }
        case 43:  // [interval-term] ::= - - [interval-addend]
            {
                int bodyIDs[] = {9, 9, 39};
                return new ProductionStub(31, bodyIDs, 3);
            }
        case 44:  // [interval-term] ::= [interval-term] + - [interval-addend]
            {
                int bodyIDs[] = {31, 10, 9, 39};
                return new ProductionStub(31, bodyIDs, 4);
            }
        case 45:  // [interval-term] ::= [interval-term] - + [interval-addend]
            {
                int bodyIDs[] = {31, 9, 10, 39};
                return new ProductionStub(31, bodyIDs, 4);
            }
        case 46:  // [interval-term] ::= [interval-term] - - [interval-addend]
            {
                int bodyIDs[] = {31, 9, 9, 39};
                return new ProductionStub(31, bodyIDs, 4);
            }
        case 47:  // [interval-term] ::= [interval-addend]
            {
                int bodyIDs[] = {39};
                return new ProductionStub(31, bodyIDs, 1);
            }
        case 48:  // [interval-addend] ::= [interval-factor]
            {
                int bodyIDs[] = {40};
                return new ProductionStub(39, bodyIDs, 1);
            }
        case 49:  // [interval-addend] ::= [interval-addend] * [interval-factor]
            {
                int bodyIDs[] = {39, 11, 40};
                return new ProductionStub(39, bodyIDs, 3);
            }
        case 50:  // [interval-addend] ::= [interval-addend] / [interval-factor]
            {
                int bodyIDs[] = {39, 12, 40};
                return new ProductionStub(39, bodyIDs, 3);
            }
        case 51:  // [interval-factor] ::= [interval-variable]
            {
                int bodyIDs[] = {41};
                return new ProductionStub(40, bodyIDs, 1);
            }
        case 52:  // [interval-factor] ::= [interval-variable] ^ [interval-factor]
            {
                int bodyIDs[] = {41, 13, 40};
                return new ProductionStub(40, bodyIDs, 3);
            }
        case 53:  // [interval-variable] ::= [interval-atom]
            {
                int bodyIDs[] = {42};
                return new ProductionStub(41, bodyIDs, 1);
            }
        case 54:  // [interval-variable] ::= [interval-bracket-term]
            {
                int bodyIDs[] = {43};
                return new ProductionStub(41, bodyIDs, 1);
            }
        case 55:  // [interval-bracket-term] ::= ( [interval-term] )
            {
                int bodyIDs[] = {1, 31, 2};
                return new ProductionStub(43, bodyIDs, 3);
            }
        case 56:  // [interval-atom] ::= [function] [interval-bracket-term]
            {
                int bodyIDs[] = {37, 43};
                return new ProductionStub(42, bodyIDs, 2);
            }
        case 57:  // [interval-atom] ::= [interval-single-value]
            {
                int bodyIDs[] = {44};
                return new ProductionStub(42, bodyIDs, 1);
            }
        case 58:  // [interval-single-value] ::= INTCONST
            {
                int bodyIDs[] = {18};
                return new ProductionStub(44, bodyIDs, 1);
            }
        case 59:  // [interval-single-value] ::= FLOATCONST
            {
                int bodyIDs[] = {19};
                return new ProductionStub(44, bodyIDs, 1);
            }
    } // switch (productionID)

    return 0;
}

std::string PlotFunctionParser::symbolID2String(const int symbolID) const {
    switch (symbolID) {
        case 0: return "[$START$]";
        case 1: return "(";
        case 2: return ")";
        case 3: return "[";
        case 4: return "]";
        case 5: return ".";
        case 6: return ",";
        case 7: return ":";
        case 8: return ";";
        case 9: return "-";
        case 10: return "+";
        case 11: return "*";
        case 12: return "/";
        case 13: return "^";
        case 14: return "|";
        case 15: return "RANGE_TERM";
        case 16: return "FUNCTION-TERM";
        case 17: return "FUNCTION";
        case 18: return "INTCONST";
        case 19: return "FLOATCONST";
        case 20: return "VARIABLE";
        case 22: return "EMPTY";
        case 23: return "[expression]";
        case 24: return "[term]";
        case 25: return "[term-list]";
        case 26: return "[termwithinterval]";
        case 27: return "[interval]";
        case 28: return "[range]";
        case 29: return "[leftrange]";
        case 30: return "[rightrange]";
        case 31: return "[interval-term]";
        case 32: return "[addend]";
        case 33: return "[factor]";
        case 34: return "[variable]";
        case 35: return "[atom]";
        case 36: return "[bracket-term]";
        case 37: return "[function]";
        case 38: return "[single-value]";
        case 39: return "[interval-addend]";
        case 40: return "[interval-factor]";
        case 41: return "[interval-variable]";
        case 42: return "[interval-atom]";
        case 43: return "[interval-bracket-term]";
        case 44: return "[interval-single-value]";
        case 45: return "$END$";
    }  // switch (symbolID)
    return "";
}

}   // namespace glslparser

}   // namespace voreen
