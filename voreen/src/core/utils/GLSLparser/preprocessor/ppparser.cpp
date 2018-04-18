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

#include "voreen/core/utils/GLSLparser/preprocessor/ppparser.h"

#include "voreen/core/utils/GLSLparser/preprocessor/pplexer.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppmacro.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppstatement.h"

namespace voreen {

namespace glslparser {

PreprocessorParser::PreprocessorParser(std::istream* is) 
    : Parser(new PreprocessorLexer(is, PreprocessorTerminals().getKeywords()), true)
    , tokens_()
    //, parseExpression_(false)
{
}

PreprocessorParser::PreprocessorParser(const std::string& fileName) 
    : Parser(new PreprocessorLexer(fileName, PreprocessorTerminals().getKeywords()), true)
    , tokens_()
    //, parseExpression_(false)
{
}

PreprocessorParser::PreprocessorParser(const std::list<Token*>& tokens)
    : Parser(0, false)
    , tokens_(tokens)
    //, parseExpression_(true)
{
}

PreprocessorParser::~PreprocessorParser() {
    for (; ! tokens_.empty(); tokens_.pop_front())
        delete tokens_.front();
}

// protected methods
//

bool PreprocessorParser::isProxySymbol(const int symbolID, const int /*original*/) const {
    // The token "TOKEN" is a wildcard for all tokens. Any production
    // which expects the terminal token, has to except any token!
    //
    if (symbolID == PreprocessorTerminals::ID_TOKEN)
        return true;

    return false;
}

Token* PreprocessorParser::nextToken() {
    if (! Parser::hasLexer()) {
        if (! tokens_.empty()) {
            Token* t = tokens_.front();
            tokens_.pop_front();
            return t;
        } else
            return 0;
    }

    return Parser::nextToken();
}

void PreprocessorParser::expandParseTree(const int productionID,
                                   const std::vector<Parser::ParserSymbol*>& reductionBody)
{
    switch (productionID) {
        case 0:  // [$START$] ::= [program]
        case 1:  // [program] ::= [particle-list]
        case 3:  // [expression-parsing] ::= pirate [constant-expression]
        case 4:  // [particle-list] ::= [particle]
        case 6:  // [particle] ::= [directive]
            break;

        case 5:  // [particle-list] ::= [particle-list] [particle]
            {
                ParseTreeNode* next = popNode();
                ParseTreeNode* list = popNode();
                if ((list != 0) && (next != 0))
                    list->addChild(next);
                pushNode(list);
            }
            break;

        case 7:  // [particle] ::= text
            if (reductionBody.size() == 1) {
                TextToken* const t
                    = dynamic_cast<TextToken* const>(reductionBody[0]->getToken());
                if (t != 0)
                    pushNode(new TextNode(*t));
            }
            break;

        case 8:  // [particle] ::= newline
        case 9:  // [directive] ::= [def-directive]
        case 10:  // [directive] ::= [undef-directive]
        case 11:  // [directive] ::= [error-directive]
        case 12:  // [directive] ::= [extension-directive]
        case 13:  // [directive] ::= [include-directive]
        case 14:  // [directive] ::= [line-directive]
        case 15:  // [directive] ::= [null-directive]
        case 16:  // [directive] ::= [pragma-directive]
        case 17:  // [directive] ::= [version-directive]
        case 18:  // [directive] ::= [conditional-directive]
            break;

        case 19:  // [def-directive] ::= #define identifier newline
            if (reductionBody.size() == 3) {
                const IdentifierToken* const id =
                    dynamic_cast<const IdentifierToken* const>(reductionBody[1]->getToken());

                if (id != 0)
                    pushNode(new DefineDirective(id->getValue(), new TokenList(), false));
            }
            break;

        case 20:  // [def-directive] ::= #define identifier [token-list] newline
            if (reductionBody.size() == 4) {
                TokenList* list = popNode<TokenList>();
                const IdentifierToken* const id =
                    dynamic_cast<const IdentifierToken* const>(reductionBody[2]->getToken());

                if ((list != 0) && (id != 0))
                    pushNode(new DefineDirective(id->getValue(), list, false));
            }
            break;

        case 21:  // [def-directive] ::= #define identifier ( ) newline
            if ((reductionBody.size() == 5)) {
                const IdentifierToken* const id =
                    dynamic_cast<const IdentifierToken* const>(reductionBody[3]->getToken());

                if (id != 0)
                    pushNode(new DefineDirective(id->getValue(), new TokenList(), true));
            }
            break;

        case 22:  // [def-directive] ::= #define identifier ( ) [token-list] newline
            if ((reductionBody.size() == 6)) {
                TokenList* list = popNode<TokenList>();
                const IdentifierToken* const id =
                    dynamic_cast<const IdentifierToken* const>(reductionBody[4]->getToken());

                if ((list != 0) && (id != 0))
                    pushNode(new DefineDirective(id->toString(), list, true));
            }
            break;

        case 23:  // [def-directive] ::= #define identifier ( [formals-list] ) newline
            if (reductionBody.size() == 6) {
                IdentifierList* formals = popNode<IdentifierList>();
                const IdentifierToken* const id =
                    dynamic_cast<const IdentifierToken* const>(reductionBody[4]->getToken());

                if ((formals != 0) && (id != 0))
                    pushNode(new DefineDirective(id->getValue(), new TokenList(), true, formals));
            }
            break;

        case 24:  // [def-directive] ::= #define identifier ( [formals-list] ) [token-list] newline
            if (reductionBody.size() == 7) {
                TokenList* list = popNode<TokenList>();
                IdentifierList* formals = popNode<IdentifierList>();
                const IdentifierToken* const id =
                    dynamic_cast<const IdentifierToken* const>(reductionBody[5]->getToken());

                if ((list != 0) && (id != 0))
                    pushNode(new DefineDirective(id->getValue(), list, true, formals));
            }
            break;

        case 25:  // [token-list] ::= token
            if (reductionBody.size() == 1) {
                Token* const t = reductionBody[0]->getToken();

                 if (t != 0)
                     pushNode(new TokenList(t));
            }
            break;

        case 26:  // [token-list] ::= [token-list] token
            if (reductionBody.size() == 2) {
                TokenList* list = popNode<TokenList>();
                Token* const t = reductionBody[0]->getToken();

                 if ((list != 0) && (t != 0)) {
                     list->addToken(t);
                     pushNode(list);
                 }
            }
            break;

        case 27:  // [formals-list] ::= identifier
            if (reductionBody.size() == 1) {
                const IdentifierToken* const id =
                    dynamic_cast<const IdentifierToken* const>(reductionBody[0]->getToken());

                if (id != 0)
                    pushNode(new IdentifierList(id->getValue()));
            }
            break;

        case 28:  // [formals-list] ::= [formals-list] , identifier
            if (reductionBody.size() == 3) {
                IdentifierList* const list = popNode<IdentifierList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if ((list != 0) && (id != 0)) {
                    list->addIdentifier(id->getValue());
                    pushNode(list);
                }
            }
            break;

        case 29:  // [error-directive] ::= #error newline
            if (reductionBody.size() == 2)
                pushNode(new ErrorDirective(new TokenList()));
            break;

        case 30:  // [error-directive] ::= #error [token-list] newline
            if (reductionBody.size() == 3) {
                TokenList* const list = popNode<TokenList>();

                if (list != 0)
                    pushNode(new ErrorDirective(list));
            }
            break;

        case 31:  // [extension-directive] ::= #extension newline
            if (reductionBody.size() == 2)
                pushNode(new ExtensionDirective(new TokenList()));
            break;

        case 32:  // [extension-directive] ::= #extension [token-list] newline
            if (reductionBody.size() == 3) {
                TokenList* const list = popNode<TokenList>();

                if (list != 0)
                    pushNode(new ExtensionDirective(list));
            }
            break;

        case 33:  // [include-directive] ::= #include " string " newline
            if (reductionBody.size() == 5) {
                StringToken* const str =
                    dynamic_cast<StringToken* const>(reductionBody[2]->getToken());
                Token* const includeToken = reductionBody[4]->getToken();

                if ((str != 0) && (includeToken != 0))
                    pushNode(new IncludeDirective(str->getValue()));
            }
            break;

        case 34:  // [line-directive] ::= #line int-constant newline
            if (reductionBody.size() == 3) {
                ConstantToken* const intConst =
                    dynamic_cast<ConstantToken* const>(reductionBody[1]->getToken());

                if ((intConst != 0) && (intConst->getType() == ConstantToken::TYPE_INT))
                    pushNode(new LineDirective(*intConst));
            }
            break;

        case 35:  // [line-directive] ::= #line int-constant " string " newline
            if (reductionBody.size() == 6) {
                StringToken* const str =
                    dynamic_cast<StringToken* const>(reductionBody[2]->getToken());
                ConstantToken* const intConst =
                    dynamic_cast<ConstantToken* const>(reductionBody[4]->getToken());

                if ((str != 0) && (intConst != 0) && (intConst->getType() == ConstantToken::TYPE_INT))
                    pushNode(new LineDirective(*intConst, *str));
            }
            break;

        case 36:  // [null-directive] ::= # newline
            break;

        case 37:  // [pragma-directive] ::= #pragma [token-list] newline
            if (reductionBody.size() == 3) {
                TokenList* const list = popNode<TokenList>();

                if (list != 0)
                    pushNode(new PragmaDirective(list));
            }
            break;

        case 38:  // [undef-directive] ::= #undef identifier newline
            if (reductionBody.size() == 2) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());
                Token* const undefToken = reductionBody[0]->getToken();

                if ((undefToken != 0) && (id != 0))
                    pushNode(new UndefineDirective(id->getValue()));
            }
            break;

        case 39:  // [version-directive] ::= #version int-constant newline
            if (reductionBody.size() == 3) {
                ConstantToken* const intConst =
                    dynamic_cast<ConstantToken* const>(reductionBody[1]->getToken());

                if (intConst != 0)
                    pushNode(new VersionDirective(*intConst));
            }
            break;

        case 40:  // [version-directive] ::= #version int-constant identifier newline
            if (reductionBody.size() == 4) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[1]->getToken());
                ConstantToken* const intConst =
                    dynamic_cast<ConstantToken* const>(reductionBody[2]->getToken());

                if ((id != 0) && (intConst != 0))
                    pushNode(new VersionDirective(*intConst, *id));
            }
            break;

        case 41:  // [conditional-directive] ::= [if-part] [endif-directive]
        case 42:  // [conditional-directive] ::= [if-part] [else-part] [endif-directive]
            break;

        case 43:  // [if-part] ::= [if-directive] [particle-list]
            if (reductionBody.size() == 2) {
                ParseTreeNode* const nodes = popNode();
                IfDirective* const ifdir = popNode<IfDirective>();

                if ((ifdir != 0) && (nodes != 0)) {
                    ifdir->setTrue(nodes);
                    pushNode(ifdir);
                }
            }
            break;

        case 44:  // [if-part] ::= [ifdef-directive] [particle-list]
            if (reductionBody.size() == 2) {
                ParseTreeNode* const nodes = popNode();
                IfdefDirective* const ifdef = popNode<IfdefDirective>();

                if ((ifdef != 0) && (nodes != 0)) {
                    ifdef->setTrue(nodes);
                    pushNode(ifdef);
                }
            }
            break;

        case 45:  // [if-part] ::= [ifndef-directive] [particle-list]
            if (reductionBody.size() == 2) {
                ParseTreeNode* const nodes = popNode();
                IfndefDirective* const ifndef = popNode<IfndefDirective>();

                if ((ifndef != 0) && (nodes != 0)) {
                    ifndef->setTrue(nodes);
                    pushNode(ifndef);
                }
            }
            break;

        case 46:  // [elif-parts] ::= [elif-part]
            break;

        case 47:  // [elif-parts] ::= [elif-parts] [elif-part]
            if (reductionBody.size() == 2) {
                IfDirective* const ifdir = popNode<IfDirective>();
                IfDirective* const list = popNode<IfDirective>();

                if ((ifdir != 0) && (list != 0)) {
                    list->setFalse(ifdir);
                    pushNode(ifdir);
                }
            }
            break;

        case 48:  // [elif-part] ::= [elif-directive] [particle-list]
             if (reductionBody.size() == 2) {
                ParseTreeNode* const list = popNode();
                IfDirective* const elif = popNode<IfDirective>();

                if ((elif != 0) && (list != 0)) {
                    elif->setTrue(list);
                    pushNode(elif);
                }
            }
            break;

        case 49:  // [else-part] ::= [elif-parts]
            break;

        case 50:  // [else-part] ::= [elif-parts] [else-directive] [particle-list]
            if (reductionBody.size() == 3) {
                ParseTreeNode* const nodes = popNode();
                ConditionalDirective* const elifs = popNode<ConditionalDirective>();
                ConditionalDirective* const cond = popNode<ConditionalDirective>();

                if ((cond != 0) && (elifs != 0) && (nodes != 0)) {
                    elifs->setFalse(nodes);
                    cond->setFalse(elifs);
                    pushNode(cond);
                }
            }
            break;

        case 51:  // [else-part] ::= [else-directive] [particle-list]
            if (reductionBody.size() == 2) {
                ParseTreeNode* const nodes = popNode();
                ConditionalDirective* const cond = popNode<ConditionalDirective>();

                if ((cond != 0) && (nodes != 0)) {
                    cond->setFalse(nodes);
                    pushNode(cond);
                }
            }
            break;

        case 52:  // [if-directive] ::= #if [constant-expression] newline
            if (reductionBody.size() == 3) {
                Expression* const exp = popNode<Expression>();

                if (exp != 0)
                    pushNode(new IfDirective(exp, 0));
            }
            break;

        case 53:  // [ifdef-directive] ::= #ifdef identifier newline
            if (reductionBody.size() == 3) {
                IdentifierToken* const id
                    = dynamic_cast<IdentifierToken* const>(reductionBody[1]->getToken());

                if (id != 0)
                    pushNode(new IfdefDirective(id, 0));
            }
            break;

        case 54:  // [ifndef-directive] ::= #ifndef identifier newline
            if (reductionBody.size() == 3) {
                IdentifierToken* const id
                    = dynamic_cast<IdentifierToken* const>(reductionBody[1]->getToken());

                if (id != 0)
                    pushNode(new IfndefDirective(id, 0));
            }
            break;

        case 55:  // [else-directive] ::= #else newline
            break;

        case 56:  // [elif-directive] ::= #elif [constant-expression] newline
            if (reductionBody.size() == 3) {
                Expression* const exp = popNode<Expression>();
                if (exp != 0)
                    pushNode(new IfDirective(exp, 0));
            }
            break;

        case 57:  // [endif-directive] ::= #endif newline
            break;

        case 58:  // [constant-expression] ::= [logical-or-expression]
            break;

        case 59:  // [macro-evaluation] ::= identifier
            if (reductionBody.size() == 1) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (id != 0)
                    pushNode(new Macro(id, false));
            }
            break;

        case 60:  // [macro-evaluation] ::= identifier ( )
            if (reductionBody.size() == 3) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());

                if (id != 0)
                    pushNode(new Macro(id, true));
            }
            break;

        case 61:  // [macro-evaluation] ::= identifier ( [parameter-list] )
            if (reductionBody.size() == 4) {
                ExpressionList* const list = popNode<ExpressionList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[3]->getToken());

                if ((list != 0) && (id != 0))
                    pushNode(new Macro(id, true, list));
            }
            break;

        case 62:  // [parameter-list] ::= [constant-expression]
            if (reductionBody.size() == 1) {
                Expression* const exp = popNode<Expression>();
                if (exp != 0)
                    pushNode(new ExpressionList(exp));
            }
            break;

        case 63:  // [parameter-list] ::= [parameter-list] , [constant-expression]
            if (reductionBody.size() == 3) {
                Expression* const exp = popNode<Expression>();
                ExpressionList* const list = popNode<ExpressionList>();

                if ((list != 0) && (exp != 0)) {
                    list->addExpression(exp);
                    pushNode(list);
                }
            }
            break;

        case 64:  // [defined-operator] ::= defined ( identifier )
            if (reductionBody.size() == 4) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[1]->getToken());

                if (id != 0)
                    pushNode( new DefinedOperator(id) );
            }
            break;

        case 65:  // [defined-operator] ::= defined identifier
            if (reductionBody.size() == 2) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (id != 0)
                    pushNode( new DefinedOperator(id) );
            }
            break;

        case 66:  // [logical-or-expression] ::= [logical-and-expression]
        case 68:  // [logical-and-expression] ::= [inclusive-or-expression]
        case 70:  // [inclusive-or-expression] ::= [exclusive-or-expression]
        case 72:  // [exclusive-or-expression] ::= [and-expression]
        case 74:  // [and-expression] ::= [equality-expression]
        case 76:  // [equality-expression] ::= [relational-expression]
        case 79:  // [relational-expression] ::= [shift-expression]
        case 84:  // [shift-expression] ::= [additive-expression]
        case 87:  // [additive-expression] ::= [multiplicative-expression]
        case 90:  // [multiplicative-expression] ::= [unary-expression]
        case 94:  // [unary-expression] ::= [primary-expression]
        case 99:  // [primary-expression] ::= [defined-operator]
        case 100:  // [primary-expression] ::= [macro-evaluation]
            break;

        case 67:  // [logical-or-expression] ::= [logical-or-expression] || [logical-and-expression]
        case 69:  // [logical-and-expression] ::= [logical-and-expression] && [inclusive-or-expression]
        case 77:  // [equality-expression] ::= [equality-expression] == [relational-expression]
        case 78:  // [equality-expression] ::= [equality-expression] != [relational-expression]
        case 80:  // [relational-expression] ::= [relational-expression] < [shift-expression]
        case 81:  // [relational-expression] ::= [relational-expression] > [shift-expression]
        case 82:  // [relational-expression] ::= [relational-expression] <= [shift-expression]
        case 83:  // [relational-expression] ::= [relational-expression] >= [shift-expression]
            if (reductionBody.size() == 3) {
                Expression* const rhs = popNode<Expression>();
                Expression* const lhs = popNode<Expression>();
                Token* const op = reductionBody[1]->getToken();

                if ((lhs != 0) && (rhs != 0) && (op != 0))
                    pushNode(new LogicalBinaryExpression(op, lhs, rhs));
            }
            break;

        case 71:  // [inclusive-or-expression] ::= [inclusive-or-expression] | [exclusive-or-expression]
        case 73:  // [exclusive-or-expression] ::= [exclusive-or-expression] ^ [and-expression]
        case 75:  // [and-expression] ::= [and-expression] & [equality-expression]
        case 85:  // [shift-expression] ::= [shift-expression] << [additive-expression]
        case 86:  // [shift-expression] ::= [shift-expression] >> [additive-expression]
        case 88:  // [additive-expression] ::= [additive-expression] + [multiplicative-expression]
        case 89:  // [additive-expression] ::= [additive-expression] - [multiplicative-expression]
        case 91:  // [multiplicative-expression] ::= [multiplicative-expression] * [unary-expression]
        case 92:  // [multiplicative-expression] ::= [multiplicative-expression] / [unary-expression]
        case 93:  // [multiplicative-expression] ::= [multiplicative-expression] % [unary-expression]
            if (reductionBody.size() == 3) {
                Expression* const rhs = popNode<Expression>();
                Expression* const lhs = popNode<Expression>();
                Token* const op = reductionBody[1]->getToken();

                if ((lhs != 0) && (rhs != 0) && (op != 0))
                    pushNode(new BinaryExpression(op, lhs, rhs));
            }
            break;

        case 95:  // [unary-expression] ::= + [unary-expression]
        case 96:  // [unary-expression] ::= - [unary-expression]
        case 97:  // [unary-expression] ::= ! [unary-expression]
        case 98:  // [unary-expression] ::= ~ [unary-expression]
            if (reductionBody.size() == 2) {
                Expression* const unexp = popNode<Expression>();
                Token* const op = reductionBody[1]->getToken();

                if ((unexp != 0) && (op != 0))
                    pushNode(new UnaryExpression(op, unexp));
            }
            break;

        case 101:  // [primary-expression] ::= int-constant
            if (reductionBody.size() == 1) {
                ConstantToken* const intConst =
                    dynamic_cast<ConstantToken* const>(reductionBody[0]->getToken());

                if ((intConst != 0) && (intConst->getType() == ConstantToken::TYPE_INT))
                    pushNode(new IntConstant(intConst));
            }
            break;

        case 102:  // [primary-expression] ::= ( [constant-expression] )
            if (reductionBody.size() == 3) {
                Token* const lparen = reductionBody[2]->getToken();
                Token* const rparen = reductionBody[0]->getToken();
                Expression* const expr = popNode<Expression>();

                if ((expr != 0) && (lparen != 0) && (rparen != 0))
                    pushNode(new ParenthesisExpression(lparen, expr));
            }
            break;
    } // switch (productionID
}   // expandParseTree();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THE FOLLOWING ~7000 LINES CONTAIN GENERATED CODE AND SHOULD *NEVER* BE    //
// MODIFIED !!!                                                              //
//                                                                           //
// The only non-generated modifications concern the treatment of the         //
// terminal symbol TOKEN which is actually a wildcard for ALL tokens emitted //
// by the lexical anlysator (lexer) driving the parser.                      //
//                                                                           //
//                                                                    (dirk) //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

ParserAction* PreprocessorParser::action(const int stateID, const int symbolID) const {
    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 50:  // pirate
                    return new ParserActionTransition(29);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 1:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(35);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 2:
            switch (symbolID) {
                case 44:  // identifier
                    return new ParserActionTransition(36);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 3:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(37);

                //case 42:  // token
                default:
                    return new ParserActionTransition(38);
            }   // switch (symbolID
            break;

        case 4:
            switch (symbolID) {
                //case 42:  // token
                default:
                    return new ParserActionTransition(38);

                case 47:  // newline
                    return new ParserActionTransition(40);
            }   // switch (symbolID
            break;

        case 5:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 6:
            switch (symbolID) {
                case 44:  // identifier
                    return new ParserActionTransition(65);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 7:
            switch (symbolID) {
                case 44:  // identifier
                    return new ParserActionTransition(66);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 8:
            switch (symbolID) {
                case 21:  // "
                    return new ParserActionTransition(67);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 9:
            switch (symbolID) {
                case 45:  // int-constant
                    return new ParserActionTransition(68);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 10:
            return new ParserActionTransition(38);
            /*switch (symbolID) {
                case 42:  // token
                    return new ParserActionTransition(38);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;*/

        case 11:
            switch (symbolID) {
                case 44:  // identifier
                    return new ParserActionTransition(70);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 12:
            switch (symbolID) {
                case 45:  // int-constant
                    return new ParserActionTransition(71);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 13:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(18);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 14:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 15:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 16:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(11);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 17:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(2);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 18:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(12);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 19:
            switch (symbolID) {
                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 20:
            switch (symbolID) {
                case 5:  // #elif
                    return new ParserActionTransition(73);

                case 6:  // #else
                    return new ParserActionTransition(74);

                case 7:  // #endif
                    return new ParserActionTransition(75);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 21:
            switch (symbolID) {
                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 22:
            switch (symbolID) {
                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 23:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(13);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 24:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(14);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 25:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 26:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(15);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 27:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(4);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 28:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(1);

                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 29:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 30:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(16);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 31:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionAccept();

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 32:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(7);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 33:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(10);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 34:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(17);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 35:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(36);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 36:
            switch (symbolID) {
                //case 42:  // token
                default:
                    return new ParserActionTransition(38);

                case 18:  // (
                    return new ParserActionTransition(86);

                case 47:  // newline
                    return new ParserActionTransition(87);
            }   // switch (symbolID
            break;

        case 37:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(29);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 38:
            return new ParserActionReduce(25);
            /*switch (symbolID) {
                case 47:  // newline
                case 42:  // token
                    return new ParserActionReduce(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;*/

        case 39:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(89);

                //case 42:  // token
                default:
                    return new ParserActionTransition(90);
            }   // switch (symbolID
            break;

        case 40:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(31);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 41:
            switch (symbolID) {
                //case 42:  // token
                default:
                    return new ParserActionTransition(90);

                case 47:  // newline
                    return new ParserActionTransition(91);
            }   // switch (symbolID
            break;

        case 42:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 43:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 44:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 45:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 46:
            switch (symbolID) {
                case 47:  // newline
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(84);

                case 22:  // +
                    return new ParserActionTransition(96);

                case 23:  // -
                    return new ParserActionTransition(97);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 47:
            switch (symbolID) {
                case 47:  // newline
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(72);

                case 37:  // &
                    return new ParserActionTransition(98);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 48:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 49:
            switch (symbolID) {
                case 18:  // (
                    return new ParserActionTransition(100);

                case 44:  // identifier
                    return new ParserActionTransition(101);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 50:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 51:
            switch (symbolID) {
                case 47:  // newline
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(74);

                case 36:  // !=
                    return new ParserActionTransition(102);

                case 35:  // ==
                    return new ParserActionTransition(103);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 52:
            switch (symbolID) {
                case 47:  // newline
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(70);

                case 38:  // ^
                    return new ParserActionTransition(104);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 53:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(59);

                case 18:  // (
                    return new ParserActionTransition(105);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 54:
            switch (symbolID) {
                case 47:  // newline
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(68);

                case 39:  // |
                    return new ParserActionTransition(106);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 55:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(101);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 56:
            switch (symbolID) {
                case 47:  // newline
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(66);

                case 40:  // &&
                    return new ParserActionTransition(107);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 57:
            switch (symbolID) {
                case 47:  // newline
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(58);

                case 41:  // ||
                    return new ParserActionTransition(108);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 58:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(100);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 59:
            switch (symbolID) {
                case 47:  // newline
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(87);

                case 28:  // %
                    return new ParserActionTransition(109);

                case 26:  // *
                    return new ParserActionTransition(110);

                case 27:  // /
                    return new ParserActionTransition(111);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 60:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(94);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 61:
            switch (symbolID) {
                case 47:  // newline
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(76);

                case 31:  // <
                    return new ParserActionTransition(112);

                case 32:  // <=
                    return new ParserActionTransition(113);

                case 34:  // >
                    return new ParserActionTransition(114);

                case 33:  // >=
                    return new ParserActionTransition(115);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 62:
            switch (symbolID) {
                case 47:  // newline
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(79);

                case 29:  // <<
                    return new ParserActionTransition(116);

                case 30:  // >>
                    return new ParserActionTransition(117);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 63:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(90);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 64:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 65:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(119);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 66:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(120);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 67:
            switch (symbolID) {
                case 46:  // string
                    return new ParserActionTransition(121);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 68:
            switch (symbolID) {
                case 21:  // "
                    return new ParserActionTransition(122);

                case 47:  // newline
                    return new ParserActionTransition(123);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 69:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(124);

                //case 42:  // token
                default:
                    return new ParserActionTransition(90);
            }   // switch (symbolID
            break;

        case 70:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(125);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 71:
            switch (symbolID) {
                case 44:  // identifier
                    return new ParserActionTransition(126);

                case 47:  // newline
                    return new ParserActionTransition(127);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 72:
            switch (symbolID) {
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(43);

                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 73:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 74:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(129);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 75:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(130);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 76:
            switch (symbolID) {
                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 77:
            switch (symbolID) {
                case 7:  // #endif
                case 6:  // #else
                case 5:  // #elif
                    return new ParserActionReduce(46);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 78:
            switch (symbolID) {
                case 7:  // #endif
                    return new ParserActionReduce(49);

                case 5:  // #elif
                    return new ParserActionTransition(73);

                case 6:  // #else
                    return new ParserActionTransition(74);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 79:
            switch (symbolID) {
                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 80:
            switch (symbolID) {
                case 7:  // #endif
                    return new ParserActionTransition(75);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 81:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(41);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 82:
            switch (symbolID) {
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(44);

                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 83:
            switch (symbolID) {
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(45);

                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 84:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(5);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 85:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(3);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 86:
            switch (symbolID) {
                case 19:  // )
                    return new ParserActionTransition(136);

                case 44:  // identifier
                    return new ParserActionTransition(138);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 87:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(19);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 88:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(139);

                //case 42:  // token
                default:
                    return new ParserActionTransition(90);
            }   // switch (symbolID
            break;

        case 89:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(30);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 90:
            return new ParserActionReduce(26);
            /*switch (symbolID) {
                case 47:  // newline
                case 42:  // token
                    return new ParserActionReduce(26);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;*/

        case 91:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(32);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 92:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(97);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 93:
            switch (symbolID) {
                case 19:  // )
                    return new ParserActionTransition(140);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 94:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(95);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 95:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(96);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 96:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 97:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 98:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 99:
            switch (symbolID) {
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                    return new ParserActionReduce(52);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 100:
            switch (symbolID) {
                case 44:  // identifier
                    return new ParserActionTransition(144);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 101:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(65);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 102:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 103:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 104:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 105:
            switch (symbolID) {
                case 19:  // )
                    return new ParserActionTransition(148);

                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 106:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 107:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 108:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 109:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 110:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 111:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 112:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 113:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 114:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 115:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 116:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 117:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 118:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(98);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 119:
            switch (symbolID) {
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                    return new ParserActionReduce(53);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 120:
            switch (symbolID) {
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                    return new ParserActionReduce(54);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 121:
            switch (symbolID) {
                case 21:  // "
                    return new ParserActionTransition(163);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 122:
            switch (symbolID) {
                case 46:  // string
                    return new ParserActionTransition(164);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 123:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(34);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 124:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(37);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 125:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(38);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 126:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(165);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 127:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(39);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 128:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(166);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 129:
            switch (symbolID) {
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                    return new ParserActionReduce(55);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 130:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(57);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 131:
            switch (symbolID) {
                case 7:  // #endif
                case 6:  // #else
                case 5:  // #elif
                    return new ParserActionReduce(48);

                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 132:
            switch (symbolID) {
                case 7:  // #endif
                case 6:  // #else
                case 5:  // #elif
                    return new ParserActionReduce(47);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 133:
            switch (symbolID) {
                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 134:
            switch (symbolID) {
                case 7:  // #endif
                    return new ParserActionReduce(51);

                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 135:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 136:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(168);

                //case 42:  // token
                default:
                    return new ParserActionTransition(38);
            }   // switch (symbolID
            break;

        case 137:
            switch (symbolID) {
                case 19:  // )
                    return new ParserActionTransition(170);

                case 20:  // ,
                    return new ParserActionTransition(171);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 138:
            switch (symbolID) {
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(27);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 139:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(20);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 140:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(102);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 141:
            switch (symbolID) {
                case 47:  // newline
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(88);

                case 28:  // %
                    return new ParserActionTransition(109);

                case 26:  // *
                    return new ParserActionTransition(110);

                case 27:  // /
                    return new ParserActionTransition(111);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 142:
            switch (symbolID) {
                case 47:  // newline
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(89);

                case 28:  // %
                    return new ParserActionTransition(109);

                case 26:  // *
                    return new ParserActionTransition(110);

                case 27:  // /
                    return new ParserActionTransition(111);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 143:
            switch (symbolID) {
                case 47:  // newline
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(75);

                case 36:  // !=
                    return new ParserActionTransition(102);

                case 35:  // ==
                    return new ParserActionTransition(103);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 144:
            switch (symbolID) {
                case 19:  // )
                    return new ParserActionTransition(172);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 145:
            switch (symbolID) {
                case 47:  // newline
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(78);

                case 31:  // <
                    return new ParserActionTransition(112);

                case 32:  // <=
                    return new ParserActionTransition(113);

                case 34:  // >
                    return new ParserActionTransition(114);

                case 33:  // >=
                    return new ParserActionTransition(115);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 146:
            switch (symbolID) {
                case 47:  // newline
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(77);

                case 31:  // <
                    return new ParserActionTransition(112);

                case 32:  // <=
                    return new ParserActionTransition(113);

                case 34:  // >
                    return new ParserActionTransition(114);

                case 33:  // >=
                    return new ParserActionTransition(115);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 147:
            switch (symbolID) {
                case 47:  // newline
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(73);

                case 37:  // &
                    return new ParserActionTransition(98);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 148:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(60);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 149:
            switch (symbolID) {
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(62);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 150:
            switch (symbolID) {
                case 19:  // )
                    return new ParserActionTransition(173);

                case 20:  // ,
                    return new ParserActionTransition(174);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 151:
            switch (symbolID) {
                case 47:  // newline
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(71);

                case 38:  // ^
                    return new ParserActionTransition(104);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 152:
            switch (symbolID) {
                case 47:  // newline
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(69);

                case 39:  // |
                    return new ParserActionTransition(106);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 153:
            switch (symbolID) {
                case 47:  // newline
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(67);

                case 40:  // &&
                    return new ParserActionTransition(107);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 154:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(93);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 155:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(91);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 156:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(92);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 157:
            switch (symbolID) {
                case 47:  // newline
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(80);

                case 29:  // <<
                    return new ParserActionTransition(116);

                case 30:  // >>
                    return new ParserActionTransition(117);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 158:
            switch (symbolID) {
                case 47:  // newline
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(82);

                case 29:  // <<
                    return new ParserActionTransition(116);

                case 30:  // >>
                    return new ParserActionTransition(117);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 159:
            switch (symbolID) {
                case 47:  // newline
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(81);

                case 29:  // <<
                    return new ParserActionTransition(116);

                case 30:  // >>
                    return new ParserActionTransition(117);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 160:
            switch (symbolID) {
                case 47:  // newline
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(83);

                case 29:  // <<
                    return new ParserActionTransition(116);

                case 30:  // >>
                    return new ParserActionTransition(117);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 161:
            switch (symbolID) {
                case 47:  // newline
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(85);

                case 22:  // +
                    return new ParserActionTransition(96);

                case 23:  // -
                    return new ParserActionTransition(97);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 162:
            switch (symbolID) {
                case 47:  // newline
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(86);

                case 22:  // +
                    return new ParserActionTransition(96);

                case 23:  // -
                    return new ParserActionTransition(97);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 163:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(175);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 164:
            switch (symbolID) {
                case 21:  // "
                    return new ParserActionTransition(176);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 165:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(40);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 166:
            switch (symbolID) {
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                    return new ParserActionReduce(56);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 167:
            switch (symbolID) {
                case 7:  // #endif
                    return new ParserActionReduce(50);

                case 1:  // #
                    return new ParserActionTransition(1);

                case 11:  // #pragma
                    return new ParserActionTransition(10);

                case 3:  // #undef
                    return new ParserActionTransition(11);

                case 12:  // #version
                    return new ParserActionTransition(12);

                case 2:  // #define
                    return new ParserActionTransition(2);

                case 47:  // newline
                    return new ParserActionTransition(25);

                case 10:  // #error
                    return new ParserActionTransition(3);

                case 43:  // text
                    return new ParserActionTransition(32);

                case 14:  // #extension
                    return new ParserActionTransition(4);

                case 4:  // #if
                    return new ParserActionTransition(5);

                case 8:  // #ifdef
                    return new ParserActionTransition(6);

                case 9:  // #ifndef
                    return new ParserActionTransition(7);

                case 15:  // #include
                    return new ParserActionTransition(8);

                case 13:  // #line
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 168:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(21);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 169:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(177);

                //case 42:  // token
                default:
                    return new ParserActionTransition(90);
            }   // switch (symbolID
            break;

        case 170:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(178);

                //case 42:  // token
                default:
                    return new ParserActionTransition(38);
            }   // switch (symbolID
            break;

        case 171:
            switch (symbolID) {
                case 44:  // identifier
                    return new ParserActionTransition(180);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 172:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 173:
            switch (symbolID) {
                case 47:  // newline
                case 26:  // *
                case 27:  // /
                case 28:  // %
                case 22:  // +
                case 23:  // -
                case 29:  // <<
                case 30:  // >>
                case 31:  // <
                case 34:  // >
                case 32:  // <=
                case 33:  // >=
                case 35:  // ==
                case 36:  // !=
                case 37:  // &
                case 38:  // ^
                case 39:  // |
                case 40:  // &&
                case 41:  // ||
                case -1:  // $END$
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(61);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 174:
            switch (symbolID) {
                case 25:  // !
                    return new ParserActionTransition(42);

                case 18:  // (
                    return new ParserActionTransition(43);

                case 22:  // +
                    return new ParserActionTransition(44);

                case 23:  // -
                    return new ParserActionTransition(45);

                case 16:  // defined
                    return new ParserActionTransition(49);

                case 44:  // identifier
                    return new ParserActionTransition(53);

                case 45:  // int-constant
                    return new ParserActionTransition(55);

                case 24:  // ~
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 175:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(33);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 176:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(182);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 177:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(22);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 178:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(23);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 179:
            switch (symbolID) {
                case 47:  // newline
                    return new ParserActionTransition(183);

                //case 42:  // token
                default:
                    return new ParserActionTransition(90);
            }   // switch (symbolID
            break;

        case 180:
            switch (symbolID) {
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(28);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 181:
            switch (symbolID) {
                case 19:  // )
                case 20:  // ,
                    return new ParserActionReduce(63);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 182:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(35);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 183:
            switch (symbolID) {
                case -1:  // $END$
                case 10:  // #error
                case 14:  // #extension
                case 4:  // #if
                case 8:  // #ifdef
                case 9:  // #ifndef
                case 15:  // #include
                case 13:  // #line
                case 11:  // #pragma
                case 3:  // #undef
                case 12:  // #version
                case 47:  // newline
                case 43:  // text
                case 1:  // #
                case 2:  // #define
                case 7:  // #endif
                case 5:  // #elif
                case 6:  // #else
                    return new ParserActionReduce(24);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        default:
            return new ParserActionError();

    } // switch (stateID

    return 0;
}

int PreprocessorParser::gotoState(const int stateID, const int symbolID) const {
    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 53:  // [expression-parsing]
                    return 17;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 55:  // [particle]
                    return 27;

                case 52:  // [particle-list]
                    return 28;

                case 64:  // [pragma-directive]
                    return 30;

                case 51:  // [program]
                    return 31;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 3:
            switch (symbolID) {
                case 67:  // [token-list]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 4:
            switch (symbolID) {
                case 67:  // [token-list]
                    return 41;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 5:
            switch (symbolID) {
                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 54:  // [constant-expression]
                    return 48;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 83:  // [inclusive-or-expression]
                    return 54;

                case 82:  // [logical-and-expression]
                    return 56;

                case 79:  // [logical-or-expression]
                    return 57;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 10:
            switch (symbolID) {
                case 67:  // [token-list]
                    return 69;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 19:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 55:  // [particle]
                    return 27;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 52:  // [particle-list]
                    return 72;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 20:
            switch (symbolID) {
                case 77:  // [elif-directive]
                    return 76;

                case 75:  // [elif-part]
                    return 77;

                case 76:  // [elif-parts]
                    return 78;

                case 78:  // [else-directive]
                    return 79;

                case 71:  // [else-part]
                    return 80;

                case 70:  // [endif-directive]
                    return 81;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 21:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 55:  // [particle]
                    return 27;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 52:  // [particle-list]
                    return 82;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 22:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 55:  // [particle]
                    return 27;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 52:  // [particle-list]
                    return 83;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 28:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 55:  // [particle]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 29:
            switch (symbolID) {
                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 83:  // [inclusive-or-expression]
                    return 54;

                case 82:  // [logical-and-expression]
                    return 56;

                case 79:  // [logical-or-expression]
                    return 57;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                case 54:  // [constant-expression]
                    return 85;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 36:
            switch (symbolID) {
                case 67:  // [token-list]
                    return 88;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 42:
            switch (symbolID) {
                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 92;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 43:
            switch (symbolID) {
                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 83:  // [inclusive-or-expression]
                    return 54;

                case 82:  // [logical-and-expression]
                    return 56;

                case 79:  // [logical-or-expression]
                    return 57;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                case 54:  // [constant-expression]
                    return 93;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 44:
            switch (symbolID) {
                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 94;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 45:
            switch (symbolID) {
                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 95;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 64:
            switch (symbolID) {
                case 91:  // [unary-expression]
                    return 118;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 72:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 55:  // [particle]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 73:
            switch (symbolID) {
                case 54:  // [constant-expression]
                    return 128;

                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 83:  // [inclusive-or-expression]
                    return 54;

                case 82:  // [logical-and-expression]
                    return 56;

                case 79:  // [logical-or-expression]
                    return 57;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 76:
            switch (symbolID) {
                case 52:  // [particle-list]
                    return 131;

                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 55:  // [particle]
                    return 27;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 78:
            switch (symbolID) {
                case 75:  // [elif-part]
                    return 132;

                case 78:  // [else-directive]
                    return 133;

                case 77:  // [elif-directive]
                    return 76;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 79:
            switch (symbolID) {
                case 52:  // [particle-list]
                    return 134;

                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 55:  // [particle]
                    return 27;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 80:
            switch (symbolID) {
                case 70:  // [endif-directive]
                    return 135;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 82:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 55:  // [particle]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 83:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 55:  // [particle]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 86:
            switch (symbolID) {
                case 68:  // [formals-list]
                    return 137;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 96:
            switch (symbolID) {
                case 90:  // [multiplicative-expression]
                    return 141;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 97:
            switch (symbolID) {
                case 90:  // [multiplicative-expression]
                    return 142;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 98:
            switch (symbolID) {
                case 86:  // [equality-expression]
                    return 143;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 102:
            switch (symbolID) {
                case 87:  // [relational-expression]
                    return 145;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 103:
            switch (symbolID) {
                case 87:  // [relational-expression]
                    return 146;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 104:
            switch (symbolID) {
                case 85:  // [and-expression]
                    return 147;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 105:
            switch (symbolID) {
                case 54:  // [constant-expression]
                    return 149;

                case 80:  // [parameter-list]
                    return 150;

                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 83:  // [inclusive-or-expression]
                    return 54;

                case 82:  // [logical-and-expression]
                    return 56;

                case 79:  // [logical-or-expression]
                    return 57;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 106:
            switch (symbolID) {
                case 84:  // [exclusive-or-expression]
                    return 151;

                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 107:
            switch (symbolID) {
                case 83:  // [inclusive-or-expression]
                    return 152;

                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 108:
            switch (symbolID) {
                case 82:  // [logical-and-expression]
                    return 153;

                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 83:  // [inclusive-or-expression]
                    return 54;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 109:
            switch (symbolID) {
                case 91:  // [unary-expression]
                    return 154;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 110:
            switch (symbolID) {
                case 91:  // [unary-expression]
                    return 155;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 111:
            switch (symbolID) {
                case 91:  // [unary-expression]
                    return 156;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 92:  // [primary-expression]
                    return 60;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 112:
            switch (symbolID) {
                case 88:  // [shift-expression]
                    return 157;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 113:
            switch (symbolID) {
                case 88:  // [shift-expression]
                    return 158;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 114:
            switch (symbolID) {
                case 88:  // [shift-expression]
                    return 159;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 115:
            switch (symbolID) {
                case 88:  // [shift-expression]
                    return 160;

                case 89:  // [additive-expression]
                    return 46;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 116:
            switch (symbolID) {
                case 89:  // [additive-expression]
                    return 161;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 117:
            switch (symbolID) {
                case 89:  // [additive-expression]
                    return 162;

                case 81:  // [defined-operator]
                    return 50;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 131:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 55:  // [particle]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 133:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 52:  // [particle-list]
                    return 167;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 55:  // [particle]
                    return 27;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 134:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 55:  // [particle]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 136:
            switch (symbolID) {
                case 67:  // [token-list]
                    return 169;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 167:
            switch (symbolID) {
                case 66:  // [conditional-directive]
                    return 13;

                case 57:  // [def-directive]
                    return 14;

                case 56:  // [directive]
                    return 15;

                case 59:  // [error-directive]
                    return 16;

                case 60:  // [extension-directive]
                    return 18;

                case 72:  // [if-directive]
                    return 19;

                case 69:  // [if-part]
                    return 20;

                case 73:  // [ifdef-directive]
                    return 21;

                case 74:  // [ifndef-directive]
                    return 22;

                case 61:  // [include-directive]
                    return 23;

                case 62:  // [line-directive]
                    return 24;

                case 63:  // [null-directive]
                    return 26;

                case 64:  // [pragma-directive]
                    return 30;

                case 58:  // [undef-directive]
                    return 33;

                case 65:  // [version-directive]
                    return 34;

                case 55:  // [particle]
                    return 84;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 170:
            switch (symbolID) {
                case 67:  // [token-list]
                    return 179;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 174:
            switch (symbolID) {
                case 54:  // [constant-expression]
                    return 181;

                case 89:  // [additive-expression]
                    return 46;

                case 85:  // [and-expression]
                    return 47;

                case 81:  // [defined-operator]
                    return 50;

                case 86:  // [equality-expression]
                    return 51;

                case 84:  // [exclusive-or-expression]
                    return 52;

                case 83:  // [inclusive-or-expression]
                    return 54;

                case 82:  // [logical-and-expression]
                    return 56;

                case 79:  // [logical-or-expression]
                    return 57;

                case 48:  // [macro-evaluation]
                    return 58;

                case 90:  // [multiplicative-expression]
                    return 59;

                case 92:  // [primary-expression]
                    return 60;

                case 87:  // [relational-expression]
                    return 61;

                case 88:  // [shift-expression]
                    return 62;

                case 91:  // [unary-expression]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

    } // switch (stateID

    return -1;
}

ProductionStub* PreprocessorParser::findProduction(const int productionID) const {
    switch (productionID) {
        case 0:  // [$START$] ::= [program]
            {
                int bodyIDs[] = {51};
                return new ProductionStub(0, bodyIDs, 1);
            }
        case 1:  // [program] ::= [particle-list]
            {
                int bodyIDs[] = {52};
                return new ProductionStub(51, bodyIDs, 1);
            }
        case 2:  // [program] ::= [expression-parsing]
            {
                int bodyIDs[] = {53};
                return new ProductionStub(51, bodyIDs, 1);
            }
        case 3:  // [expression-parsing] ::= pirate [constant-expression]
            {
                int bodyIDs[] = {50, 54};
                return new ProductionStub(53, bodyIDs, 2);
            }
        case 4:  // [particle-list] ::= [particle]
            {
                int bodyIDs[] = {55};
                return new ProductionStub(52, bodyIDs, 1);
            }
        case 5:  // [particle-list] ::= [particle-list] [particle]
            {
                int bodyIDs[] = {52, 55};
                return new ProductionStub(52, bodyIDs, 2);
            }
        case 6:  // [particle] ::= [directive]
            {
                int bodyIDs[] = {56};
                return new ProductionStub(55, bodyIDs, 1);
            }
        case 7:  // [particle] ::= text
            {
                int bodyIDs[] = {43};
                return new ProductionStub(55, bodyIDs, 1);
            }
        case 8:  // [particle] ::= newline
            {
                int bodyIDs[] = {47};
                return new ProductionStub(55, bodyIDs, 1);
            }
        case 9:  // [directive] ::= [def-directive]
            {
                int bodyIDs[] = {57};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 10:  // [directive] ::= [undef-directive]
            {
                int bodyIDs[] = {58};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 11:  // [directive] ::= [error-directive]
            {
                int bodyIDs[] = {59};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 12:  // [directive] ::= [extension-directive]
            {
                int bodyIDs[] = {60};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 13:  // [directive] ::= [include-directive]
            {
                int bodyIDs[] = {61};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 14:  // [directive] ::= [line-directive]
            {
                int bodyIDs[] = {62};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 15:  // [directive] ::= [null-directive]
            {
                int bodyIDs[] = {63};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 16:  // [directive] ::= [pragma-directive]
            {
                int bodyIDs[] = {64};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 17:  // [directive] ::= [version-directive]
            {
                int bodyIDs[] = {65};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 18:  // [directive] ::= [conditional-directive]
            {
                int bodyIDs[] = {66};
                return new ProductionStub(56, bodyIDs, 1);
            }
        case 19:  // [def-directive] ::= #define identifier newline
            {
                int bodyIDs[] = {2, 44, 47};
                return new ProductionStub(57, bodyIDs, 3);
            }
        case 20:  // [def-directive] ::= #define identifier [token-list] newline
            {
                int bodyIDs[] = {2, 44, 67, 47};
                return new ProductionStub(57, bodyIDs, 4);
            }
        case 21:  // [def-directive] ::= #define identifier ( ) newline
            {
                int bodyIDs[] = {2, 44, 18, 19, 47};
                return new ProductionStub(57, bodyIDs, 5);
            }
        case 22:  // [def-directive] ::= #define identifier ( ) [token-list] newline
            {
                int bodyIDs[] = {2, 44, 18, 19, 67, 47};
                return new ProductionStub(57, bodyIDs, 6);
            }
        case 23:  // [def-directive] ::= #define identifier ( [formals-list] ) newline
            {
                int bodyIDs[] = {2, 44, 18, 68, 19, 47};
                return new ProductionStub(57, bodyIDs, 6);
            }
        case 24:  // [def-directive] ::= #define identifier ( [formals-list] ) [token-list] newline
            {
                int bodyIDs[] = {2, 44, 18, 68, 19, 67, 47};
                return new ProductionStub(57, bodyIDs, 7);
            }
        case 25:  // [token-list] ::= token
            {
                int bodyIDs[] = {42};
                return new ProductionStub(67, bodyIDs, 1);
            }
        case 26:  // [token-list] ::= [token-list] token
            {
                int bodyIDs[] = {67, 42};
                return new ProductionStub(67, bodyIDs, 2);
            }
        case 27:  // [formals-list] ::= identifier
            {
                int bodyIDs[] = {44};
                return new ProductionStub(68, bodyIDs, 1);
            }
        case 28:  // [formals-list] ::= [formals-list] , identifier
            {
                int bodyIDs[] = {68, 20, 44};
                return new ProductionStub(68, bodyIDs, 3);
            }
        case 29:  // [error-directive] ::= #error newline
            {
                int bodyIDs[] = {10, 47};
                return new ProductionStub(59, bodyIDs, 2);
            }
        case 30:  // [error-directive] ::= #error [token-list] newline
            {
                int bodyIDs[] = {10, 67, 47};
                return new ProductionStub(59, bodyIDs, 3);
            }
        case 31:  // [extension-directive] ::= #extension newline
            {
                int bodyIDs[] = {14, 47};
                return new ProductionStub(60, bodyIDs, 2);
            }
        case 32:  // [extension-directive] ::= #extension [token-list] newline
            {
                int bodyIDs[] = {14, 67, 47};
                return new ProductionStub(60, bodyIDs, 3);
            }
        case 33:  // [include-directive] ::= #include " string " newline
            {
                int bodyIDs[] = {15, 21, 46, 21, 47};
                return new ProductionStub(61, bodyIDs, 5);
            }
        case 34:  // [line-directive] ::= #line int-constant newline
            {
                int bodyIDs[] = {13, 45, 47};
                return new ProductionStub(62, bodyIDs, 3);
            }
        case 35:  // [line-directive] ::= #line int-constant " string " newline
            {
                int bodyIDs[] = {13, 45, 21, 46, 21, 47};
                return new ProductionStub(62, bodyIDs, 6);
            }
        case 36:  // [null-directive] ::= # newline
            {
                int bodyIDs[] = {1, 47};
                return new ProductionStub(63, bodyIDs, 2);
            }
        case 37:  // [pragma-directive] ::= #pragma [token-list] newline
            {
                int bodyIDs[] = {11, 67, 47};
                return new ProductionStub(64, bodyIDs, 3);
            }
        case 38:  // [undef-directive] ::= #undef identifier newline
            {
                int bodyIDs[] = {3, 44, 47};
                return new ProductionStub(58, bodyIDs, 3);
            }
        case 39:  // [version-directive] ::= #version int-constant newline
            {
                int bodyIDs[] = {12, 45, 47};
                return new ProductionStub(65, bodyIDs, 3);
            }
        case 40:  // [version-directive] ::= #version int-constant identifier newline
            {
                int bodyIDs[] = {12, 45, 44, 47};
                return new ProductionStub(65, bodyIDs, 4);
            }
        case 41:  // [conditional-directive] ::= [if-part] [endif-directive]
            {
                int bodyIDs[] = {69, 70};
                return new ProductionStub(66, bodyIDs, 2);
            }
        case 42:  // [conditional-directive] ::= [if-part] [else-part] [endif-directive]
            {
                int bodyIDs[] = {69, 71, 70};
                return new ProductionStub(66, bodyIDs, 3);
            }
        case 43:  // [if-part] ::= [if-directive] [particle-list]
            {
                int bodyIDs[] = {72, 52};
                return new ProductionStub(69, bodyIDs, 2);
            }
        case 44:  // [if-part] ::= [ifdef-directive] [particle-list]
            {
                int bodyIDs[] = {73, 52};
                return new ProductionStub(69, bodyIDs, 2);
            }
        case 45:  // [if-part] ::= [ifndef-directive] [particle-list]
            {
                int bodyIDs[] = {74, 52};
                return new ProductionStub(69, bodyIDs, 2);
            }
        case 46:  // [elif-parts] ::= [elif-part]
            {
                int bodyIDs[] = {75};
                return new ProductionStub(76, bodyIDs, 1);
            }
        case 47:  // [elif-parts] ::= [elif-parts] [elif-part]
            {
                int bodyIDs[] = {76, 75};
                return new ProductionStub(76, bodyIDs, 2);
            }
        case 48:  // [elif-part] ::= [elif-directive] [particle-list]
            {
                int bodyIDs[] = {77, 52};
                return new ProductionStub(75, bodyIDs, 2);
            }
        case 49:  // [else-part] ::= [elif-parts]
            {
                int bodyIDs[] = {76};
                return new ProductionStub(71, bodyIDs, 1);
            }
        case 50:  // [else-part] ::= [elif-parts] [else-directive] [particle-list]
            {
                int bodyIDs[] = {76, 78, 52};
                return new ProductionStub(71, bodyIDs, 3);
            }
        case 51:  // [else-part] ::= [else-directive] [particle-list]
            {
                int bodyIDs[] = {78, 52};
                return new ProductionStub(71, bodyIDs, 2);
            }
        case 52:  // [if-directive] ::= #if [constant-expression] newline
            {
                int bodyIDs[] = {4, 54, 47};
                return new ProductionStub(72, bodyIDs, 3);
            }
        case 53:  // [ifdef-directive] ::= #ifdef identifier newline
            {
                int bodyIDs[] = {8, 44, 47};
                return new ProductionStub(73, bodyIDs, 3);
            }
        case 54:  // [ifndef-directive] ::= #ifndef identifier newline
            {
                int bodyIDs[] = {9, 44, 47};
                return new ProductionStub(74, bodyIDs, 3);
            }
        case 55:  // [else-directive] ::= #else newline
            {
                int bodyIDs[] = {6, 47};
                return new ProductionStub(78, bodyIDs, 2);
            }
        case 56:  // [elif-directive] ::= #elif [constant-expression] newline
            {
                int bodyIDs[] = {5, 54, 47};
                return new ProductionStub(77, bodyIDs, 3);
            }
        case 57:  // [endif-directive] ::= #endif newline
            {
                int bodyIDs[] = {7, 47};
                return new ProductionStub(70, bodyIDs, 2);
            }
        case 58:  // [constant-expression] ::= [logical-or-expression]
            {
                int bodyIDs[] = {79};
                return new ProductionStub(54, bodyIDs, 1);
            }
        case 59:  // [macro-evaluation] ::= identifier
            {
                int bodyIDs[] = {44};
                return new ProductionStub(48, bodyIDs, 1);
            }
        case 60:  // [macro-evaluation] ::= identifier ( )
            {
                int bodyIDs[] = {44, 18, 19};
                return new ProductionStub(48, bodyIDs, 3);
            }
        case 61:  // [macro-evaluation] ::= identifier ( [parameter-list] )
            {
                int bodyIDs[] = {44, 18, 80, 19};
                return new ProductionStub(48, bodyIDs, 4);
            }
        case 62:  // [parameter-list] ::= [constant-expression]
            {
                int bodyIDs[] = {54};
                return new ProductionStub(80, bodyIDs, 1);
            }
        case 63:  // [parameter-list] ::= [parameter-list] , [constant-expression]
            {
                int bodyIDs[] = {80, 20, 54};
                return new ProductionStub(80, bodyIDs, 3);
            }
        case 64:  // [defined-operator] ::= defined ( identifier )
            {
                int bodyIDs[] = {16, 18, 44, 19};
                return new ProductionStub(81, bodyIDs, 4);
            }
        case 65:  // [defined-operator] ::= defined identifier
            {
                int bodyIDs[] = {16, 44};
                return new ProductionStub(81, bodyIDs, 2);
            }
        case 66:  // [logical-or-expression] ::= [logical-and-expression]
            {
                int bodyIDs[] = {82};
                return new ProductionStub(79, bodyIDs, 1);
            }
        case 67:  // [logical-or-expression] ::= [logical-or-expression] || [logical-and-expression]
            {
                int bodyIDs[] = {79, 41, 82};
                return new ProductionStub(79, bodyIDs, 3);
            }
        case 68:  // [logical-and-expression] ::= [inclusive-or-expression]
            {
                int bodyIDs[] = {83};
                return new ProductionStub(82, bodyIDs, 1);
            }
        case 69:  // [logical-and-expression] ::= [logical-and-expression] && [inclusive-or-expression]
            {
                int bodyIDs[] = {82, 40, 83};
                return new ProductionStub(82, bodyIDs, 3);
            }
        case 70:  // [inclusive-or-expression] ::= [exclusive-or-expression]
            {
                int bodyIDs[] = {84};
                return new ProductionStub(83, bodyIDs, 1);
            }
        case 71:  // [inclusive-or-expression] ::= [inclusive-or-expression] | [exclusive-or-expression]
            {
                int bodyIDs[] = {83, 39, 84};
                return new ProductionStub(83, bodyIDs, 3);
            }
        case 72:  // [exclusive-or-expression] ::= [and-expression]
            {
                int bodyIDs[] = {85};
                return new ProductionStub(84, bodyIDs, 1);
            }
        case 73:  // [exclusive-or-expression] ::= [exclusive-or-expression] ^ [and-expression]
            {
                int bodyIDs[] = {84, 38, 85};
                return new ProductionStub(84, bodyIDs, 3);
            }
        case 74:  // [and-expression] ::= [equality-expression]
            {
                int bodyIDs[] = {86};
                return new ProductionStub(85, bodyIDs, 1);
            }
        case 75:  // [and-expression] ::= [and-expression] & [equality-expression]
            {
                int bodyIDs[] = {85, 37, 86};
                return new ProductionStub(85, bodyIDs, 3);
            }
        case 76:  // [equality-expression] ::= [relational-expression]
            {
                int bodyIDs[] = {87};
                return new ProductionStub(86, bodyIDs, 1);
            }
        case 77:  // [equality-expression] ::= [equality-expression] == [relational-expression]
            {
                int bodyIDs[] = {86, 35, 87};
                return new ProductionStub(86, bodyIDs, 3);
            }
        case 78:  // [equality-expression] ::= [equality-expression] != [relational-expression]
            {
                int bodyIDs[] = {86, 36, 87};
                return new ProductionStub(86, bodyIDs, 3);
            }
        case 79:  // [relational-expression] ::= [shift-expression]
            {
                int bodyIDs[] = {88};
                return new ProductionStub(87, bodyIDs, 1);
            }
        case 80:  // [relational-expression] ::= [relational-expression] < [shift-expression]
            {
                int bodyIDs[] = {87, 31, 88};
                return new ProductionStub(87, bodyIDs, 3);
            }
        case 81:  // [relational-expression] ::= [relational-expression] > [shift-expression]
            {
                int bodyIDs[] = {87, 34, 88};
                return new ProductionStub(87, bodyIDs, 3);
            }
        case 82:  // [relational-expression] ::= [relational-expression] <= [shift-expression]
            {
                int bodyIDs[] = {87, 32, 88};
                return new ProductionStub(87, bodyIDs, 3);
            }
        case 83:  // [relational-expression] ::= [relational-expression] >= [shift-expression]
            {
                int bodyIDs[] = {87, 33, 88};
                return new ProductionStub(87, bodyIDs, 3);
            }
        case 84:  // [shift-expression] ::= [additive-expression]
            {
                int bodyIDs[] = {89};
                return new ProductionStub(88, bodyIDs, 1);
            }
        case 85:  // [shift-expression] ::= [shift-expression] << [additive-expression]
            {
                int bodyIDs[] = {88, 29, 89};
                return new ProductionStub(88, bodyIDs, 3);
            }
        case 86:  // [shift-expression] ::= [shift-expression] >> [additive-expression]
            {
                int bodyIDs[] = {88, 30, 89};
                return new ProductionStub(88, bodyIDs, 3);
            }
        case 87:  // [additive-expression] ::= [multiplicative-expression]
            {
                int bodyIDs[] = {90};
                return new ProductionStub(89, bodyIDs, 1);
            }
        case 88:  // [additive-expression] ::= [additive-expression] + [multiplicative-expression]
            {
                int bodyIDs[] = {89, 22, 90};
                return new ProductionStub(89, bodyIDs, 3);
            }
        case 89:  // [additive-expression] ::= [additive-expression] - [multiplicative-expression]
            {
                int bodyIDs[] = {89, 23, 90};
                return new ProductionStub(89, bodyIDs, 3);
            }
        case 90:  // [multiplicative-expression] ::= [unary-expression]
            {
                int bodyIDs[] = {91};
                return new ProductionStub(90, bodyIDs, 1);
            }
        case 91:  // [multiplicative-expression] ::= [multiplicative-expression] * [unary-expression]
            {
                int bodyIDs[] = {90, 26, 91};
                return new ProductionStub(90, bodyIDs, 3);
            }
        case 92:  // [multiplicative-expression] ::= [multiplicative-expression] / [unary-expression]
            {
                int bodyIDs[] = {90, 27, 91};
                return new ProductionStub(90, bodyIDs, 3);
            }
        case 93:  // [multiplicative-expression] ::= [multiplicative-expression] % [unary-expression]
            {
                int bodyIDs[] = {90, 28, 91};
                return new ProductionStub(90, bodyIDs, 3);
            }
        case 94:  // [unary-expression] ::= [primary-expression]
            {
                int bodyIDs[] = {92};
                return new ProductionStub(91, bodyIDs, 1);
            }
        case 95:  // [unary-expression] ::= + [unary-expression]
            {
                int bodyIDs[] = {22, 91};
                return new ProductionStub(91, bodyIDs, 2);
            }
        case 96:  // [unary-expression] ::= - [unary-expression]
            {
                int bodyIDs[] = {23, 91};
                return new ProductionStub(91, bodyIDs, 2);
            }
        case 97:  // [unary-expression] ::= ! [unary-expression]
            {
                int bodyIDs[] = {25, 91};
                return new ProductionStub(91, bodyIDs, 2);
            }
        case 98:  // [unary-expression] ::= ~ [unary-expression]
            {
                int bodyIDs[] = {24, 91};
                return new ProductionStub(91, bodyIDs, 2);
            }
        case 99:  // [primary-expression] ::= [defined-operator]
            {
                int bodyIDs[] = {81};
                return new ProductionStub(92, bodyIDs, 1);
            }
        case 100:  // [primary-expression] ::= [macro-evaluation]
            {
                int bodyIDs[] = {48};
                return new ProductionStub(92, bodyIDs, 1);
            }
        case 101:  // [primary-expression] ::= int-constant
            {
                int bodyIDs[] = {45};
                return new ProductionStub(92, bodyIDs, 1);
            }
        case 102:  // [primary-expression] ::= ( [constant-expression] )
            {
                int bodyIDs[] = {18, 54, 19};
                return new ProductionStub(92, bodyIDs, 3);
            }
    } // switch (productionID

    return 0;
}

std::string PreprocessorParser::symbolID2String(const int symbolID) const {
    switch (symbolID) {
        case 0: return "[$START$]";
        case 1: return "#";
        case 2: return "#define";
        case 3: return "#undef";
        case 4: return "#if";
        case 5: return "#elif";
        case 6: return "#else";
        case 7: return "#endif";
        case 8: return "#ifdef";
        case 9: return "#ifndef";
        case 10: return "#error";
        case 11: return "#pragma";
        case 12: return "#version";
        case 13: return "#line";
        case 14: return "#extension";
        case 15: return "#include";
        case 16: return "defined";
        case 17: return "##";
        case 18: return "(";
        case 19: return ")";
        case 20: return ",";
        case 21: return "\"";
        case 22: return "+";
        case 23: return "-";
        case 24: return "~";
        case 25: return "!";
        case 26: return "*";
        case 27: return "/";
        case 28: return "%";
        case 29: return "<<";
        case 30: return ">>";
        case 31: return "<";
        case 32: return "<=";
        case 33: return ">=";
        case 34: return ">";
        case 35: return "==";
        case 36: return "!=";
        case 37: return "&";
        case 38: return "^";
        case 39: return "|";
        case 40: return "&&";
        case 41: return "||";
        case 42: return "token";
        case 43: return "text";
        case 44: return "identifier";
        case 45: return "int-constant";
        case 46: return "string";
        case 47: return "newline";
        case 48: return "[macro-evaluation]";
        case 50: return "pirate";
        case 51: return "[program]";
        case 52: return "[particle-list]";
        case 53: return "[expression-parsing]";
        case 54: return "[constant-expression]";
        case 55: return "[particle]";
        case 56: return "[directive]";
        case 57: return "[def-directive]";
        case 58: return "[undef-directive]";
        case 59: return "[error-directive]";
        case 60: return "[extension-directive]";
        case 61: return "[include-directive]";
        case 62: return "[line-directive]";
        case 63: return "[null-directive]";
        case 64: return "[pragma-directive]";
        case 65: return "[version-directive]";
        case 66: return "[conditional-directive]";
        case 67: return "[token-list]";
        case 68: return "[formals-list]";
        case 69: return "[if-part]";
        case 70: return "[endif-directive]";
        case 71: return "[else-part]";
        case 72: return "[if-directive]";
        case 73: return "[ifdef-directive]";
        case 74: return "[ifndef-directive]";
        case 75: return "[elif-part]";
        case 76: return "[elif-parts]";
        case 77: return "[elif-directive]";
        case 78: return "[else-directive]";
        case 79: return "[logical-or-expression]";
        case 80: return "[parameter-list]";
        case 81: return "[defined-operator]";
        case 82: return "[logical-and-expression]";
        case 83: return "[inclusive-or-expression]";
        case 84: return "[exclusive-or-expression]";
        case 85: return "[and-expression]";
        case 86: return "[equality-expression]";
        case 87: return "[relational-expression]";
        case 88: return "[shift-expression]";
        case 89: return "[additive-expression]";
        case 90: return "[multiplicative-expression]";
        case 91: return "[unary-expression]";
        case 92: return "[primary-expression]";
        case 93: return "$END$";
    }  // switch (symbolID
    return "";
}

}   // namespace glslparser

}   // namespace voreen
