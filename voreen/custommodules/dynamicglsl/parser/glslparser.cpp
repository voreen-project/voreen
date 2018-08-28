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

#include "glslparser.h"

#include "glsllexer.h"

#include "voreen/core/utils/GLSLparser/glsl/glslcondition.h"
#include "voreen/core/utils/GLSLparser/glsl/glslfunctioncall.h"
#include "voreen/core/utils/GLSLparser/glsl/glslfunctiondefinition.h"
#include "voreen/core/utils/GLSLparser/glsl/glslstructspecifier.h"
#include "voreen/core/utils/GLSLparser/glsl/glsltranslation.h"

namespace voreen {

namespace glslparser {

GLSLParser::GLSLParser(std::istream* const is) 
    : Parser(new GLSLLexer(is, GLSLTerminals().getKeywords()), true)
{
}

GLSLParser::GLSLParser(const std::string& fileName) 
    : Parser(new GLSLLexer(fileName, GLSLTerminals().getKeywords()), true)
{
}

GLSLParser::~GLSLParser() {
}

// protected methods
//

bool GLSLParser::isProxySymbol(const int symbolID, const int originalID) const {
    // The token "NATIVE-TYPE" is a wildcard for all tokens in the numeric range
    // from VOID to USAMPLER2DMSARRAY.
    // Any production which expects that terminal token, has to except any of those token
    // within this range.
    //
    if ((symbolID == GLSLTerminals::ID_NATIVE_TYPE)
        && (originalID >= GLSLTerminals::ID_BOOL) && (originalID <= GLSLTerminals::ID_USAMPLER2DMSARRAY))
    {
        return true;
    }

    return false;
}

void GLSLParser::expandParseTree(const int productionID,
                                 const std::vector<Parser::ParserSymbol*>& reductionBody)
{
    switch (productionID) {
        case 0:  // [$START$] ::= [program]
        case 1:  // [program] ::= [translation-unit]
            // Nothing to do: top-most node on stack will remain an object
            // wrapping the [translation-unit]
            break;

        case 2:  // [translation-unit] ::= [external-declaration]
            if (reductionBody.size() == 1)
            {
                GLSLExternalDeclaration* const extDecl = popNode<GLSLExternalDeclaration>();
                if (extDecl)
                    pushNode(new GLSLTranslation(extDecl));
            }
            break;

        case 3:  // [translation-unit] ::= [translation-unit] [external-declaration]
            if (reductionBody.size() == 2)
            {
                GLSLExternalDeclaration* const extDecl = popNode<GLSLExternalDeclaration>();
                GLSLTranslation* const trans = popNode<GLSLTranslation>();
                if ((extDecl) && (trans))
                {
                    trans->addExternalDeclaration(extDecl);
                    pushNode(trans);
                }
            }
            break;

        case 4:  // [external-declaration] ::= [function-definition]
        case 5:  // [external-declaration] ::= [declaration]
        case 6:  // [external-declaration] ::= [enhanced-declaration]
            // No modification on parse tree required, because the
            // production bodies correspond to C++ classes which inherit
            // from the C++ class representing the production head.
            break;

        case 7:  // [function-definition] ::= [function-prototype] [compound-statement-no-new-scope]
            if (reductionBody.size() == 2) {
                GLSLCompoundStatement* const statements = popNode<GLSLCompoundStatement>();
                GLSLFunctionPrototype* const function = popNode<GLSLFunctionPrototype>();

                if ((function) && (statements))
                    pushNode(new GLSLFunctionDefinition(function, statements));
                else
                    std::cout << "GLSLParser::expandParseTree(): error in case 7!\n";
            }
            break;

        case 8:  // [statement] ::= [compound-statement]
        case 9:  // [statement] ::= [simple-statement]
            // No further modification on parse tree required, because the
            // production bodies correspond to C++ classes which inherit
            // from the C++ class representing the production head.
            break;

        case 10:  // [simple-statement] ::= [declaration-statement]
        case 11:  // [simple-statement] ::= [expression-statement]
        case 12:  // [simple-statement] ::= [selection-statement]
        case 13:  // [simple-statement] ::= [switch-statement]
        case 14:  // [simple-statement] ::= [case-label]
        case 15:  // [simple-statement] ::= [iteration-statement]
        case 16:  // [simple-statement] ::= [jump-statement]
            // No modification on parse tree required, because the
            // production bodies correspond to C++ classes which inherit
            // from the C++ class representing the production head.
            break;

        case 17:  // [compound-statement] ::= { }
        case 23:  // [compound-statement-no-new-scope] ::= { }
            // Do not ignor empty statments, because the top-most symbol for the
            // parse tree must be a GLSLCompoundStatement.
            //
            if (reductionBody.size() == 2)
                pushNode(new GLSLCompoundStatement(0, (productionID == 17)));
            break;

        case 18:  // [compound-statement] ::= { [statement-list] }
        case 24:  // [compound-statement-no-new-scope] ::= { [statement-list] }
            if (reductionBody.size() == 3) {
                GLSLStatementList* const statements = popNode<GLSLStatementList>();
                if (statements)
                    pushNode(new GLSLCompoundStatement(statements, (productionID == 18)));
                else
                    std::cout << "GLSLParser::expandParseTree(): error in case 18, 24!\n";
            }
            break;

        case 19:  // [declaration-statement] ::= [declaration]
        case 20:  // [declaration-statement] ::= [enhanced-declaration]
            if (reductionBody.size() == 1) {
                GLSLDeclaration* const decl = popNode<GLSLDeclaration>();

                if (decl)
                    pushNode(new GLSLDeclarationStatement(decl));
                else
                    std::cout << "GLSLParser::expandParseTree(): error in case 19, 20!\n";
            }
            break;

        case 21:  // [statement-no-new-scope] ::= [compound-statement-no-new-scope]
        case 22:  // [statement-no-new-scope] ::= [simple-statement]
            // Nothing else to do: production bodies correspond to classes which
            // inherit form classes representing the production head.
            break;

        case 25:  // [statement-list] ::= [statement]
            if (reductionBody.size() == 1) {
                GLSLStatement* const statement = popNode<GLSLStatement>();

                if (statement)
                    pushNode(new GLSLStatementList(statement));
                else
                    std::cout << "GLSLParser::expandParseTree(): error in case 25!\n";
            }
            break;

        case 26:  // [statement-list] ::= [statement-list] [statement]
            if (reductionBody.size() == 2) {
                GLSLStatement* const statement = popNode<GLSLStatement>();
                GLSLStatementList* const list = popNode<GLSLStatementList>();

                if ((list) && (statement)) {
                    list->addStatement(statement);
                    pushNode(list);
                }
            }
            break;

        case 27:  // [expression-statement] ::= ;
            if (reductionBody.size() == 1)
                pushNode(new GLSLExpressionStatement(0));
            break;

        case 28:  // [expression-statement] ::= [expression] ;
            if (reductionBody.size() == 2) {
                GLSLExpression* const expr = popNode<GLSLExpression>();

                if (expr)
                    pushNode(new GLSLExpressionStatement(expr));
                else
                    std::cout << "GLSLParser::expandParseTree(): error in case 28!\n";
            }
            break;

        case 29:  // [selection-statement] ::= if ( [expression] ) [statement]
            if (reductionBody.size() == 5) {
                GLSLStatement* const stmt = popNode<GLSLStatement>();
                GLSLExpression* const cond = popNode<GLSLExpression>();

                if ((cond) && (stmt))
                    pushNode(new GLSLSelectionStatement(cond, stmt));
            }
            break;

        case 30:  // [selection-statement] ::= if ( [expression] ) [statement] else [statement]
            if (reductionBody.size() == 7) {
                GLSLStatement* const falseStmt = popNode<GLSLStatement>();
                GLSLStatement* const trueStmt = popNode<GLSLStatement>();
                GLSLExpression* const cond = popNode<GLSLExpression>();

                if ((cond) && (trueStmt) && (falseStmt))
                    pushNode(new GLSLSelectionStatement(cond, trueStmt, falseStmt));
            }
            break;

        case 31:  // [condition] ::= [expression]
            if (reductionBody.size() == 1) {
                GLSLExpression* const expr = popNode<GLSLExpression>();

                if (expr)
                    pushNode(new GLSLCondition(expr));
            }
            break;

        case 32:  // [condition] ::= [fully-specified-type] IDENTIFIER = [initializer]
            if (reductionBody.size() == 4) {
                GLSLExpression* const init = popNode<GLSLExpression>();
                IdentifierToken* const identifier =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();

                if ((spec) && (identifier) && (init))
                    pushNode(new GLSLCondition(spec, identifier, init));
            }
            break;

        case 33:  // [switch-statement] ::= switch ( [expression] ) { [statement-list] }
            // !!! Attention: There are special semantical rules for [statement-list] !!!
            if (reductionBody.size() == 7) {
                GLSLStatementList* const stmts = popNode<GLSLStatementList>();
                GLSLExpression* const expr = popNode<GLSLExpression>();

                if ((expr) && (stmts))
                    pushNode(new GLSLSwitchStatement(expr, stmts));
            }
            break;

        case 34:  // [switch-statement] ::= switch ( [expression] ) { }
            if (reductionBody.size() == 6) {
                GLSLExpression* const expr = popNode<GLSLExpression>();

                if (expr)
                    pushNode(new GLSLSwitchStatement(expr));
            }
            break;

        case 35:  // [case-label] ::= case [expression] :
            if (reductionBody.size() == 3) {
                GLSLExpression* const cond = popNode<GLSLExpression>();

                if (cond)
                    pushNode(new GLSLCaseLabel(cond));
            }
            break;

        case 36:  // [case-label] ::= default :
            if (reductionBody.size() == 2)
                pushNode(new GLSLCaseLabel());
            break;

        case 37:  // [iteration-statement] ::= while ( [condition] ) [statement-no-new-scope]
            if (reductionBody.size() == 5) {
                GLSLStatement* const stmt = popNode<GLSLStatement>();
                GLSLCondition* const cond = popNode<GLSLCondition>();

                if ((stmt) && (cond))
                    pushNode(new GLSLWhileStatement(cond, stmt));
            }
            break;

        case 38:  // [iteration-statement] ::= do [statement] while ( [expression] ) ;
            if (reductionBody.size() == 7) {
                GLSLExpression* const cond = popNode<GLSLExpression>();
                GLSLStatement* const stmt = popNode<GLSLStatement>();

                if ((cond) && (stmt))
                    pushNode(new GLSLDoWhileStatement(cond, stmt));
            }
            break;

        case 39:  // [iteration-statement] ::= for ( [for-init-statement] [for-rest-statement] ) [statement-no-new-scope]
            if (reductionBody.size() == 6) {
                GLSLStatement* const stmt = popNode<GLSLStatement>();
                GLSLForStatement* const fr = popNode<GLSLForStatement>();

                // NOTE: to be more precise, init may only be an instance of the
                // GLSLSimpleStatement subclasses GLSLExpressionStatement or
                // GLSLDeclarationStatement
                GLSLSimpleStatement* const init = popNode<GLSLSimpleStatement>();

                if ((init) && (fr) && (stmt))
                {
                    fr->setInitStatement(init);
                    fr->setBody(stmt);
                    pushNode(fr);
                }
            }
            break;

        case 40:  // [for-init-statement] ::= [expression-statement]
        case 41:  // [for-init-statement] ::= [declaration-statement]
            // Simply keep the statement nodes on the top of the stack and pop
            // them as GLSLSimpleStatements for both subclasses inherit from it
            break;

        case 42:  // [for-rest-statement] ::= ;
            if (reductionBody.size() == 1)
                pushNode(new GLSLForStatement());
            break;

        case 43:  // [for-rest-statement] ::= [condition] ;
            if (reductionBody.size() == 2) {
                GLSLCondition* const cond = popNode<GLSLCondition>();

                if (cond)
                    pushNode(new GLSLForStatement(cond, 0));
            }
            break;

        case 44:  // [for-rest-statement] ::= ; [expression]
            if (reductionBody.size() == 2) {
                GLSLExpression* const expr = popNode<GLSLExpression>();

                if (expr)
                    pushNode(new GLSLForStatement(0, expr));
            }
            break;

        case 45:  // [for-rest-statement] ::= [condition] ; [expression]
            if (reductionBody.size() == 3) {
                GLSLExpression* const expr = popNode<GLSLExpression>();
                GLSLCondition* const cond = popNode<GLSLCondition>();

                if ((cond) && (expr))
                    pushNode(new GLSLForStatement(cond, expr));
            }
            break;

        case 46:  // [jump-statement] ::= continue ;
        case 47:  // [jump-statement] ::= break ;
        case 48:  // [jump-statement] ::= return ;
        case 50:  // [jump-statement] ::= discard ;
            if (reductionBody.size() == 2) {
                Token* const token = reductionBody[1]->getToken();

                if (token)
                    pushNode(new GLSLJumpStatement(token));
            }
            break;

        case 49:  // [jump-statement] ::= return [expression] ;
            if (reductionBody.size() == 3) {
                GLSLExpression* const expr = popNode<GLSLExpression>();

                if (expr)
                    pushNode(new GLSLJumpStatement(expr));
            }
            break;

        case 51:  // [primary-expression] ::= IDENTIFIER
        case 52:  // [primary-expression] ::= INTCONSTANT
        case 53:  // [primary-expression] ::= UINTCONSTANT
        case 54:  // [primary-expression] ::= FLOATCONSTANT
        case 55:  // [primary-expression] ::= false
        case 56:  // [primary-expression] ::= true
            if (reductionBody.size() == 1) {
                Token* const token = reductionBody[0]->getToken();

                if (token != 0)
                    pushNode(new GLSLPrimaryExpression(*token));
            }
            break;

        case 57:  // [primary-expression] ::= ( [expression] )
            if (reductionBody.size() == 3) {
                GLSLExpression* const expr = popNode<GLSLExpression>();

                if (expr != 0)
                    pushNode(new GLSLPrimaryExpression(expr));
            }
            break;

        case 58:  // [postfix-expression] ::= [primary-expression]
            if (reductionBody.size() == 1) {
                GLSLPrimaryExpression* const expr = popNode<GLSLPrimaryExpression>();

                if (expr != 0)
                    pushNode(new GLSLPostfixExpression(expr));
            }
            break;

        case 59:  // [postfix-expression] ::= [postfix-expression] [ [integer-expression] ]
            if (reductionBody.size() == 4) {
                GLSLExpression* const intExpr = popNode<GLSLExpression>();
                GLSLPostfixExpression* const arr = popNode<GLSLPostfixExpression>();

                if ((arr != 0) && (intExpr != 0))
                    pushNode(new GLSLPostfixExpression(arr, intExpr));
            }
            break;

        case 60:  // [postfix-expression] ::= [function-call-or-method]
            if (reductionBody.size() == 1) {
                GLSLExpression* const funCall = popNode<GLSLFunctionCall>();

                if (funCall != 0)
                    pushNode(new GLSLPostfixExpression(funCall));
            }
            break;

        case 61:  // [postfix-expression] ::= [postfix-expression] . FIELD-SELECTION
            if (reductionBody.size() == 3) {
                GLSLPostfixExpression* const operand = popNode<GLSLPostfixExpression>();
                IdentifierToken* const field =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if ((operand != 0) && (field != 0))
                    pushNode(new GLSLPostfixExpression(operand, field));
            }
            break;

        case 62:  // [postfix-expression] ::= [postfix-expression] ++
        case 63:  // [postfix-expression] ::= [postfix-expression] --
            if (reductionBody.size() == 2) {
                GLSLPostfixExpression* const expr = popNode<GLSLPostfixExpression>();
                Token* const token = reductionBody[0]->getToken();

                if ((token != 0) && (expr != 0))
                    pushNode(new GLSLPostfixExpression(token, expr));
            }
            break;

        case 64:  // [integer-expression] ::= [expression]
            break;

        case 65:  // [function-call-or-method] ::= [function-call-generic]
            break;

        case 66:  // [function-call-or-method] ::= [postfix-expression] . [function-call-generic]
            if (reductionBody.size() == 3) {
                GLSLFunctionCall* const func = popNode<GLSLFunctionCall>();
                GLSLPostfixExpression* const postfix = popNode<GLSLPostfixExpression>();

                if ((func != 0) && (postfix != 0)) {
                    func->setPostfixExpression(postfix);
                    pushNode(func);
                }
            }
            break;

        case 67:  // [function-call-generic] ::= [function-call-header] void )
        case 68:  // [function-call-generic] ::= [function-call-header] )
        case 69:  // [function-call-generic] ::= [function-call-header-with-parameters] )
            break;

        case 70:  // [function-call-header-with-parameters] ::= [function-call-header] [assignment-expression]
        case 71:  // [function-call-header-with-parameters] ::= [function-call-header-with-parameters] , [assignment-expression]
            {
                GLSLExpression* const param = popNode<GLSLExpression>();
                GLSLFunctionCall* const func = popNode<GLSLFunctionCall>();

                if ((func != 0) && (param != 0)) {
                    func->addParameter(param);
                    pushNode(func);
                }
            }
            break;

        case 72:  // [function-call-header] ::= [type-specifier] (
            if (reductionBody.size() == 2) {
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();

                if (typeSpec != 0)
                    pushNode(new GLSLFunctionCall(typeSpec));
            }
            break;

        case 73:  // [function-call-header] ::= IDENTIFIER (
        case 74:  // [function-call-header] ::= FIELD-SELECTION (
            if (reductionBody.size() == 2) {
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[1]->getToken());

                if (name != 0)
                    pushNode(new GLSLFunctionCall(name));
            }
            break;

        case 75:  // [unary-expression] ::= [postfix-expression]
            if (reductionBody.size() == 1) {
                GLSLPostfixExpression* const expr = popNode<GLSLPostfixExpression>();

                if (expr != 0)
                    pushNode(new GLSLUnaryExpression(expr));
            }
            break;

        case 76:  // [unary-expression] ::= ++ [unary-expression]
        case 77:  // [unary-expression] ::= -- [unary-expression]
        case 78:  // [unary-expression] ::= + [unary-expression]
        case 79:  // [unary-expression] ::= - [unary-expression]
        case 80:  // [unary-expression] ::= ! [unary-expression]
        case 81:  // [unary-expression] ::= ~ [unary-expression]
            if (reductionBody.size() == 2) {
                GLSLUnaryExpression* const expr = popNode<GLSLUnaryExpression>();
                Token* const token = reductionBody[1]->getToken();

                if ((expr != 0) && (token != 0))
                    pushNode(new GLSLUnaryExpression(token, expr));
            }
            break;

        case 82:  // [multiplicative-expression] ::= [unary-expression]
        case 86:  // [additive-expression] ::= [multiplicative-expression]
        case 89:  // [shift-expression] ::= [additive-expression]
        case 92:  // [relational-expression] ::= [shift-expression]
        case 97:  // [equality-expression] ::= [relational-expression]
        case 100:  // [and-expression] ::= [equality-expression]
        case 102:  // [exclusive-or-expression] ::= [and-expression]
        case 104:  // [inclusive-or-expression] ::= [exclusive-or-expression]
        case 106:  // [logical-and-expression] ::= [inclusive-or-expression]
        case 108:  // [logical-xor-expression] ::= [logical-and-expression]
        case 110:  // [logical-or-expression] ::= [logical-xor-expression]
            break;

        case 83:  // [multiplicative-expression] ::= [multiplicative-expression] * [unary-expression]
        case 84:  // [multiplicative-expression] ::= [multiplicative-expression] / [unary-expression]
        case 85:  // [multiplicative-expression] ::= [multiplicative-expression] % [unary-expression]
        case 87:  // [additive-expression] ::= [additive-expression] + [multiplicative-expression]
        case 88:  // [additive-expression] ::= [additive-expression] - [multiplicative-expression]
        case 90:  // [shift-expression] ::= [shift-expression] << [additive-expression]
        case 91:  // [shift-expression] ::= [shift-expression] >> [additive-expression]
        case 93:  // [relational-expression] ::= [relational-expression] < [shift-expression]
        case 94:  // [relational-expression] ::= [relational-expression] > [shift-expression]
        case 95:  // [relational-expression] ::= [relational-expression] <= [shift-expression]
        case 96:  // [relational-expression] ::= [relational-expression] >= [shift-expression]
        case 98:  // [equality-expression] ::= [equality-expression] == [relational-expression]
        case 99:  // [equality-expression] ::= [equality-expression] != [relational-expression]
        case 101:  // [and-expression] ::= [and-expression] & [equality-expression]
        case 103:  // [exclusive-or-expression] ::= [exclusive-or-expression] ^ [and-expression]
        case 105:  // [inclusive-or-expression] ::= [inclusive-or-expression] | [exclusive-or-expression]
        case 107:  // [logical-and-expression] ::= [logical-and-expression] && [inclusive-or-expression]
        case 109:  // [logical-xor-expression] ::= [logical-xor-expression] ^^ [logical-and-expression]
        case 111:  // [logical-or-expression] ::= [logical-or-expression] || [logical-xor-expression]
            if (reductionBody.size() == 3) {
                GLSLExpression* const rhs = popNode<GLSLExpression>();
                GLSLExpression* const lhs = popNode<GLSLExpression>();
                Token* const token = reductionBody[1]->getToken();

                if ((token != 0) && (lhs != 0) && (rhs != 0))
                    pushNode(new GLSLBinaryExpression(token, lhs, rhs));
            }
            break;

        case 112:  // [conditional-expression] ::= [logical-or-expression]
            break;

        case 113:  // [conditional-expression] ::= [logical-or-expression] ? [expression] : [assignment-expression]
            if (reductionBody.size() == 5) {
                GLSLExpression* const falseExpr = popNode<GLSLExpression>();
                GLSLExpression* const trueExpr = popNode<GLSLExpression>();
                GLSLExpression* const condition = popNode<GLSLExpression>();

                if ((condition != 0) && (trueExpr != 0) && (falseExpr != 0))
                    pushNode(new GLSLConditionalExpression(condition, trueExpr, falseExpr));
            }
            break;

        case 114:  // [assignment-expression] ::= [conditional-expression]
            break;

        case 115:  // [assignment-expression] ::= [unary-expression] = [assignment-expression]
        case 116:  // [assignment-expression] ::= [unary-expression] *= [assignment-expression]
        case 117:  // [assignment-expression] ::= [unary-expression] /= [assignment-expression]
        case 118:  // [assignment-expression] ::= [unary-expression] %= [assignment-expression]
        case 119:  // [assignment-expression] ::= [unary-expression] += [assignment-expression]
        case 120:  // [assignment-expression] ::= [unary-expression] -= [assignment-expression]
        case 121:  // [assignment-expression] ::= [unary-expression] <<= [assignment-expression]
        case 122:  // [assignment-expression] ::= [unary-expression] >>= [assignment-expression]
        case 123:  // [assignment-expression] ::= [unary-expression] &= [assignment-expression]
        case 124:  // [assignment-expression] ::= [unary-expression] ^= [assignment-expression]
        case 125:  // [assignment-expression] ::= [unary-expression] |= [assignment-expression]
            if (reductionBody.size() == 3) {
                GLSLExpression* const assign = popNode<GLSLExpression>();
                GLSLUnaryExpression* const unary = popNode<GLSLUnaryExpression>();
                Token* const op = reductionBody[1]->getToken();

                if ((op != 0) && (unary != 0) && (assign != 0))
                    pushNode(new GLSLAssignmentExpression(op, unary, assign));
            }
            break;

        case 126:  // [expression] ::= [assignment-expression]
            {
                GLSLExpression* const assign = popNode<GLSLExpression>();

                if (assign != 0)
                    pushNode(new GLSLExpressionList(assign));
            }
            break;

        case 127:  // [expression] ::= [expression] , [assignment-expression]
            {
                GLSLExpression* const assign = popNode<GLSLExpression>();
                GLSLExpressionList* const list = popNode<GLSLExpressionList>();

                if ((list != 0) && (assign != 0)) {
                    list->addExpression(assign);
                    pushNode(list);
                }
            }
            break;

        case 128:  // [constant-expression] ::= [conditional-expression]
            break;

        case 129:  // [enhanced-declaration] ::= LEADING-ANNOTATION [declaration]
            if (reductionBody.size() == 2) {
                GLSLDeclaration* const decl = popNode<GLSLDeclaration>();
                AnnotationToken* const annotation =
                    dynamic_cast<AnnotationToken* const>(reductionBody[1]->getToken());

                if ((decl != 0) && (annotation != 0)) {
                    decl->setLeadingAnnotation(annotation);
                    pushNode(decl);
                }
            }
            break;

        case 130:  // [enhanced-declaration] ::= [declaration] TRAILING-ANNOTATION
            if (reductionBody.size() == 2) {
                GLSLDeclaration* const decl = popNode<GLSLDeclaration>();
                AnnotationToken* const annotation =
                    dynamic_cast<AnnotationToken* const>(reductionBody[0]->getToken());

                if ((decl != 0) && (annotation != 0)) {
                    decl->setTrailingAnnotation(annotation);
                    pushNode(decl);
                }
            }
            break;

        case 131:  // [enhanced-declaration] ::= LEADING-ANNOTATION [declaration] TRAILING-ANNOTATION
            if (reductionBody.size() == 3) {
                GLSLDeclaration* const decl = popNode<GLSLDeclaration>();
                AnnotationToken* const trailing =
                    dynamic_cast<AnnotationToken* const>(reductionBody[0]->getToken());
                AnnotationToken* const leading =
                    dynamic_cast<AnnotationToken* const>(reductionBody[2]->getToken());

                if ((decl != 0) && (leading != 0) && (trailing)) {
                    decl->setLeadingAnnotation(leading);
                    decl->setTrailingAnnotation(trailing);
                    pushNode(decl);
                }
            }
            break;

        case 132:  // [declaration] ::= [function-prototype] ;
            if (reductionBody.size() == 2) {
                GLSLFunctionPrototype* const funcProto = popNode<GLSLFunctionPrototype>();
                if (funcProto)
                    pushNode(new GLSLFunctionDeclaration(funcProto));
            }
            break;

        case 133:  // [declaration] ::= [init-declarator-list] ;
            break;

        case 134:  // [declaration] ::= precision [precision-qualifier] [type-specifier-no-prec] ;
            if (reductionBody.size() == 4) {
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                GLSLPrecisionQualifier* const prec = popNode<GLSLPrecisionQualifier>();

                if ((spec != 0) && (prec != 0)) {
                    spec->setPrecision(prec);
                    IdentifierToken id(GLSLTerminals::ID_IDENTIFIER, "anonym");
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(&id, false, 0)));
                }
            }
            break;

        case 135:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } ;
            if (reductionBody.size() == 6) {
                GLSLStructDeclaratorList* const structDecls = popNode<GLSLStructDeclaratorList>();
                GLSLTypeQualifier* const typeQuali = popNode<GLSLTypeQualifier>();
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[4]->getToken());

                if ((name != 0) && (typeQuali != 0) && (structDecls != 0))
                    pushNode(new GLSLFieldDeclaration(typeQuali, name, structDecls));
            }
            break;

        case 136:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } IDENTIFIER ;
            if (reductionBody.size() == 7) {
                GLSLStructDeclaratorList* const structDecls = popNode<GLSLStructDeclaratorList>();
                GLSLTypeQualifier* const typeQuali = popNode<GLSLTypeQualifier>();
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[5]->getToken());
                IdentifierToken* const varName =
                    dynamic_cast<IdentifierToken* const>(reductionBody[1]->getToken());

                if ((name != 0) && (varName != 0) && (typeQuali != 0) && (structDecls != 0)) {
                    GLSLFieldDeclaration* const fieldDecl =
                        new GLSLFieldDeclaration(typeQuali, name, structDecls);
                    fieldDecl->setVariable(new GLSLVariable(varName, false, 0));
                    pushNode(fieldDecl);
                }
            }
            break;

        case 137:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } IDENTIFIER [ ] ;
            if (reductionBody.size() == 9) {
                GLSLStructDeclaratorList* const structDecls = popNode<GLSLStructDeclaratorList>();
                GLSLTypeQualifier* const typeQuali = popNode<GLSLTypeQualifier>();
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[7]->getToken());
                IdentifierToken* const varName =
                    dynamic_cast<IdentifierToken* const>(reductionBody[3]->getToken());

                if ((name != 0) && (varName != 0) && (typeQuali != 0) && (structDecls != 0)) {
                    GLSLFieldDeclaration* const fieldDecl =
                        new GLSLFieldDeclaration(typeQuali, name, structDecls);
                    fieldDecl->setVariable(new GLSLVariable(varName, true, 0));
                    pushNode(fieldDecl);
                }
            }
            break;

        case 138:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } IDENTIFIER [ [constant-expression] ] ;
            if (reductionBody.size() == 10) {
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                GLSLStructDeclaratorList* const structDecls = popNode<GLSLStructDeclaratorList>();
                GLSLTypeQualifier* const typeQuali = popNode<GLSLTypeQualifier>();
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[8]->getToken());
                IdentifierToken* const varName =
                    dynamic_cast<IdentifierToken* const>(reductionBody[4]->getToken());

                if ((name != 0) && (varName != 0) && (typeQuali != 0) && (structDecls != 0)
                    && (constExpr != 0))
                {
                    GLSLFieldDeclaration* const fieldDecl =
                        new GLSLFieldDeclaration(typeQuali, name, structDecls);
                    fieldDecl->setVariable(new GLSLVariable(varName, false, constExpr));
                    pushNode(fieldDecl);
                }
            }
            break;

        case 139:  // [declaration] ::= [type-qualifier] ;
            if (reductionBody.size() == 2) {
                GLSLTypeQualifier* const typeQuali = popNode<GLSLTypeQualifier>();

                if (typeQuali != 0) {
                    IdentifierToken id(GLSLTerminals::ID_IDENTIFIER, "anonym");
                    GLSLTypeSpecifier* const spec = new GLSLTypeSpecifier(0, 0, false, 0);
                    spec->setTypeQualifier(typeQuali);
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(&id,false, 0)));
                }
            }
            break;

        case 140:  // [function-prototype] ::= [function-header] )
        case 141:  // [function-prototype] ::= [function-header-with-parameters] )
            // Top node is a GLSLFunctionPrototype node and no new node for a
            // semantic action is required as closing brackets are implicitly
            // assumed to be present.
            {
            GLSLFunctionPrototype* const funcProto = popNode<GLSLFunctionPrototype>();
            pushNode(funcProto);
            }
            break;

        case 142:  // [function-header-with-parameters] ::= [function-header] [parameter-declaration]
        case 143:  // [function-header-with-parameters] ::= [function-header-with-parameters] , [parameter-declaration]
            {
                GLSLParameter* const param = popNode<GLSLParameter>();
                GLSLFunctionPrototype* const funcProto = popNode<GLSLFunctionPrototype>();

                if ((funcProto != 0) && (param != 0)) {
                    funcProto->addParameter(param);
                    pushNode(funcProto);
                }
            }
            break;

        case 144:  // [function-header] ::= [fully-specified-type] IDENTIFIER (
            if (reductionBody.size() == 3) {
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[1]->getToken());

                if ((typeSpec != 0) && (name != 0))
                    pushNode(new GLSLFunctionPrototype(name, typeSpec));
            }
            break;

        case 145:  // [parameter-declarator] ::= [type-specifier] IDENTIFIER
            if (reductionBody.size() == 2) {
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if ((name != 0) && (typeSpec != 0)) {
                    GLSLVariable* const var = new GLSLVariable(name, true, 0);
                    pushNode(new GLSLParameter(typeSpec, 0, false, var));
                }
            }
            break;

        case 146:  // [parameter-declarator] ::= [type-specifier] IDENTIFIER [ [constant-expression] ]
            if (reductionBody.size() == 5) {
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[5]->getToken());

                if ((name != 0) && (typeSpec != 0) && (constExpr != 0)) {
                    GLSLVariable* const var = new GLSLVariable(name, true, constExpr);
                    pushNode(new GLSLParameter(typeSpec, 0, false, var));
                }
            }
            break;

        case 147:  // [parameter-declaration] ::= const [parameter-declarator]
            if (reductionBody.size() == 3) {
                GLSLParameter* const param = popNode<GLSLParameter>();

                if (param != 0) {
                    param->setIsConst(true);
                    pushNode(param);
                }
            }
            break;

        case 148:  // [parameter-declaration] ::= const [parameter-qualifier] [parameter-declarator]
            if (reductionBody.size() == 3) {
                GLSLParameter* const param = popNode<GLSLParameter>();
                GLSLParameterQualifier* const paramQuali = popNode<GLSLParameterQualifier>();

                if ((param != 0) && (paramQuali != 0)) {
                    param->setParameterQulifier(paramQuali);
                    param->setIsConst(true);
                    pushNode(param);
                }
            }
            break;

        case 149:  // [parameter-declaration] ::= [parameter-declarator]
            break;

        case 150:  // [parameter-declaration] ::= [parameter-qualifier] [parameter-declarator]
            if (reductionBody.size() == 2) {
                GLSLParameter* const param = popNode<GLSLParameter>();
                GLSLParameterQualifier* const paramQuali = popNode<GLSLParameterQualifier>();

                if ((param != 0) && (paramQuali != 0)) {
                    param->setParameterQulifier(paramQuali);
                    pushNode(param);
                }
            }
            break;

        case 151:  // [parameter-declaration] ::= const [parameter-type-specifier]
            if (reductionBody.size() == 2) {
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();

                if (typeSpec != 0)
                    pushNode(new GLSLParameter(typeSpec, 0, true, 0));
            }
            break;

        case 152:  // [parameter-declaration] ::= const [parameter-qualifier] [parameter-type-specifier]
            if (reductionBody.size() == 3) {
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();
                GLSLParameterQualifier* const paramQuali = popNode<GLSLParameterQualifier>();

                if ((typeSpec != 0) && (paramQuali != 0))
                    pushNode(new GLSLParameter(typeSpec, paramQuali, true, 0));
            }
            break;

        case 153:  // [parameter-declaration] ::= [parameter-type-specifier]
            if (reductionBody.size() == 1) {
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();

                if (typeSpec != 0)
                    pushNode(new GLSLParameter(typeSpec, 0, false, 0));
            }
            break;

        case 154:  // [parameter-declaration] ::= [parameter-qualifier] [parameter-type-specifier]
            if (reductionBody.size() == 2) {
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();
                GLSLParameterQualifier* const paramQuali = popNode<GLSLParameterQualifier>();

                if ((typeSpec != 0) && (paramQuali != 0))
                    pushNode(new GLSLParameter(typeSpec, paramQuali, false, 0));
            }
            break;

        case 155:  // [parameter-qualifier] ::= in
        case 156:  // [parameter-qualifier] ::= out
        case 157:  // [parameter-qualifier] ::= inout
            if (reductionBody.size() == 1) {
                Token* const paramQuali = reductionBody[0]->getToken();

                if (paramQuali != 0)
                    pushNode(new GLSLParameterQualifier(paramQuali));
            }
            break;

        case 158:  // [parameter-type-specifier] ::= [type-specifier]
            break;

        case 159:  // [init-declarator-list] ::= [single-declaration]
            break;

        case 160:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER
            if (reductionBody.size() == 3) {
                GLSLDeclarationList* const decl = popNode<GLSLDeclarationList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if ((decl != 0) && (id != 0)) {
                    decl->addVariable(new GLSLVariable(id, false, 0));
                    pushNode(decl);
                }
            }
            break;

        case 161:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ ]
            if (reductionBody.size() == 5) {
                GLSLDeclarationList* const decl = popNode<GLSLDeclarationList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());

                if ((decl != 0) && (id != 0)) {
                    decl->addVariable(new GLSLVariable(id, true, 0));
                    pushNode(decl);
                }
            }
            break;

        case 162:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ [constant-expression] ]
            if (reductionBody.size() == 6) {
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                GLSLDeclarationList* const decl = popNode<GLSLDeclarationList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[3]->getToken());

                if ((decl != 0) && (id != 0) && (constExpr != 0)) {
                    decl->addVariable(new GLSLVariable(id, true, constExpr, 0));
                    pushNode(decl);
                }
            }
            break;

        case 163:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ ] = [initializer]
            if (reductionBody.size() == 7) {
                GLSLExpression* const init = popNode<GLSLExpression>();
                GLSLDeclarationList* const decl = popNode<GLSLDeclarationList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[4]->getToken());

                if ((decl != 0) && (id != 0) && (init != 0)) {
                    decl->addVariable(new GLSLVariable(id, true, 0, init));
                    pushNode(decl);
                }
            }
            break;

        case 164:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ [constant-expression] ] = [initializer]
            if (reductionBody.size() == 8) {
                GLSLExpression* const init = popNode<GLSLExpression>();
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                GLSLDeclarationList* const decl = popNode<GLSLDeclarationList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[4]->getToken());

                if ((decl != 0) && (id != 0) && (constExpr != 0) && (init != 0)) {
                    decl->addVariable(new GLSLVariable(id, true, constExpr, init));
                    pushNode(decl);
                }
            }
            break;

        case 165:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER = [initializer]
            if (reductionBody.size() == 5) {
                GLSLExpression* const init = popNode<GLSLExpression>();
                GLSLDeclarationList* const decl = popNode<GLSLDeclarationList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());

                if ((decl != 0) && (id != 0) && (init != 0)) {
                    decl->addVariable(new GLSLVariable(id, false, 0, init));
                    pushNode(decl);
                }
            }
            break;

        case 166:  // [single-declaration] ::= [fully-specified-type]
            if (reductionBody.size() == 1) {
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();

                // The specifier can either be a struct specifier or a "normal"
                // type specifier
                //
                if (spec != 0) {
                    GLSLStructSpecifier* const structSpec =
                        dynamic_cast<GLSLStructSpecifier* const>(spec);

                    if (structSpec == 0) {
                        IdentifierToken anonym(GLSLTerminals::ID_IDENTIFIER, "anonym");
                        pushNode(new GLSLDeclarationList(spec, new GLSLVariable(&(anonym), false, 0)));
                    } else {
                        GLSLVariable* const var = new GLSLVariable(structSpec->getStructName(), false, 0);
                        pushNode(new GLSLDeclarationList(structSpec, var));
                    }
                }
            }
            return;

        case 167:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER
            if (reductionBody.size() == 2) {
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if ((spec != 0) && (id != 0))
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(id, false, 0)));
            }
            return;

        case 168:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ ]
            if (reductionBody.size() == 4) {
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());

                if ((spec != 0) && (id != 0))
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(id, true, 0)));
            }
            break;

        case 169:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ [constant-expression] ]
            if (reductionBody.size() == 5) {
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[3]->getToken());

                if ((spec != 0) && (id != 0) && (constExpr != 0))
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(id, true, constExpr)));
            }
            break;

        case 170:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ ] = [initializer]
            if (reductionBody.size() == 6) {
                GLSLExpression* const init = popNode<GLSLExpression>();
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[4]->getToken());

                if ((spec != 0) && (id != 0) && (init != 0))
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(id, true, 0, init)));
            }
            break;

        case 171:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ [constant-expression] ] = [initializer]
            if (reductionBody.size() == 7) {
                GLSLExpression* const init = popNode<GLSLExpression>();
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[5]->getToken());

                if ((spec != 0) && (id != 0) && (constExpr != 0) && (init != 0))
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(id, true, constExpr, init)));
            }
            break;

        case 172:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER = [initializer]
            if (reductionBody.size() == 4) {
                GLSLExpression* const init = popNode<GLSLExpression>();
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());

                if ((spec != 0) && (id != 0) && (init != 0))
                    pushNode(new GLSLDeclarationList(spec, new GLSLVariable(id, false, 0, init)));
            }
            break;


        case 173:  // [fully-specified-type] ::= [type-specifier]
            break;

        case 174:  // [fully-specified-type] ::= [type-qualifier] [type-specifier]
            if (reductionBody.size() == 2) {
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                GLSLTypeQualifier* const qual = popNode<GLSLTypeQualifier>();

                if ((spec != 0) && (qual != 0)) {
                    spec->setTypeQualifier(qual);
                    pushNode(spec);
                }
            }
            break;

        case 175:  // [interpolation-qualifier] ::= smooth
        case 176:  // [interpolation-qualifier] ::= flat
        case 177:  // [interpolation-qualifier] ::= noperspective
            if (reductionBody.size() == 1) {
                Token* const token = reductionBody[0]->getToken();
                if (token != 0)
                    pushNode(new GLSLInterpolationQualifier(token));
            }
            break;

        case 178:  // [layout-qualifier] ::= layout ( [layout-qualifier-id-list] )
            break;

        case 179:  // [layout-qualifier-id-list] ::= IDENTIFIER
            if (reductionBody.size() == 1) {
                IdentifierToken* const identifier =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (identifier != 0)
                    pushNode(new GLSLLayoutQualifier(identifier, 0));
            }
            break;

        case 180:  // [layout-qualifier-id-list] ::= IDENTIFIER = INTCONSTANT
            if (reductionBody.size() == 3) {
                ConstantToken* const constToken =
                    dynamic_cast<ConstantToken* const>(reductionBody[0]->getToken());
                IdentifierToken* const identifier =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());

                if ((identifier != 0) && (constToken != 0))
                    pushNode(new GLSLLayoutQualifier(identifier, constToken));
            }
            break;

        case 181:  // [layout-qualifier-id-list] ::= [layout-qualifier-id-list] , IDENTIFIER
            if (reductionBody.size() == 3) {
                IdentifierToken* const identifier =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());
                GLSLLayoutQualifier* const layout = popNode<GLSLLayoutQualifier>();

                if ((layout != 0) && (identifier != 0)) {
                    layout->addElement(identifier, 0);
                    pushNode(layout);
                }
            }
            break;

        case 182:  // [layout-qualifier-id-list] ::= [layout-qualifier-id-list] , IDENTIFIER = INTCONSTANT
            if (reductionBody.size() == 5) {
                ConstantToken* const constToken =
                    dynamic_cast<ConstantToken* const>(reductionBody[0]->getToken());
                IdentifierToken* const identifier =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());
                GLSLLayoutQualifier* const layout = popNode<GLSLLayoutQualifier>();

                if ((layout != 0) && (identifier != 0) && (constToken != 0)) {
                    layout->addElement(identifier, constToken);
                    pushNode(layout);
                }
            }
            break;

        case 183:  // [type-qualifier] ::= [storage-qualifier]
            if (reductionBody.size() == 1) {
                GLSLStorageQualifier* const stor = popNode<GLSLStorageQualifier>();
                if (stor != 0)
                    pushNode(new GLSLTypeQualifier(stor, 0, 0, false));
            }
            break;

        case 184:  // [type-qualifier] ::= [layout-qualifier]
            {
                GLSLLayoutQualifier* const layout = popNode<GLSLLayoutQualifier>();

                if (layout != 0)
                    pushNode(new GLSLTypeQualifier(0, 0, layout, false));
            }
            break;

        case 185:  // [type-qualifier] ::= [layout-qualifier] [storage-qualifier]
            {
                GLSLStorageQualifier* const storage = popNode<GLSLStorageQualifier>();
                GLSLLayoutQualifier* const layout = popNode<GLSLLayoutQualifier>();

                if ((storage != 0) && (layout != 0))
                    pushNode(new GLSLTypeQualifier(storage, 0, layout, false));
            }
            break;

        case 186:  // [type-qualifier] ::= [interpolation-qualifier] [storage-qualifier]
            if (reductionBody.size() == 2) {
                GLSLStorageQualifier* const stor = popNode<GLSLStorageQualifier>();
                GLSLInterpolationQualifier* const inter = popNode<GLSLInterpolationQualifier>();
                if ((stor != 0) && (inter != 0))
                    pushNode(new GLSLTypeQualifier(stor, inter, 0, false));
            }
            break;

        case 187:  // [type-qualifier] ::= [interpolation-qualifier]
            if (reductionBody.size() == 1) {
                GLSLInterpolationQualifier* const inter = popNode<GLSLInterpolationQualifier>();
                if (inter != 0)
                    pushNode(new GLSLTypeQualifier(0, inter, 0, false));
            }
            break;

        case 188:  // [type-qualifier] ::= invariant [storage-qualifier]
            if (reductionBody.size() == 2) {
                GLSLStorageQualifier* const stor = popNode<GLSLStorageQualifier>();
                if (stor != 0)
                    pushNode(new GLSLTypeQualifier(stor, 0, 0, true));
            }
            break;

        case 189:  // [type-qualifier] ::= invariant [interpolation-qualifier] [storage-qualifier]
            if (reductionBody.size() == 3) {
                GLSLStorageQualifier* const stor = popNode<GLSLStorageQualifier>();
                GLSLInterpolationQualifier* const inter = popNode<GLSLInterpolationQualifier>();
                if ((stor != 0) && (inter != 0))
                    pushNode(new GLSLTypeQualifier(stor, inter, 0, true));
            }
            break;

        case 190:  // [type-qualifier] ::= invariant
            if (reductionBody.size() == 1)
                pushNode(new GLSLTypeQualifier(0, 0, 0, true));
            break;


        case 191:  // [storage-qualifier] ::= const
        case 192:  // [storage-qualifier] ::= attribute
        case 193:  // [storage-qualifier] ::= varying
        case 194:  // [storage-qualifier] ::= centroid varying
        case 195:  // [storage-qualifier] ::= in
        case 196:  // [storage-qualifier] ::= out
        case 197:  // [storage-qualifier] ::= centroid in
        case 198:  // [storage-qualifier] ::= centroid out
        case 199:  // [storage-qualifier] ::= uniform
            if ((reductionBody.size() >= 1) && (reductionBody.size() <= 2)) {
                Token* const token = reductionBody[0]->getToken();
                Token* const centroidToken = ((reductionBody.size() == 2) ? reductionBody[1]->getToken() : 0);

                if (token != 0) {
                    bool centroid = ((centroidToken != 0) ?
                        (centroidToken->getTokenID() == GLSLTerminals::ID_CENTROID) : false);
                    pushNode(new GLSLStorageQualifier(token, centroid));
                }
            }
            break;

        case 202:  // [precision-qualifier] ::= highp
        case 203:  // [precision-qualifier] ::= mediump
        case 204:  // [precision-qualifier] ::= lowp
            if (reductionBody.size() == 1) {
                Token* const token = reductionBody[0]->getToken();
                if (token != 0)
                    pushNode(new GLSLPrecisionQualifier(token));
            }
            break;

        case 200:  // [type-specifier] ::= [type-specifier-no-prec]
        case 205:  // [type-specifier-no-prec] ::= [type-specifier-nonarray]
            break;

        case 201:  // [type-specifier] ::= [precision-qualifier] [type-specifier-no-prec]
            if (reductionBody.size() == 2) {
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();
                GLSLPrecisionQualifier* const prec = popNode<GLSLPrecisionQualifier>();

                if ((spec != 0) && (prec != 0)) {
                    spec->setPrecision(prec);
                    pushNode(spec);
                }
            }
            break;

        case 206:  // [type-specifier-no-prec] ::= [type-specifier-nonarray] [ ]
            if (reductionBody.size() == 3) {
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();

                if (spec != 0) {
                    spec->setIsArray(true);
                    pushNode(spec);
                }
            }
            break;

        case 207:  // [type-specifier-no-prec] ::= [type-specifier-nonarray] [ [constant-expression] ]
            if (reductionBody.size() == 4) {
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                GLSLTypeSpecifier* const spec = popNode<GLSLTypeSpecifier>();

                if ((spec != 0) && (constExpr)) {
                    spec->setIsArray(true);
                    spec->setNumArrayElements(constExpr);
                    pushNode(spec);
                }
            }
            break;

        case 208:  // [type-specifier-nonarray] ::= void
        case 209:  // [type-specifier-nonarray] ::= NATIVE-TYPE
        case 211:  // [type-specifier-nonarray] ::= TYPE-NAME
            if (reductionBody.size() == 1) {
                Token* const token = reductionBody[0]->getToken();
                if (token != 0)
                    pushNode(new GLSLTypeSpecifier(*token, 0, false, 0));
            }
            break;

        case 210:  // [type-specifier-nonarray] ::= [struct-specifier]
            break;

        case 212:  // [struct-specifier] ::= struct IDENTIFIER { [struct-declaration-list] }
            if (reductionBody.size() == 5) {
                GLSLStructDeclaratorList* const decls = popNode<GLSLStructDeclaratorList>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[3]->getToken());

                if ((decls != 0) && (id != 0))
                    pushNode(new GLSLStructSpecifier(id, decls));
            }
            break;

        case 213:  // [struct-specifier] ::= struct { [struct-declaration-list] }
            if (reductionBody.size() == 4) {
                GLSLStructDeclaratorList* const decls = popNode<GLSLStructDeclaratorList>();
                IdentifierToken id(GLSLTerminals::ID_IDENTIFIER, "anonym");

                if (decls != 0)
                    pushNode(new GLSLStructSpecifier(&id, decls));
            }
            break;

        case 214:  // [struct-declaration-list] ::= [struct-declaration]
            if (reductionBody.size() == 1) {
                GLSLDeclarationList* const decls = popNode<GLSLDeclarationList>();

                if (decls != 0)
                    pushNode(new GLSLStructDeclaratorList(decls));
            }
            break;

        case 215:  // [struct-declaration-list] ::= [struct-declaration-list] [struct-declaration]
            if (reductionBody.size() == 2) {
                GLSLDeclarationList* const decls = popNode<GLSLDeclarationList>();
                GLSLStructDeclaratorList* const list = popNode<GLSLStructDeclaratorList>();

                if ((list != 0) && (decls != 0)) {
                    list->addDeclarator(decls);
                    pushNode(list);
                }
            }
            break;

        case 216:  // [struct-declaration] ::= [type-specifier] [struct-declarator-list] ;
            if (reductionBody.size() == 3) {
                GLSLDeclarationList* const decls = popNode<GLSLDeclarationList>();
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();

                if ((decls != 0) && (typeSpec != 0)) {
                    decls->setTypeSpecifier(typeSpec);
                    pushNode(decls);
                }
            }
            break;

        case 217:  // [struct-declaration] ::= [type-qualifier] [type-specifier] [struct-declarator-list] ;
            if (reductionBody.size() == 4) {
                GLSLDeclarationList* const decls = popNode<GLSLDeclarationList>();
                GLSLTypeSpecifier* const typeSpec = popNode<GLSLTypeSpecifier>();
                GLSLTypeQualifier* const typeQuali = popNode<GLSLTypeQualifier>();

                if ((decls != 0) && (typeSpec != 0) && (typeQuali != 0)) {
                    typeSpec->setTypeQualifier(typeQuali);
                    decls->setTypeSpecifier(typeSpec);
                    pushNode(decls);
                }
            }
            break;

        case 218:  // [struct-declarator-list] ::= [struct-declarator]
            if (reductionBody.size() == 1) {
                GLSLVariable* const var = popNode<GLSLVariable>();

                if (var != 0)
                    pushNode(new GLSLDeclarationList(0, var));
            }
            break;

        case 219:  // [struct-declarator-list] ::= [struct-declarator-list] , [struct-declarator]
            if (reductionBody.size() == 3) {
                GLSLVariable* const var = popNode<GLSLVariable>();
                GLSLDeclarationList* const decls = popNode<GLSLDeclarationList>();

                if ((var != 0) && (decls != 0)) {
                    decls->addVariable(var);
                    pushNode(decls);
                }
            }
            break;

        case 220:  // [struct-declarator] ::= IDENTIFIER
            if (reductionBody.size() == 1) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (id != 0)
                    pushNode(new GLSLVariable(id, false, 0));
            }
            break;

        case 221:  // [struct-declarator] ::= IDENTIFIER [ ]
            if (reductionBody.size() == 3) {
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[2]->getToken());

                if (id != 0)
                    pushNode(new GLSLVariable(id, true, 0));
            }
            break;

        case 222:  // [struct-declarator] ::= IDENTIFIER [ [constant-expression] ]
            if (reductionBody.size() == 4) {
                GLSLConditionalExpression* const constExpr = popNode<GLSLConditionalExpression>();
                IdentifierToken* const id =
                    dynamic_cast<IdentifierToken* const>(reductionBody[3]->getToken());

                if ((id != 0) && (constExpr != 0))
                    pushNode(new GLSLVariable(id, true, constExpr));
            }
            break;

        case 223:  // [initializer] ::= [assignment-expression]
            // Only production, so leave assignment expression on top of the
            // stack and use GLSLAssignmentExpression in all cases where
            // [intializer] is exprected
            break;

    }   // switch (productionID)
}   // expandParseTree()

// The following lines are generated code, except for some few modifications which are marked.
//
// DO NOT MODIFY THEM MANUALLY!!!

ParserAction* GLSLParser::action(const int stateID, const int realSymbolID) const {
    // added manually
    //
    int symbolID = realSymbolID;
    if ((symbolID >= GLSLTerminals::ID_BOOL) && (symbolID <= GLSLTerminals::ID_USAMPLER2DMSARRAY))
        symbolID = GLSLTerminals::ID_NATIVE_TYPE;
    //
    // added manually

    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 1:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 2:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 24:  // [
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(209);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 3:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 24:  // [
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(211);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 4:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(192);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 5:
            switch (symbolID) {
                case 48:  // in
                    return new ParserActionTransition(45);

                case 50:  // out
                    return new ParserActionTransition(46);

                case 52:  // varying
                    return new ParserActionTransition(47);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 6:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(191);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 7:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(5);

                case 150:  // TRAILING-ANNOTATION
                    return new ParserActionTransition(48);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 8:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 9:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(2);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 10:
            switch (symbolID) {
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 48:  // in
                case 50:  // out
                case 51:  // uniform
                case 52:  // varying
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(176);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 11:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(166);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(49);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 12:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(4);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 13:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 23:  // )
                    return new ParserActionTransition(50);

                case 47:  // const
                    return new ParserActionTransition(51);

                case 48:  // in
                    return new ParserActionTransition(52);

                case 49:  // inout
                    return new ParserActionTransition(53);

                case 50:  // out
                    return new ParserActionTransition(54);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 14:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionTransition(60);

                case 31:  // ,
                    return new ParserActionTransition(61);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 15:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(62);

                case 26:  // {
                    return new ParserActionTransition(64);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 16:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(202);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 17:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(195);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 18:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(65);

                case 33:  // ;
                    return new ParserActionTransition(66);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 19:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(187);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 20:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(190);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 21:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(70);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 22:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(184);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 23:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(204);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 24:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(203);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 25:
            switch (symbolID) {
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 48:  // in
                case 50:  // out
                case 51:  // uniform
                case 52:  // varying
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(177);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 26:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(196);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 27:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 28:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 29:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionAccept();

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 30:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(159);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 31:
            switch (symbolID) {
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 48:  // in
                case 50:  // out
                case 51:  // uniform
                case 52:  // varying
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(175);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 32:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(183);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 33:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(74);

                case 26:  // {
                    return new ParserActionTransition(75);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 34:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 24:  // [
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(210);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 35:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(1);

                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 36:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 33:  // ;
                    return new ParserActionTransition(77);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(78);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 37:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(173);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 38:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(200);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 39:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(205);

                case 24:  // [
                    return new ParserActionTransition(80);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 40:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(199);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 41:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(193);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 42:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 24:  // [
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(208);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 43:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(129);

                case 150:  // TRAILING-ANNOTATION
                    return new ParserActionTransition(81);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 44:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(62);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 45:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(197);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 46:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(198);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 47:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(194);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 48:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(130);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 49:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(167);

                case 22:  // (
                    return new ParserActionTransition(82);

                case 34:  // =
                    return new ParserActionTransition(83);

                case 24:  // [
                    return new ParserActionTransition(84);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 50:
            switch (symbolID) {
                case 26:  // {
                case 33:  // ;
                    return new ParserActionReduce(140);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 51:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 48:  // in
                    return new ParserActionTransition(52);

                case 49:  // inout
                    return new ParserActionTransition(53);

                case 50:  // out
                    return new ParserActionTransition(54);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 52:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(155);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 53:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(157);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 54:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(156);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 55:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(142);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 56:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(149);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 57:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 58:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(153);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 59:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(158);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(90);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 60:
            switch (symbolID) {
                case 26:  // {
                case 33:  // ;
                    return new ParserActionReduce(141);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 61:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 47:  // const
                    return new ParserActionTransition(51);

                case 48:  // in
                    return new ParserActionTransition(52);

                case 49:  // inout
                    return new ParserActionTransition(53);

                case 50:  // out
                    return new ParserActionTransition(54);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 62:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(132);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 63:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(7);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 64:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 27:  // }
                    return new ParserActionTransition(153);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 65:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(155);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 66:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(133);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 67:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(186);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 68:
            switch (symbolID) {
                case 48:  // in
                    return new ParserActionTransition(17);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 69:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(188);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 70:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(157);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 71:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(185);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 72:
            switch (symbolID) {
                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 73:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(201);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 74:
            switch (symbolID) {
                case 26:  // {
                    return new ParserActionTransition(160);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 75:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 76:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(3);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 77:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(139);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 78:
            switch (symbolID) {
                case 26:  // {
                    return new ParserActionTransition(165);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 79:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(174);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 80:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 25:  // ]
                    return new ParserActionTransition(166);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 81:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(131);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 82:
            switch (symbolID) {
                case 23:  // )
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 47:  // const
                case 60:  // highp
                case 48:  // in
                case 49:  // inout
                case 59:  // lowp
                case 61:  // mediump
                case 50:  // out
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(144);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 83:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 84:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 25:  // ]
                    return new ParserActionTransition(173);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 85:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(147);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 86:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 87:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(151);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 88:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(150);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 89:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(154);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 90:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(145);

                case 24:  // [
                    return new ParserActionTransition(177);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 91:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(143);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 92:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 93:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 94:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 95:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 96:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 97:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 98:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(27);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 99:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(184);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 100:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(54);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 101:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(51);

                case 22:  // (
                    return new ParserActionTransition(185);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 102:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(52);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 103:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(53);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 104:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(89);

                case 38:  // +
                    return new ParserActionTransition(186);

                case 36:  // -
                    return new ParserActionTransition(187);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 105:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(102);

                case 44:  // &
                    return new ParserActionTransition(188);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 106:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(126);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 107:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(189);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 108:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 109:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(14);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 110:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 111:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(114);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 112:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(191);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 113:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(19);

                case 150:  // TRAILING-ANNOTATION
                    return new ParserActionTransition(48);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 114:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(10);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 115:
            switch (symbolID) {
                case 32:  // :
                    return new ParserActionTransition(192);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 116:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(193);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 117:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 118:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(20);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 119:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(100);

                case 8:  // !=
                    return new ParserActionTransition(195);

                case 7:  // ==
                    return new ParserActionTransition(196);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 120:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(104);

                case 43:  // ^
                    return new ParserActionTransition(197);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 121:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 33:  // ;
                    return new ParserActionTransition(199);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 122:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(11);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 123:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(55);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 124:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(200);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 125:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(65);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 126:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 23:  // )
                    return new ParserActionTransition(201);

                case 78:  // void
                    return new ParserActionTransition(203);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 127:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionTransition(204);

                case 31:  // ,
                    return new ParserActionTransition(205);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 128:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(60);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 129:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(206);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 130:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(106);

                case 42:  // |
                    return new ParserActionTransition(207);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 131:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(15);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 132:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(16);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 133:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(108);

                case 9:  // &&
                    return new ParserActionTransition(208);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 134:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(112);

                case 45:  // ?
                    return new ParserActionTransition(209);

                case 10:  // ||
                    return new ParserActionTransition(210);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 135:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(110);

                case 11:  // ^^
                    return new ParserActionTransition(211);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 136:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(86);

                case 41:  // %
                    return new ParserActionTransition(212);

                case 39:  // *
                    return new ParserActionTransition(213);

                case 40:  // /
                    return new ParserActionTransition(214);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 137:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(75);

                case 3:  // ++
                    return new ParserActionTransition(215);

                case 4:  // --
                    return new ParserActionTransition(216);

                case 30:  // .
                    return new ParserActionTransition(217);

                case 24:  // [
                    return new ParserActionTransition(218);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 138:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(58);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 139:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(97);

                case 28:  // <
                    return new ParserActionTransition(219);

                case 5:  // <=
                    return new ParserActionTransition(220);

                case 29:  // >
                    return new ParserActionTransition(221);

                case 6:  // >=
                    return new ParserActionTransition(222);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 140:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 33:  // ;
                    return new ParserActionTransition(223);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 141:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(12);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 142:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(92);

                case 1:  // <<
                    return new ParserActionTransition(225);

                case 2:  // >>
                    return new ParserActionTransition(226);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 143:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 144:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                    return new ParserActionReduce(25);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 145:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 27:  // }
                    return new ParserActionTransition(228);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 146:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(229);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 147:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(13);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 148:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(56);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 149:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(173);

                case 22:  // (
                    return new ParserActionTransition(230);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 150:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(82);

                case 18:  // %=
                    return new ParserActionTransition(231);

                case 19:  // &=
                    return new ParserActionTransition(232);

                case 16:  // *=
                    return new ParserActionTransition(233);

                case 14:  // +=
                    return new ParserActionTransition(234);

                case 17:  // -=
                    return new ParserActionTransition(235);

                case 15:  // /=
                    return new ParserActionTransition(236);

                case 12:  // <<=
                    return new ParserActionTransition(237);

                case 34:  // =
                    return new ParserActionTransition(238);

                case 13:  // >>=
                    return new ParserActionTransition(239);

                case 21:  // ^=
                    return new ParserActionTransition(240);

                case 20:  // |=
                    return new ParserActionTransition(241);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 151:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(242);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 152:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 27:  // }
                    return new ParserActionTransition(244);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 153:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(23);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 154:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 155:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(160);

                case 34:  // =
                    return new ParserActionTransition(246);

                case 24:  // [
                    return new ParserActionTransition(247);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 156:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(189);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 157:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(179);

                case 34:  // =
                    return new ParserActionTransition(248);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 158:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionTransition(249);

                case 31:  // ,
                    return new ParserActionTransition(250);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 159:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(251);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 160:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 161:
            switch (symbolID) {
                case 27:  // }
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(214);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 162:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 27:  // }
                    return new ParserActionTransition(254);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 163:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 164:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(256);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 165:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 166:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(206);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 167:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionReduce(128);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 168:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionTransition(260);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 169:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(230);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 170:
            switch (symbolID) {
                case 25:  // ]
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(82);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 171:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                    return new ParserActionReduce(223);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 172:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(172);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 173:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(168);

                case 34:  // =
                    return new ParserActionTransition(261);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 174:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionTransition(262);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 175:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(148);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 176:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(152);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 177:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 178:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(80);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 179:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 23:  // )
                    return new ParserActionTransition(264);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 180:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(78);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 181:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(76);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 182:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(79);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 183:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(77);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 184:
            switch (symbolID) {
                case 78:  // void
                case 23:  // )
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 76:  // false
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 75:  // struct
                case 77:  // true
                case 37:  // ~
                    return new ParserActionReduce(74);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 185:
            switch (symbolID) {
                case 78:  // void
                case 23:  // )
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 76:  // false
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 75:  // struct
                case 77:  // true
                case 37:  // ~
                    return new ParserActionReduce(73);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 186:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 187:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 188:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 189:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(47);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 190:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 32:  // :
                    return new ParserActionTransition(268);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 191:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(46);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 192:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(36);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 193:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(50);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 194:
            switch (symbolID) {
                case 67:  // while
                    return new ParserActionTransition(269);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 195:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 196:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 197:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 198:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 199:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(28);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 200:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 201:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(68);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 202:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(70);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 203:
            switch (symbolID) {
                case 22:  // (
                case 24:  // [
                    return new ParserActionReduce(208);

                case 23:  // )
                    return new ParserActionTransition(277);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 204:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(69);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 205:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 206:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 207:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 208:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 209:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 210:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 211:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 212:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 213:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 214:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 215:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(62);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 216:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(63);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 217:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(288);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(289);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 218:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 219:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 220:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 221:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 222:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 223:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(48);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 224:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 33:  // ;
                    return new ParserActionTransition(297);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 225:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 226:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 227:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                    return new ParserActionReduce(26);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 228:
            switch (symbolID) {
                case -1:  // $END$
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(24);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 229:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 230:
            switch (symbolID) {
                case 78:  // void
                case 23:  // )
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 76:  // false
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 75:  // struct
                case 77:  // true
                case 37:  // ~
                    return new ParserActionReduce(72);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 231:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 232:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 233:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 234:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 235:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 236:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 237:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 238:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 239:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 240:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 241:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 242:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 243:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 27:  // }
                    return new ParserActionTransition(316);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 244:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(17);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 245:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(81);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 246:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 247:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 25:  // ]
                    return new ParserActionTransition(318);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 248:
            switch (symbolID) {
                case 145:  // INTCONSTANT
                    return new ParserActionTransition(320);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 249:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 48:  // in
                case 50:  // out
                case 51:  // uniform
                case 52:  // varying
                case 33:  // ;
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 60:  // highp
                case 59:  // lowp
                case 61:  // mediump
                case 75:  // struct
                case 78:  // void
                    return new ParserActionReduce(178);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 250:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(321);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 251:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(134);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 252:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 27:  // }
                    return new ParserActionTransition(322);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 253:
            switch (symbolID) {
                case 27:  // }
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(215);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 254:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 24:  // [
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(213);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 255:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(256);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 256:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(220);

                case 24:  // [
                    return new ParserActionTransition(324);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 257:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(218);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 258:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(325);

                case 33:  // ;
                    return new ParserActionTransition(326);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 259:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 27:  // }
                    return new ParserActionTransition(327);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 260:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(207);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 261:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 262:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(169);

                case 34:  // =
                    return new ParserActionTransition(329);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 263:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionTransition(330);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 264:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(57);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 265:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(87);

                case 41:  // %
                    return new ParserActionTransition(212);

                case 39:  // *
                    return new ParserActionTransition(213);

                case 40:  // /
                    return new ParserActionTransition(214);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 266:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(88);

                case 41:  // %
                    return new ParserActionTransition(212);

                case 39:  // *
                    return new ParserActionTransition(213);

                case 40:  // /
                    return new ParserActionTransition(214);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 267:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(101);

                case 8:  // !=
                    return new ParserActionTransition(195);

                case 7:  // ==
                    return new ParserActionTransition(196);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 268:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(35);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 269:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(331);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 270:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(99);

                case 28:  // <
                    return new ParserActionTransition(219);

                case 5:  // <=
                    return new ParserActionTransition(220);

                case 29:  // >
                    return new ParserActionTransition(221);

                case 6:  // >=
                    return new ParserActionTransition(222);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 271:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(98);

                case 28:  // <
                    return new ParserActionTransition(219);

                case 5:  // <=
                    return new ParserActionTransition(220);

                case 29:  // >
                    return new ParserActionTransition(221);

                case 6:  // >=
                    return new ParserActionTransition(222);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 272:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(103);

                case 44:  // &
                    return new ParserActionTransition(188);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 273:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(127);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 274:
            switch (symbolID) {
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 76:  // false
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 57:  // smooth
                case 75:  // struct
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 37:  // ~
                    return new ParserActionReduce(41);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 275:
            switch (symbolID) {
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 76:  // false
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 57:  // smooth
                case 75:  // struct
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 37:  // ~
                    return new ParserActionReduce(40);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 276:
            switch (symbolID) {
                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 33:  // ;
                    return new ParserActionTransition(332);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 277:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(67);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 278:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(71);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 279:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 23:  // )
                    return new ParserActionTransition(335);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 280:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(105);

                case 43:  // ^
                    return new ParserActionTransition(197);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 281:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(107);

                case 42:  // |
                    return new ParserActionTransition(207);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 282:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 32:  // :
                    return new ParserActionTransition(336);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 283:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(111);

                case 11:  // ^^
                    return new ParserActionTransition(211);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 284:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(109);

                case 9:  // &&
                    return new ParserActionTransition(208);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 285:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(85);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 286:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(83);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 287:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(84);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 288:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(61);

                case 22:  // (
                    return new ParserActionTransition(184);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 289:
            switch (symbolID) {
                case 22:  // (
                    return new ParserActionTransition(185);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 290:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(66);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 291:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionReduce(64);

                case 31:  // ,
                    return new ParserActionTransition(198);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 292:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionTransition(337);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 293:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(93);

                case 1:  // <<
                    return new ParserActionTransition(225);

                case 2:  // >>
                    return new ParserActionTransition(226);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 294:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(95);

                case 1:  // <<
                    return new ParserActionTransition(225);

                case 2:  // >>
                    return new ParserActionTransition(226);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 295:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(94);

                case 1:  // <<
                    return new ParserActionTransition(225);

                case 2:  // >>
                    return new ParserActionTransition(226);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 296:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(96);

                case 1:  // <<
                    return new ParserActionTransition(225);

                case 2:  // >>
                    return new ParserActionTransition(226);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 297:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(49);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 298:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(90);

                case 38:  // +
                    return new ParserActionTransition(186);

                case 36:  // -
                    return new ParserActionTransition(187);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 299:
            switch (symbolID) {
                case 33:  // ;
                case 45:  // ?
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(91);

                case 38:  // +
                    return new ParserActionTransition(186);

                case 36:  // -
                    return new ParserActionTransition(187);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 300:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 23:  // )
                    return new ParserActionTransition(338);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 301:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(118);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 302:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(123);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 303:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(116);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 304:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(119);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 305:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(120);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 306:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(117);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 307:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(121);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 308:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(115);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 309:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(122);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 310:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(124);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 311:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 32:  // :
                case 25:  // ]
                    return new ParserActionReduce(125);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 312:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionTransition(339);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 313:
            switch (symbolID) {
                case 23:  // )
                case 33:  // ;
                    return new ParserActionReduce(31);

                case 31:  // ,
                    return new ParserActionTransition(198);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 314:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(340);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 315:
            switch (symbolID) {
                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 316:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(18);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 317:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(165);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 318:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(161);

                case 34:  // =
                    return new ParserActionTransition(341);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 319:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionTransition(342);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 320:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(180);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 321:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(181);

                case 34:  // =
                    return new ParserActionTransition(343);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 322:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                case 24:  // [
                case 33:  // ;
                case 31:  // ,
                case 23:  // )
                case 22:  // (
                    return new ParserActionReduce(212);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 323:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(325);

                case 33:  // ;
                    return new ParserActionTransition(344);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 324:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 25:  // ]
                    return new ParserActionTransition(345);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 325:
            switch (symbolID) {
                case 142:  // IDENTIFIER
                    return new ParserActionTransition(256);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 326:
            switch (symbolID) {
                case 27:  // }
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(216);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 327:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(348);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(349);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 328:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(170);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 329:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 330:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(146);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 331:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 332:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionReduce(42);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 333:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(353);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 334:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionTransition(354);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 335:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 336:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 337:
            switch (symbolID) {
                case 34:  // =
                case 16:  // *=
                case 15:  // /=
                case 18:  // %=
                case 14:  // +=
                case 17:  // -=
                case 12:  // <<=
                case 13:  // >>=
                case 19:  // &=
                case 21:  // ^=
                case 20:  // |=
                case 24:  // [
                case 30:  // .
                case 3:  // ++
                case 4:  // --
                case 33:  // ;
                case 45:  // ?
                case 39:  // *
                case 40:  // /
                case 41:  // %
                case 38:  // +
                case 36:  // -
                case 1:  // <<
                case 2:  // >>
                case 28:  // <
                case 29:  // >
                case 5:  // <=
                case 6:  // >=
                case 7:  // ==
                case 8:  // !=
                case 44:  // &
                case 43:  // ^
                case 42:  // |
                case 9:  // &&
                case 11:  // ^^
                case 10:  // ||
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(59);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 338:
            switch (symbolID) {
                case 26:  // {
                    return new ParserActionTransition(357);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 339:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 26:  // {
                    return new ParserActionTransition(64);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 340:
            switch (symbolID) {
                case 34:  // =
                    return new ParserActionTransition(361);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 341:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 342:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(162);

                case 34:  // =
                    return new ParserActionTransition(363);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 343:
            switch (symbolID) {
                case 145:  // INTCONSTANT
                    return new ParserActionTransition(364);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 344:
            switch (symbolID) {
                case 27:  // }
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                    return new ParserActionReduce(217);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 345:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(221);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 346:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionTransition(365);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 347:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(219);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 348:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(135);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 349:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(366);

                case 24:  // [
                    return new ParserActionTransition(367);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 350:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(171);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 351:
            switch (symbolID) {
                case 31:  // ,
                    return new ParserActionTransition(198);

                case 23:  // )
                    return new ParserActionTransition(368);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 352:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionReduce(44);

                case 31:  // ,
                    return new ParserActionTransition(198);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 353:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionReduce(43);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 354:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 26:  // {
                    return new ParserActionTransition(64);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 355:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                    return new ParserActionReduce(29);

                case 64:  // else
                    return new ParserActionTransition(371);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 356:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                case 25:  // ]
                case 23:  // )
                case 32:  // :
                    return new ParserActionReduce(113);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 357:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 27:  // }
                    return new ParserActionTransition(373);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 358:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(21);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 359:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(22);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 360:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(37);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 361:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 362:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(163);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 363:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 364:
            switch (symbolID) {
                case 23:  // )
                case 31:  // ,
                    return new ParserActionReduce(182);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 365:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(222);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 366:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(136);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 367:
            switch (symbolID) {
                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 25:  // ]
                    return new ParserActionTransition(376);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 368:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(378);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 369:
            switch (symbolID) {
                case 23:  // )
                    return new ParserActionReduce(45);

                case 31:  // ,
                    return new ParserActionTransition(198);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 370:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(39);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 371:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 372:
            switch (symbolID) {
                case 149:  // LEADING-ANNOTATION
                    return new ParserActionTransition(1);

                case 56:  // flat
                    return new ParserActionTransition(10);

                case 144:  // FLOATCONSTANT
                    return new ParserActionTransition(100);

                case 142:  // IDENTIFIER
                    return new ParserActionTransition(101);

                case 145:  // INTCONSTANT
                    return new ParserActionTransition(102);

                case 146:  // UINTCONSTANT
                    return new ParserActionTransition(103);

                case 71:  // break
                    return new ParserActionTransition(107);

                case 68:  // case
                    return new ParserActionTransition(108);

                case 72:  // continue
                    return new ParserActionTransition(112);

                case 69:  // default
                    return new ParserActionTransition(115);

                case 73:  // discard
                    return new ParserActionTransition(116);

                case 63:  // do
                    return new ParserActionTransition(117);

                case 76:  // false
                    return new ParserActionTransition(123);

                case 65:  // for
                    return new ParserActionTransition(124);

                case 66:  // if
                    return new ParserActionTransition(129);

                case 74:  // return
                    return new ParserActionTransition(140);

                case 70:  // switch
                    return new ParserActionTransition(146);

                case 77:  // true
                    return new ParserActionTransition(148);

                case 67:  // while
                    return new ParserActionTransition(151);

                case 26:  // {
                    return new ParserActionTransition(152);

                case 37:  // ~
                    return new ParserActionTransition(154);

                case 60:  // highp
                    return new ParserActionTransition(16);

                case 48:  // in
                    return new ParserActionTransition(17);

                case 148:  // NATIVE-TYPE
                    return new ParserActionTransition(2);

                case 54:  // invariant
                    return new ParserActionTransition(20);

                case 55:  // layout
                    return new ParserActionTransition(21);

                case 59:  // lowp
                    return new ParserActionTransition(23);

                case 61:  // mediump
                    return new ParserActionTransition(24);

                case 58:  // noperspective
                    return new ParserActionTransition(25);

                case 50:  // out
                    return new ParserActionTransition(26);

                case 62:  // precision
                    return new ParserActionTransition(27);

                case 143:  // TYPE-NAME
                    return new ParserActionTransition(3);

                case 57:  // smooth
                    return new ParserActionTransition(31);

                case 75:  // struct
                    return new ParserActionTransition(33);

                case 27:  // }
                    return new ParserActionTransition(380);

                case 46:  // attribute
                    return new ParserActionTransition(4);

                case 51:  // uniform
                    return new ParserActionTransition(40);

                case 52:  // varying
                    return new ParserActionTransition(41);

                case 78:  // void
                    return new ParserActionTransition(42);

                case 53:  // centroid
                    return new ParserActionTransition(5);

                case 47:  // const
                    return new ParserActionTransition(6);

                case 35:  // !
                    return new ParserActionTransition(92);

                case 22:  // (
                    return new ParserActionTransition(93);

                case 38:  // +
                    return new ParserActionTransition(94);

                case 3:  // ++
                    return new ParserActionTransition(95);

                case 36:  // -
                    return new ParserActionTransition(96);

                case 4:  // --
                    return new ParserActionTransition(97);

                case 33:  // ;
                    return new ParserActionTransition(98);

                case 147:  // FIELD-SELECTION
                    return new ParserActionTransition(99);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 373:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(34);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 374:
            switch (symbolID) {
                case 23:  // )
                case 33:  // ;
                    return new ParserActionReduce(32);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 375:
            switch (symbolID) {
                case 33:  // ;
                case 31:  // ,
                    return new ParserActionReduce(164);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 376:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(381);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 377:
            switch (symbolID) {
                case 25:  // ]
                    return new ParserActionTransition(382);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 378:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(38);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 379:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(30);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 380:
            switch (symbolID) {
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 146:  // UINTCONSTANT
                case 46:  // attribute
                case 71:  // break
                case 68:  // case
                case 53:  // centroid
                case 47:  // const
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 56:  // flat
                case 65:  // for
                case 60:  // highp
                case 66:  // if
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 57:  // smooth
                case 75:  // struct
                case 70:  // switch
                case 77:  // true
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(33);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 381:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(137);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 382:
            switch (symbolID) {
                case 33:  // ;
                    return new ParserActionTransition(383);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 383:
            switch (symbolID) {
                case -1:  // $END$
                case 150:  // TRAILING-ANNOTATION
                case 149:  // LEADING-ANNOTATION
                case 148:  // NATIVE-TYPE
                case 143:  // TYPE-NAME
                case 46:  // attribute
                case 53:  // centroid
                case 47:  // const
                case 56:  // flat
                case 60:  // highp
                case 48:  // in
                case 54:  // invariant
                case 55:  // layout
                case 59:  // lowp
                case 61:  // mediump
                case 58:  // noperspective
                case 50:  // out
                case 62:  // precision
                case 57:  // smooth
                case 75:  // struct
                case 51:  // uniform
                case 52:  // varying
                case 78:  // void
                case 27:  // }
                case 147:  // FIELD-SELECTION
                case 144:  // FLOATCONSTANT
                case 142:  // IDENTIFIER
                case 145:  // INTCONSTANT
                case 146:  // UINTCONSTANT
                case 71:  // break
                case 68:  // case
                case 72:  // continue
                case 69:  // default
                case 73:  // discard
                case 63:  // do
                case 76:  // false
                case 65:  // for
                case 66:  // if
                case 35:  // !
                case 22:  // (
                case 38:  // +
                case 3:  // ++
                case 36:  // -
                case 4:  // --
                case 33:  // ;
                case 74:  // return
                case 70:  // switch
                case 77:  // true
                case 67:  // while
                case 26:  // {
                case 37:  // ~
                case 64:  // else
                    return new ParserActionReduce(138);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        default:
            return new ParserActionError();

    } // switch (stateID

    return 0;
}   // action()

int GLSLParser::gotoState(const int stateID, const int symbolID) const {
    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 174:  // [fully-specified-type]
                    return 11;

                case 155:  // [function-definition]
                    return 12;

                case 206:  // [function-header]
                    return 13;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 158:  // [function-prototype]
                    return 15;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 152:  // [program]
                    return 29;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 153:  // [translation-unit]
                    return 35;

                case 204:  // [type-qualifier]
                    return 36;

                case 186:  // [type-specifier]
                    return 37;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 156:  // [declaration]
                    return 7;

                case 157:  // [enhanced-declaration]
                    return 8;

                case 154:  // [external-declaration]
                    return 9;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 1:
            switch (symbolID) {
                case 174:  // [fully-specified-type]
                    return 11;

                case 206:  // [function-header]
                    return 13;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 186:  // [type-specifier]
                    return 37;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 156:  // [declaration]
                    return 43;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 13:
            switch (symbolID) {
                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 208:  // [parameter-declaration]
                    return 55;

                case 209:  // [parameter-declarator]
                    return 56;

                case 210:  // [parameter-qualifier]
                    return 57;

                case 211:  // [parameter-type-specifier]
                    return 58;

                case 186:  // [type-specifier]
                    return 59;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 15:
            switch (symbolID) {
                case 159:  // [compound-statement-no-new-scope]
                    return 63;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 19:
            switch (symbolID) {
                case 216:  // [storage-qualifier]
                    return 67;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 20:
            switch (symbolID) {
                case 213:  // [interpolation-qualifier]
                    return 68;

                case 216:  // [storage-qualifier]
                    return 69;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 22:
            switch (symbolID) {
                case 216:  // [storage-qualifier]
                    return 71;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 27:
            switch (symbolID) {
                case 202:  // [precision-qualifier]
                    return 72;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 28:
            switch (symbolID) {
                case 218:  // [struct-specifier]
                    return 34;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 203:  // [type-specifier-no-prec]
                    return 73;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 35:
            switch (symbolID) {
                case 174:  // [fully-specified-type]
                    return 11;

                case 155:  // [function-definition]
                    return 12;

                case 206:  // [function-header]
                    return 13;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 158:  // [function-prototype]
                    return 15;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 186:  // [type-specifier]
                    return 37;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 154:  // [external-declaration]
                    return 76;

                case 156:  // [declaration]
                    return 7;

                case 157:  // [enhanced-declaration]
                    return 8;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 36:
            switch (symbolID) {
                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 186:  // [type-specifier]
                    return 79;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 51:
            switch (symbolID) {
                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 186:  // [type-specifier]
                    return 59;

                case 209:  // [parameter-declarator]
                    return 85;

                case 210:  // [parameter-qualifier]
                    return 86;

                case 211:  // [parameter-type-specifier]
                    return 87;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 57:
            switch (symbolID) {
                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 186:  // [type-specifier]
                    return 59;

                case 209:  // [parameter-declarator]
                    return 88;

                case 211:  // [parameter-type-specifier]
                    return 89;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 61:
            switch (symbolID) {
                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 209:  // [parameter-declarator]
                    return 56;

                case 210:  // [parameter-qualifier]
                    return 57;

                case 211:  // [parameter-type-specifier]
                    return 58;

                case 186:  // [type-specifier]
                    return 59;

                case 208:  // [parameter-declaration]
                    return 91;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 64:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 161:  // [statement]
                    return 144;

                case 170:  // [statement-list]
                    return 145;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 68:
            switch (symbolID) {
                case 216:  // [storage-qualifier]
                    return 156;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 70:
            switch (symbolID) {
                case 214:  // [layout-qualifier-id-list]
                    return 158;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 72:
            switch (symbolID) {
                case 203:  // [type-specifier-no-prec]
                    return 159;

                case 218:  // [struct-specifier]
                    return 34;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 75:
            switch (symbolID) {
                case 219:  // [struct-declaration]
                    return 161;

                case 205:  // [struct-declaration-list]
                    return 162;

                case 204:  // [type-qualifier]
                    return 163;

                case 186:  // [type-specifier]
                    return 164;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 80:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 199:  // [conditional-expression]
                    return 167;

                case 200:  // [constant-expression]
                    return 168;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 83:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 171;

                case 175:  // [initializer]
                    return 172;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 84:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 199:  // [conditional-expression]
                    return 167;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 200:  // [constant-expression]
                    return 174;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 86:
            switch (symbolID) {
                case 209:  // [parameter-declarator]
                    return 175;

                case 211:  // [parameter-type-specifier]
                    return 176;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 186:  // [type-specifier]
                    return 59;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 92:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 178;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 93:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 172:  // [expression]
                    return 179;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 94:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 180;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 95:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 181;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 96:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 182;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 97:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 183;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 108:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 172:  // [expression]
                    return 190;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 117:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 161:  // [statement]
                    return 194;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 126:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 202;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 140:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 172:  // [expression]
                    return 224;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 145:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 161:  // [statement]
                    return 227;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 152:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 161:  // [statement]
                    return 144;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 170:  // [statement-list]
                    return 243;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 154:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 245;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 160:
            switch (symbolID) {
                case 219:  // [struct-declaration]
                    return 161;

                case 204:  // [type-qualifier]
                    return 163;

                case 186:  // [type-specifier]
                    return 164;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 205:  // [struct-declaration-list]
                    return 252;

                case 202:  // [precision-qualifier]
                    return 28;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 162:
            switch (symbolID) {
                case 204:  // [type-qualifier]
                    return 163;

                case 186:  // [type-specifier]
                    return 164;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 219:  // [struct-declaration]
                    return 253;

                case 202:  // [precision-qualifier]
                    return 28;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 163:
            switch (symbolID) {
                case 186:  // [type-specifier]
                    return 255;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 164:
            switch (symbolID) {
                case 221:  // [struct-declarator]
                    return 257;

                case 220:  // [struct-declarator-list]
                    return 258;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 165:
            switch (symbolID) {
                case 219:  // [struct-declaration]
                    return 161;

                case 204:  // [type-qualifier]
                    return 163;

                case 186:  // [type-specifier]
                    return 164;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 205:  // [struct-declaration-list]
                    return 259;

                case 202:  // [precision-qualifier]
                    return 28;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 177:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 199:  // [conditional-expression]
                    return 167;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 200:  // [constant-expression]
                    return 263;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 186:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 188:  // [multiplicative-expression]
                    return 265;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 187:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 188:  // [multiplicative-expression]
                    return 266;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 188:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 192:  // [equality-expression]
                    return 267;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 195:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 191:  // [relational-expression]
                    return 270;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 196:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 191:  // [relational-expression]
                    return 271;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 197:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 192:  // [equality-expression]
                    return 119;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 193:  // [and-expression]
                    return 272;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 198:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 273;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 200:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 163:  // [declaration-statement]
                    return 274;

                case 164:  // [expression-statement]
                    return 275;

                case 176:  // [for-init-statement]
                    return 276;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 205:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 278;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 206:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 172:  // [expression]
                    return 279;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 207:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 194:  // [exclusive-or-expression]
                    return 280;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 208:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 195:  // [inclusive-or-expression]
                    return 281;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 209:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 172:  // [expression]
                    return 282;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 210:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 197:  // [logical-xor-expression]
                    return 283;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 211:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 196:  // [logical-and-expression]
                    return 284;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 212:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 285;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 213:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 286;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 214:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 287;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 217:
            switch (symbolID) {
                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 182:  // [function-call-generic]
                    return 290;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 218:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 172:  // [expression]
                    return 291;

                case 180:  // [integer-expression]
                    return 292;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 219:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 190:  // [shift-expression]
                    return 293;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 220:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 190:  // [shift-expression]
                    return 294;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 221:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 190:  // [shift-expression]
                    return 295;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 222:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 190:  // [shift-expression]
                    return 296;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 225:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 189:  // [additive-expression]
                    return 298;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 226:
            switch (symbolID) {
                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 189:  // [additive-expression]
                    return 299;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 229:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 172:  // [expression]
                    return 300;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 231:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 301;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 232:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 302;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 233:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 303;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 234:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 304;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 235:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 305;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 236:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 306;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 237:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 307;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 238:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 308;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 239:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 309;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 240:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 310;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 241:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 185:  // [assignment-expression]
                    return 311;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 242:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 149;

                case 187:  // [unary-expression]
                    return 150;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 173:  // [condition]
                    return 312;

                case 172:  // [expression]
                    return 313;

                case 174:  // [fully-specified-type]
                    return 314;

                case 204:  // [type-qualifier]
                    return 315;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 243:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 161:  // [statement]
                    return 227;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 246:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 171;

                case 202:  // [precision-qualifier]
                    return 28;

                case 175:  // [initializer]
                    return 317;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 247:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 199:  // [conditional-expression]
                    return 167;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 200:  // [constant-expression]
                    return 319;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 252:
            switch (symbolID) {
                case 204:  // [type-qualifier]
                    return 163;

                case 186:  // [type-specifier]
                    return 164;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 219:  // [struct-declaration]
                    return 253;

                case 202:  // [precision-qualifier]
                    return 28;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 255:
            switch (symbolID) {
                case 221:  // [struct-declarator]
                    return 257;

                case 220:  // [struct-declarator-list]
                    return 323;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 259:
            switch (symbolID) {
                case 204:  // [type-qualifier]
                    return 163;

                case 186:  // [type-specifier]
                    return 164;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 219:  // [struct-declaration]
                    return 253;

                case 202:  // [precision-qualifier]
                    return 28;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 261:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 171;

                case 202:  // [precision-qualifier]
                    return 28;

                case 175:  // [initializer]
                    return 328;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 276:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 186:  // [type-specifier]
                    return 149;

                case 187:  // [unary-expression]
                    return 150;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 172:  // [expression]
                    return 313;

                case 174:  // [fully-specified-type]
                    return 314;

                case 204:  // [type-qualifier]
                    return 315;

                case 216:  // [storage-qualifier]
                    return 32;

                case 173:  // [condition]
                    return 333;

                case 177:  // [for-rest-statement]
                    return 334;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 315:
            switch (symbolID) {
                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 186:  // [type-specifier]
                    return 79;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 324:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 199:  // [conditional-expression]
                    return 167;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 200:  // [constant-expression]
                    return 346;

                case 218:  // [struct-specifier]
                    return 34;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 325:
            switch (symbolID) {
                case 221:  // [struct-declarator]
                    return 347;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 329:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 171;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 175:  // [initializer]
                    return 350;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 331:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 172:  // [expression]
                    return 351;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 332:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 172:  // [expression]
                    return 352;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 335:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 161:  // [statement]
                    return 355;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 336:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 185:  // [assignment-expression]
                    return 356;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 339:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 159:  // [compound-statement-no-new-scope]
                    return 358;

                case 162:  // [simple-statement]
                    return 359;

                case 171:  // [statement-no-new-scope]
                    return 360;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 341:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 171;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 175:  // [initializer]
                    return 362;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 353:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 172:  // [expression]
                    return 369;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 354:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 159:  // [compound-statement-no-new-scope]
                    return 358;

                case 162:  // [simple-statement]
                    return 359;

                case 204:  // [type-qualifier]
                    return 36;

                case 171:  // [statement-no-new-scope]
                    return 370;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 357:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 161:  // [statement]
                    return 144;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 170:  // [statement-list]
                    return 372;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 361:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 171;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 175:  // [initializer]
                    return 374;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 363:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 199:  // [conditional-expression]
                    return 111;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 187:  // [unary-expression]
                    return 150;

                case 186:  // [type-specifier]
                    return 169;

                case 185:  // [assignment-expression]
                    return 171;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 175:  // [initializer]
                    return 375;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 367:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 192:  // [equality-expression]
                    return 119;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 190:  // [shift-expression]
                    return 142;

                case 199:  // [conditional-expression]
                    return 167;

                case 186:  // [type-specifier]
                    return 169;

                case 187:  // [unary-expression]
                    return 170;

                case 202:  // [precision-qualifier]
                    return 28;

                case 218:  // [struct-specifier]
                    return 34;

                case 200:  // [constant-expression]
                    return 377;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 371:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 161:  // [statement]
                    return 379;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 372:
            switch (symbolID) {
                case 189:  // [additive-expression]
                    return 104;

                case 193:  // [and-expression]
                    return 105;

                case 185:  // [assignment-expression]
                    return 106;

                case 167:  // [case-label]
                    return 109;

                case 160:  // [compound-statement]
                    return 110;

                case 199:  // [conditional-expression]
                    return 111;

                case 156:  // [declaration]
                    return 113;

                case 163:  // [declaration-statement]
                    return 114;

                case 157:  // [enhanced-declaration]
                    return 118;

                case 192:  // [equality-expression]
                    return 119;

                case 174:  // [fully-specified-type]
                    return 11;

                case 194:  // [exclusive-or-expression]
                    return 120;

                case 172:  // [expression]
                    return 121;

                case 164:  // [expression-statement]
                    return 122;

                case 182:  // [function-call-generic]
                    return 125;

                case 183:  // [function-call-header]
                    return 126;

                case 184:  // [function-call-header-with-parameters]
                    return 127;

                case 181:  // [function-call-or-method]
                    return 128;

                case 195:  // [inclusive-or-expression]
                    return 130;

                case 168:  // [iteration-statement]
                    return 131;

                case 169:  // [jump-statement]
                    return 132;

                case 196:  // [logical-and-expression]
                    return 133;

                case 198:  // [logical-or-expression]
                    return 134;

                case 197:  // [logical-xor-expression]
                    return 135;

                case 188:  // [multiplicative-expression]
                    return 136;

                case 179:  // [postfix-expression]
                    return 137;

                case 178:  // [primary-expression]
                    return 138;

                case 191:  // [relational-expression]
                    return 139;

                case 206:  // [function-header]
                    return 13;

                case 165:  // [selection-statement]
                    return 141;

                case 190:  // [shift-expression]
                    return 142;

                case 162:  // [simple-statement]
                    return 143;

                case 166:  // [switch-statement]
                    return 147;

                case 186:  // [type-specifier]
                    return 149;

                case 207:  // [function-header-with-parameters]
                    return 14;

                case 187:  // [unary-expression]
                    return 150;

                case 201:  // [init-declarator-list]
                    return 18;

                case 213:  // [interpolation-qualifier]
                    return 19;

                case 161:  // [statement]
                    return 227;

                case 215:  // [layout-qualifier]
                    return 22;

                case 202:  // [precision-qualifier]
                    return 28;

                case 212:  // [single-declaration]
                    return 30;

                case 216:  // [storage-qualifier]
                    return 32;

                case 218:  // [struct-specifier]
                    return 34;

                case 204:  // [type-qualifier]
                    return 36;

                case 203:  // [type-specifier-no-prec]
                    return 38;

                case 217:  // [type-specifier-nonarray]
                    return 39;

                case 158:  // [function-prototype]
                    return 44;

                default:
                    return -1;
            }   // switch (symbolID
            break;

    } // switch (stateID

    return -1;
}   // gotoState()

ProductionStub* GLSLParser::findProduction(const int productionID) const {
    switch (productionID) {
        case 0:  // [$START$] ::= [program]
            {
                int bodyIDs[] = {152};
                return new ProductionStub(0, bodyIDs, 1);
            }
        case 1:  // [program] ::= [translation-unit]
            {
                int bodyIDs[] = {153};
                return new ProductionStub(152, bodyIDs, 1);
            }
        case 2:  // [translation-unit] ::= [external-declaration]
            {
                int bodyIDs[] = {154};
                return new ProductionStub(153, bodyIDs, 1);
            }
        case 3:  // [translation-unit] ::= [translation-unit] [external-declaration]
            {
                int bodyIDs[] = {153, 154};
                return new ProductionStub(153, bodyIDs, 2);
            }
        case 4:  // [external-declaration] ::= [function-definition]
            {
                int bodyIDs[] = {155};
                return new ProductionStub(154, bodyIDs, 1);
            }
        case 5:  // [external-declaration] ::= [declaration]
            {
                int bodyIDs[] = {156};
                return new ProductionStub(154, bodyIDs, 1);
            }
        case 6:  // [external-declaration] ::= [enhanced-declaration]
            {
                int bodyIDs[] = {157};
                return new ProductionStub(154, bodyIDs, 1);
            }
        case 7:  // [function-definition] ::= [function-prototype] [compound-statement-no-new-scope]
            {
                int bodyIDs[] = {158, 159};
                return new ProductionStub(155, bodyIDs, 2);
            }
        case 8:  // [statement] ::= [compound-statement]
            {
                int bodyIDs[] = {160};
                return new ProductionStub(161, bodyIDs, 1);
            }
        case 9:  // [statement] ::= [simple-statement]
            {
                int bodyIDs[] = {162};
                return new ProductionStub(161, bodyIDs, 1);
            }
        case 10:  // [simple-statement] ::= [declaration-statement]
            {
                int bodyIDs[] = {163};
                return new ProductionStub(162, bodyIDs, 1);
            }
        case 11:  // [simple-statement] ::= [expression-statement]
            {
                int bodyIDs[] = {164};
                return new ProductionStub(162, bodyIDs, 1);
            }
        case 12:  // [simple-statement] ::= [selection-statement]
            {
                int bodyIDs[] = {165};
                return new ProductionStub(162, bodyIDs, 1);
            }
        case 13:  // [simple-statement] ::= [switch-statement]
            {
                int bodyIDs[] = {166};
                return new ProductionStub(162, bodyIDs, 1);
            }
        case 14:  // [simple-statement] ::= [case-label]
            {
                int bodyIDs[] = {167};
                return new ProductionStub(162, bodyIDs, 1);
            }
        case 15:  // [simple-statement] ::= [iteration-statement]
            {
                int bodyIDs[] = {168};
                return new ProductionStub(162, bodyIDs, 1);
            }
        case 16:  // [simple-statement] ::= [jump-statement]
            {
                int bodyIDs[] = {169};
                return new ProductionStub(162, bodyIDs, 1);
            }
        case 17:  // [compound-statement] ::= { }
            {
                int bodyIDs[] = {26, 27};
                return new ProductionStub(160, bodyIDs, 2);
            }
        case 18:  // [compound-statement] ::= { [statement-list] }
            {
                int bodyIDs[] = {26, 170, 27};
                return new ProductionStub(160, bodyIDs, 3);
            }
        case 19:  // [declaration-statement] ::= [declaration]
            {
                int bodyIDs[] = {156};
                return new ProductionStub(163, bodyIDs, 1);
            }
        case 20:  // [declaration-statement] ::= [enhanced-declaration]
            {
                int bodyIDs[] = {157};
                return new ProductionStub(163, bodyIDs, 1);
            }
        case 21:  // [statement-no-new-scope] ::= [compound-statement-no-new-scope]
            {
                int bodyIDs[] = {159};
                return new ProductionStub(171, bodyIDs, 1);
            }
        case 22:  // [statement-no-new-scope] ::= [simple-statement]
            {
                int bodyIDs[] = {162};
                return new ProductionStub(171, bodyIDs, 1);
            }
        case 23:  // [compound-statement-no-new-scope] ::= { }
            {
                int bodyIDs[] = {26, 27};
                return new ProductionStub(159, bodyIDs, 2);
            }
        case 24:  // [compound-statement-no-new-scope] ::= { [statement-list] }
            {
                int bodyIDs[] = {26, 170, 27};
                return new ProductionStub(159, bodyIDs, 3);
            }
        case 25:  // [statement-list] ::= [statement]
            {
                int bodyIDs[] = {161};
                return new ProductionStub(170, bodyIDs, 1);
            }
        case 26:  // [statement-list] ::= [statement-list] [statement]
            {
                int bodyIDs[] = {170, 161};
                return new ProductionStub(170, bodyIDs, 2);
            }
        case 27:  // [expression-statement] ::= ;
            {
                int bodyIDs[] = {33};
                return new ProductionStub(164, bodyIDs, 1);
            }
        case 28:  // [expression-statement] ::= [expression] ;
            {
                int bodyIDs[] = {172, 33};
                return new ProductionStub(164, bodyIDs, 2);
            }
        case 29:  // [selection-statement] ::= if ( [expression] ) [statement]
            {
                int bodyIDs[] = {66, 22, 172, 23, 161};
                return new ProductionStub(165, bodyIDs, 5);
            }
        case 30:  // [selection-statement] ::= if ( [expression] ) [statement] else [statement]
            {
                int bodyIDs[] = {66, 22, 172, 23, 161, 64, 161};
                return new ProductionStub(165, bodyIDs, 7);
            }
        case 31:  // [condition] ::= [expression]
            {
                int bodyIDs[] = {172};
                return new ProductionStub(173, bodyIDs, 1);
            }
        case 32:  // [condition] ::= [fully-specified-type] IDENTIFIER = [initializer]
            {
                int bodyIDs[] = {174, 142, 34, 175};
                return new ProductionStub(173, bodyIDs, 4);
            }
        case 33:  // [switch-statement] ::= switch ( [expression] ) { [statement-list] }
            {
                int bodyIDs[] = {70, 22, 172, 23, 26, 170, 27};
                return new ProductionStub(166, bodyIDs, 7);
            }
        case 34:  // [switch-statement] ::= switch ( [expression] ) { }
            {
                int bodyIDs[] = {70, 22, 172, 23, 26, 27};
                return new ProductionStub(166, bodyIDs, 6);
            }
        case 35:  // [case-label] ::= case [expression] :
            {
                int bodyIDs[] = {68, 172, 32};
                return new ProductionStub(167, bodyIDs, 3);
            }
        case 36:  // [case-label] ::= default :
            {
                int bodyIDs[] = {69, 32};
                return new ProductionStub(167, bodyIDs, 2);
            }
        case 37:  // [iteration-statement] ::= while ( [condition] ) [statement-no-new-scope]
            {
                int bodyIDs[] = {67, 22, 173, 23, 171};
                return new ProductionStub(168, bodyIDs, 5);
            }
        case 38:  // [iteration-statement] ::= do [statement] while ( [expression] ) ;
            {
                int bodyIDs[] = {63, 161, 67, 22, 172, 23, 33};
                return new ProductionStub(168, bodyIDs, 7);
            }
        case 39:  // [iteration-statement] ::= for ( [for-init-statement] [for-rest-statement] ) [statement-no-new-scope]
            {
                int bodyIDs[] = {65, 22, 176, 177, 23, 171};
                return new ProductionStub(168, bodyIDs, 6);
            }
        case 40:  // [for-init-statement] ::= [expression-statement]
            {
                int bodyIDs[] = {164};
                return new ProductionStub(176, bodyIDs, 1);
            }
        case 41:  // [for-init-statement] ::= [declaration-statement]
            {
                int bodyIDs[] = {163};
                return new ProductionStub(176, bodyIDs, 1);
            }
        case 42:  // [for-rest-statement] ::= ;
            {
                int bodyIDs[] = {33};
                return new ProductionStub(177, bodyIDs, 1);
            }
        case 43:  // [for-rest-statement] ::= [condition] ;
            {
                int bodyIDs[] = {173, 33};
                return new ProductionStub(177, bodyIDs, 2);
            }
        case 44:  // [for-rest-statement] ::= ; [expression]
            {
                int bodyIDs[] = {33, 172};
                return new ProductionStub(177, bodyIDs, 2);
            }
        case 45:  // [for-rest-statement] ::= [condition] ; [expression]
            {
                int bodyIDs[] = {173, 33, 172};
                return new ProductionStub(177, bodyIDs, 3);
            }
        case 46:  // [jump-statement] ::= continue ;
            {
                int bodyIDs[] = {72, 33};
                return new ProductionStub(169, bodyIDs, 2);
            }
        case 47:  // [jump-statement] ::= break ;
            {
                int bodyIDs[] = {71, 33};
                return new ProductionStub(169, bodyIDs, 2);
            }
        case 48:  // [jump-statement] ::= return ;
            {
                int bodyIDs[] = {74, 33};
                return new ProductionStub(169, bodyIDs, 2);
            }
        case 49:  // [jump-statement] ::= return [expression] ;
            {
                int bodyIDs[] = {74, 172, 33};
                return new ProductionStub(169, bodyIDs, 3);
            }
        case 50:  // [jump-statement] ::= discard ;
            {
                int bodyIDs[] = {73, 33};
                return new ProductionStub(169, bodyIDs, 2);
            }
        case 51:  // [primary-expression] ::= IDENTIFIER
            {
                int bodyIDs[] = {142};
                return new ProductionStub(178, bodyIDs, 1);
            }
        case 52:  // [primary-expression] ::= INTCONSTANT
            {
                int bodyIDs[] = {145};
                return new ProductionStub(178, bodyIDs, 1);
            }
        case 53:  // [primary-expression] ::= UINTCONSTANT
            {
                int bodyIDs[] = {146};
                return new ProductionStub(178, bodyIDs, 1);
            }
        case 54:  // [primary-expression] ::= FLOATCONSTANT
            {
                int bodyIDs[] = {144};
                return new ProductionStub(178, bodyIDs, 1);
            }
        case 55:  // [primary-expression] ::= false
            {
                int bodyIDs[] = {76};
                return new ProductionStub(178, bodyIDs, 1);
            }
        case 56:  // [primary-expression] ::= true
            {
                int bodyIDs[] = {77};
                return new ProductionStub(178, bodyIDs, 1);
            }
        case 57:  // [primary-expression] ::= ( [expression] )
            {
                int bodyIDs[] = {22, 172, 23};
                return new ProductionStub(178, bodyIDs, 3);
            }
        case 58:  // [postfix-expression] ::= [primary-expression]
            {
                int bodyIDs[] = {178};
                return new ProductionStub(179, bodyIDs, 1);
            }
        case 59:  // [postfix-expression] ::= [postfix-expression] [ [integer-expression] ]
            {
                int bodyIDs[] = {179, 24, 180, 25};
                return new ProductionStub(179, bodyIDs, 4);
            }
        case 60:  // [postfix-expression] ::= [function-call-or-method]
            {
                int bodyIDs[] = {181};
                return new ProductionStub(179, bodyIDs, 1);
            }
        case 61:  // [postfix-expression] ::= [postfix-expression] . FIELD-SELECTION
            {
                int bodyIDs[] = {179, 30, 147};
                return new ProductionStub(179, bodyIDs, 3);
            }
        case 62:  // [postfix-expression] ::= [postfix-expression] ++
            {
                int bodyIDs[] = {179, 3};
                return new ProductionStub(179, bodyIDs, 2);
            }
        case 63:  // [postfix-expression] ::= [postfix-expression] --
            {
                int bodyIDs[] = {179, 4};
                return new ProductionStub(179, bodyIDs, 2);
            }
        case 64:  // [integer-expression] ::= [expression]
            {
                int bodyIDs[] = {172};
                return new ProductionStub(180, bodyIDs, 1);
            }
        case 65:  // [function-call-or-method] ::= [function-call-generic]
            {
                int bodyIDs[] = {182};
                return new ProductionStub(181, bodyIDs, 1);
            }
        case 66:  // [function-call-or-method] ::= [postfix-expression] . [function-call-generic]
            {
                int bodyIDs[] = {179, 30, 182};
                return new ProductionStub(181, bodyIDs, 3);
            }
        case 67:  // [function-call-generic] ::= [function-call-header] void )
            {
                int bodyIDs[] = {183, 78, 23};
                return new ProductionStub(182, bodyIDs, 3);
            }
        case 68:  // [function-call-generic] ::= [function-call-header] )
            {
                int bodyIDs[] = {183, 23};
                return new ProductionStub(182, bodyIDs, 2);
            }
        case 69:  // [function-call-generic] ::= [function-call-header-with-parameters] )
            {
                int bodyIDs[] = {184, 23};
                return new ProductionStub(182, bodyIDs, 2);
            }
        case 70:  // [function-call-header-with-parameters] ::= [function-call-header] [assignment-expression]
            {
                int bodyIDs[] = {183, 185};
                return new ProductionStub(184, bodyIDs, 2);
            }
        case 71:  // [function-call-header-with-parameters] ::= [function-call-header-with-parameters] , [assignment-expression]
            {
                int bodyIDs[] = {184, 31, 185};
                return new ProductionStub(184, bodyIDs, 3);
            }
        case 72:  // [function-call-header] ::= [type-specifier] (
            {
                int bodyIDs[] = {186, 22};
                return new ProductionStub(183, bodyIDs, 2);
            }
        case 73:  // [function-call-header] ::= IDENTIFIER (
            {
                int bodyIDs[] = {142, 22};
                return new ProductionStub(183, bodyIDs, 2);
            }
        case 74:  // [function-call-header] ::= FIELD-SELECTION (
            {
                int bodyIDs[] = {147, 22};
                return new ProductionStub(183, bodyIDs, 2);
            }
        case 75:  // [unary-expression] ::= [postfix-expression]
            {
                int bodyIDs[] = {179};
                return new ProductionStub(187, bodyIDs, 1);
            }
        case 76:  // [unary-expression] ::= ++ [unary-expression]
            {
                int bodyIDs[] = {3, 187};
                return new ProductionStub(187, bodyIDs, 2);
            }
        case 77:  // [unary-expression] ::= -- [unary-expression]
            {
                int bodyIDs[] = {4, 187};
                return new ProductionStub(187, bodyIDs, 2);
            }
        case 78:  // [unary-expression] ::= + [unary-expression]
            {
                int bodyIDs[] = {38, 187};
                return new ProductionStub(187, bodyIDs, 2);
            }
        case 79:  // [unary-expression] ::= - [unary-expression]
            {
                int bodyIDs[] = {36, 187};
                return new ProductionStub(187, bodyIDs, 2);
            }
        case 80:  // [unary-expression] ::= ! [unary-expression]
            {
                int bodyIDs[] = {35, 187};
                return new ProductionStub(187, bodyIDs, 2);
            }
        case 81:  // [unary-expression] ::= ~ [unary-expression]
            {
                int bodyIDs[] = {37, 187};
                return new ProductionStub(187, bodyIDs, 2);
            }
        case 82:  // [multiplicative-expression] ::= [unary-expression]
            {
                int bodyIDs[] = {187};
                return new ProductionStub(188, bodyIDs, 1);
            }
        case 83:  // [multiplicative-expression] ::= [multiplicative-expression] * [unary-expression]
            {
                int bodyIDs[] = {188, 39, 187};
                return new ProductionStub(188, bodyIDs, 3);
            }
        case 84:  // [multiplicative-expression] ::= [multiplicative-expression] / [unary-expression]
            {
                int bodyIDs[] = {188, 40, 187};
                return new ProductionStub(188, bodyIDs, 3);
            }
        case 85:  // [multiplicative-expression] ::= [multiplicative-expression] % [unary-expression]
            {
                int bodyIDs[] = {188, 41, 187};
                return new ProductionStub(188, bodyIDs, 3);
            }
        case 86:  // [additive-expression] ::= [multiplicative-expression]
            {
                int bodyIDs[] = {188};
                return new ProductionStub(189, bodyIDs, 1);
            }
        case 87:  // [additive-expression] ::= [additive-expression] + [multiplicative-expression]
            {
                int bodyIDs[] = {189, 38, 188};
                return new ProductionStub(189, bodyIDs, 3);
            }
        case 88:  // [additive-expression] ::= [additive-expression] - [multiplicative-expression]
            {
                int bodyIDs[] = {189, 36, 188};
                return new ProductionStub(189, bodyIDs, 3);
            }
        case 89:  // [shift-expression] ::= [additive-expression]
            {
                int bodyIDs[] = {189};
                return new ProductionStub(190, bodyIDs, 1);
            }
        case 90:  // [shift-expression] ::= [shift-expression] << [additive-expression]
            {
                int bodyIDs[] = {190, 1, 189};
                return new ProductionStub(190, bodyIDs, 3);
            }
        case 91:  // [shift-expression] ::= [shift-expression] >> [additive-expression]
            {
                int bodyIDs[] = {190, 2, 189};
                return new ProductionStub(190, bodyIDs, 3);
            }
        case 92:  // [relational-expression] ::= [shift-expression]
            {
                int bodyIDs[] = {190};
                return new ProductionStub(191, bodyIDs, 1);
            }
        case 93:  // [relational-expression] ::= [relational-expression] < [shift-expression]
            {
                int bodyIDs[] = {191, 28, 190};
                return new ProductionStub(191, bodyIDs, 3);
            }
        case 94:  // [relational-expression] ::= [relational-expression] > [shift-expression]
            {
                int bodyIDs[] = {191, 29, 190};
                return new ProductionStub(191, bodyIDs, 3);
            }
        case 95:  // [relational-expression] ::= [relational-expression] <= [shift-expression]
            {
                int bodyIDs[] = {191, 5, 190};
                return new ProductionStub(191, bodyIDs, 3);
            }
        case 96:  // [relational-expression] ::= [relational-expression] >= [shift-expression]
            {
                int bodyIDs[] = {191, 6, 190};
                return new ProductionStub(191, bodyIDs, 3);
            }
        case 97:  // [equality-expression] ::= [relational-expression]
            {
                int bodyIDs[] = {191};
                return new ProductionStub(192, bodyIDs, 1);
            }
        case 98:  // [equality-expression] ::= [equality-expression] == [relational-expression]
            {
                int bodyIDs[] = {192, 7, 191};
                return new ProductionStub(192, bodyIDs, 3);
            }
        case 99:  // [equality-expression] ::= [equality-expression] != [relational-expression]
            {
                int bodyIDs[] = {192, 8, 191};
                return new ProductionStub(192, bodyIDs, 3);
            }
        case 100:  // [and-expression] ::= [equality-expression]
            {
                int bodyIDs[] = {192};
                return new ProductionStub(193, bodyIDs, 1);
            }
        case 101:  // [and-expression] ::= [and-expression] & [equality-expression]
            {
                int bodyIDs[] = {193, 44, 192};
                return new ProductionStub(193, bodyIDs, 3);
            }
        case 102:  // [exclusive-or-expression] ::= [and-expression]
            {
                int bodyIDs[] = {193};
                return new ProductionStub(194, bodyIDs, 1);
            }
        case 103:  // [exclusive-or-expression] ::= [exclusive-or-expression] ^ [and-expression]
            {
                int bodyIDs[] = {194, 43, 193};
                return new ProductionStub(194, bodyIDs, 3);
            }
        case 104:  // [inclusive-or-expression] ::= [exclusive-or-expression]
            {
                int bodyIDs[] = {194};
                return new ProductionStub(195, bodyIDs, 1);
            }
        case 105:  // [inclusive-or-expression] ::= [inclusive-or-expression] | [exclusive-or-expression]
            {
                int bodyIDs[] = {195, 42, 194};
                return new ProductionStub(195, bodyIDs, 3);
            }
        case 106:  // [logical-and-expression] ::= [inclusive-or-expression]
            {
                int bodyIDs[] = {195};
                return new ProductionStub(196, bodyIDs, 1);
            }
        case 107:  // [logical-and-expression] ::= [logical-and-expression] && [inclusive-or-expression]
            {
                int bodyIDs[] = {196, 9, 195};
                return new ProductionStub(196, bodyIDs, 3);
            }
        case 108:  // [logical-xor-expression] ::= [logical-and-expression]
            {
                int bodyIDs[] = {196};
                return new ProductionStub(197, bodyIDs, 1);
            }
        case 109:  // [logical-xor-expression] ::= [logical-xor-expression] ^^ [logical-and-expression]
            {
                int bodyIDs[] = {197, 11, 196};
                return new ProductionStub(197, bodyIDs, 3);
            }
        case 110:  // [logical-or-expression] ::= [logical-xor-expression]
            {
                int bodyIDs[] = {197};
                return new ProductionStub(198, bodyIDs, 1);
            }
        case 111:  // [logical-or-expression] ::= [logical-or-expression] || [logical-xor-expression]
            {
                int bodyIDs[] = {198, 10, 197};
                return new ProductionStub(198, bodyIDs, 3);
            }
        case 112:  // [conditional-expression] ::= [logical-or-expression]
            {
                int bodyIDs[] = {198};
                return new ProductionStub(199, bodyIDs, 1);
            }
        case 113:  // [conditional-expression] ::= [logical-or-expression] ? [expression] : [assignment-expression]
            {
                int bodyIDs[] = {198, 45, 172, 32, 185};
                return new ProductionStub(199, bodyIDs, 5);
            }
        case 114:  // [assignment-expression] ::= [conditional-expression]
            {
                int bodyIDs[] = {199};
                return new ProductionStub(185, bodyIDs, 1);
            }
        case 115:  // [assignment-expression] ::= [unary-expression] = [assignment-expression]
            {
                int bodyIDs[] = {187, 34, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 116:  // [assignment-expression] ::= [unary-expression] *= [assignment-expression]
            {
                int bodyIDs[] = {187, 16, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 117:  // [assignment-expression] ::= [unary-expression] /= [assignment-expression]
            {
                int bodyIDs[] = {187, 15, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 118:  // [assignment-expression] ::= [unary-expression] %= [assignment-expression]
            {
                int bodyIDs[] = {187, 18, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 119:  // [assignment-expression] ::= [unary-expression] += [assignment-expression]
            {
                int bodyIDs[] = {187, 14, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 120:  // [assignment-expression] ::= [unary-expression] -= [assignment-expression]
            {
                int bodyIDs[] = {187, 17, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 121:  // [assignment-expression] ::= [unary-expression] <<= [assignment-expression]
            {
                int bodyIDs[] = {187, 12, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 122:  // [assignment-expression] ::= [unary-expression] >>= [assignment-expression]
            {
                int bodyIDs[] = {187, 13, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 123:  // [assignment-expression] ::= [unary-expression] &= [assignment-expression]
            {
                int bodyIDs[] = {187, 19, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 124:  // [assignment-expression] ::= [unary-expression] ^= [assignment-expression]
            {
                int bodyIDs[] = {187, 21, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 125:  // [assignment-expression] ::= [unary-expression] |= [assignment-expression]
            {
                int bodyIDs[] = {187, 20, 185};
                return new ProductionStub(185, bodyIDs, 3);
            }
        case 126:  // [expression] ::= [assignment-expression]
            {
                int bodyIDs[] = {185};
                return new ProductionStub(172, bodyIDs, 1);
            }
        case 127:  // [expression] ::= [expression] , [assignment-expression]
            {
                int bodyIDs[] = {172, 31, 185};
                return new ProductionStub(172, bodyIDs, 3);
            }
        case 128:  // [constant-expression] ::= [conditional-expression]
            {
                int bodyIDs[] = {199};
                return new ProductionStub(200, bodyIDs, 1);
            }
        case 129:  // [enhanced-declaration] ::= LEADING-ANNOTATION [declaration]
            {
                int bodyIDs[] = {149, 156};
                return new ProductionStub(157, bodyIDs, 2);
            }
        case 130:  // [enhanced-declaration] ::= [declaration] TRAILING-ANNOTATION
            {
                int bodyIDs[] = {156, 150};
                return new ProductionStub(157, bodyIDs, 2);
            }
        case 131:  // [enhanced-declaration] ::= LEADING-ANNOTATION [declaration] TRAILING-ANNOTATION
            {
                int bodyIDs[] = {149, 156, 150};
                return new ProductionStub(157, bodyIDs, 3);
            }
        case 132:  // [declaration] ::= [function-prototype] ;
            {
                int bodyIDs[] = {158, 33};
                return new ProductionStub(156, bodyIDs, 2);
            }
        case 133:  // [declaration] ::= [init-declarator-list] ;
            {
                int bodyIDs[] = {201, 33};
                return new ProductionStub(156, bodyIDs, 2);
            }
        case 134:  // [declaration] ::= precision [precision-qualifier] [type-specifier-no-prec] ;
            {
                int bodyIDs[] = {62, 202, 203, 33};
                return new ProductionStub(156, bodyIDs, 4);
            }
        case 135:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } ;
            {
                int bodyIDs[] = {204, 142, 26, 205, 27, 33};
                return new ProductionStub(156, bodyIDs, 6);
            }
        case 136:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } IDENTIFIER ;
            {
                int bodyIDs[] = {204, 142, 26, 205, 27, 142, 33};
                return new ProductionStub(156, bodyIDs, 7);
            }
        case 137:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } IDENTIFIER [ ] ;
            {
                int bodyIDs[] = {204, 142, 26, 205, 27, 142, 24, 25, 33};
                return new ProductionStub(156, bodyIDs, 9);
            }
        case 138:  // [declaration] ::= [type-qualifier] IDENTIFIER { [struct-declaration-list] } IDENTIFIER [ [constant-expression] ] ;
            {
                int bodyIDs[] = {204, 142, 26, 205, 27, 142, 24, 200, 25, 33};
                return new ProductionStub(156, bodyIDs, 10);
            }
        case 139:  // [declaration] ::= [type-qualifier] ;
            {
                int bodyIDs[] = {204, 33};
                return new ProductionStub(156, bodyIDs, 2);
            }
        case 140:  // [function-prototype] ::= [function-header] )
            {
                int bodyIDs[] = {206, 23};
                return new ProductionStub(158, bodyIDs, 2);
            }
        case 141:  // [function-prototype] ::= [function-header-with-parameters] )
            {
                int bodyIDs[] = {207, 23};
                return new ProductionStub(158, bodyIDs, 2);
            }
        case 142:  // [function-header-with-parameters] ::= [function-header] [parameter-declaration]
            {
                int bodyIDs[] = {206, 208};
                return new ProductionStub(207, bodyIDs, 2);
            }
        case 143:  // [function-header-with-parameters] ::= [function-header-with-parameters] , [parameter-declaration]
            {
                int bodyIDs[] = {207, 31, 208};
                return new ProductionStub(207, bodyIDs, 3);
            }
        case 144:  // [function-header] ::= [fully-specified-type] IDENTIFIER (
            {
                int bodyIDs[] = {174, 142, 22};
                return new ProductionStub(206, bodyIDs, 3);
            }
        case 145:  // [parameter-declarator] ::= [type-specifier] IDENTIFIER
            {
                int bodyIDs[] = {186, 142};
                return new ProductionStub(209, bodyIDs, 2);
            }
        case 146:  // [parameter-declarator] ::= [type-specifier] IDENTIFIER [ [constant-expression] ]
            {
                int bodyIDs[] = {186, 142, 24, 200, 25};
                return new ProductionStub(209, bodyIDs, 5);
            }
        case 147:  // [parameter-declaration] ::= const [parameter-declarator]
            {
                int bodyIDs[] = {47, 209};
                return new ProductionStub(208, bodyIDs, 2);
            }
        case 148:  // [parameter-declaration] ::= const [parameter-qualifier] [parameter-declarator]
            {
                int bodyIDs[] = {47, 210, 209};
                return new ProductionStub(208, bodyIDs, 3);
            }
        case 149:  // [parameter-declaration] ::= [parameter-declarator]
            {
                int bodyIDs[] = {209};
                return new ProductionStub(208, bodyIDs, 1);
            }
        case 150:  // [parameter-declaration] ::= [parameter-qualifier] [parameter-declarator]
            {
                int bodyIDs[] = {210, 209};
                return new ProductionStub(208, bodyIDs, 2);
            }
        case 151:  // [parameter-declaration] ::= const [parameter-type-specifier]
            {
                int bodyIDs[] = {47, 211};
                return new ProductionStub(208, bodyIDs, 2);
            }
        case 152:  // [parameter-declaration] ::= const [parameter-qualifier] [parameter-type-specifier]
            {
                int bodyIDs[] = {47, 210, 211};
                return new ProductionStub(208, bodyIDs, 3);
            }
        case 153:  // [parameter-declaration] ::= [parameter-type-specifier]
            {
                int bodyIDs[] = {211};
                return new ProductionStub(208, bodyIDs, 1);
            }
        case 154:  // [parameter-declaration] ::= [parameter-qualifier] [parameter-type-specifier]
            {
                int bodyIDs[] = {210, 211};
                return new ProductionStub(208, bodyIDs, 2);
            }
        case 155:  // [parameter-qualifier] ::= in
            {
                int bodyIDs[] = {48};
                return new ProductionStub(210, bodyIDs, 1);
            }
        case 156:  // [parameter-qualifier] ::= out
            {
                int bodyIDs[] = {50};
                return new ProductionStub(210, bodyIDs, 1);
            }
        case 157:  // [parameter-qualifier] ::= inout
            {
                int bodyIDs[] = {49};
                return new ProductionStub(210, bodyIDs, 1);
            }
        case 158:  // [parameter-type-specifier] ::= [type-specifier]
            {
                int bodyIDs[] = {186};
                return new ProductionStub(211, bodyIDs, 1);
            }
        case 159:  // [init-declarator-list] ::= [single-declaration]
            {
                int bodyIDs[] = {212};
                return new ProductionStub(201, bodyIDs, 1);
            }
        case 160:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER
            {
                int bodyIDs[] = {201, 31, 142};
                return new ProductionStub(201, bodyIDs, 3);
            }
        case 161:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ ]
            {
                int bodyIDs[] = {201, 31, 142, 24, 25};
                return new ProductionStub(201, bodyIDs, 5);
            }
        case 162:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ [constant-expression] ]
            {
                int bodyIDs[] = {201, 31, 142, 24, 200, 25};
                return new ProductionStub(201, bodyIDs, 6);
            }
        case 163:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ ] = [initializer]
            {
                int bodyIDs[] = {201, 31, 142, 24, 25, 34, 175};
                return new ProductionStub(201, bodyIDs, 7);
            }
        case 164:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER [ [constant-expression] ] = [initializer]
            {
                int bodyIDs[] = {201, 31, 142, 24, 200, 25, 34, 175};
                return new ProductionStub(201, bodyIDs, 8);
            }
        case 165:  // [init-declarator-list] ::= [init-declarator-list] , IDENTIFIER = [initializer]
            {
                int bodyIDs[] = {201, 31, 142, 34, 175};
                return new ProductionStub(201, bodyIDs, 5);
            }
        case 166:  // [single-declaration] ::= [fully-specified-type]
            {
                int bodyIDs[] = {174};
                return new ProductionStub(212, bodyIDs, 1);
            }
        case 167:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER
            {
                int bodyIDs[] = {174, 142};
                return new ProductionStub(212, bodyIDs, 2);
            }
        case 168:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ ]
            {
                int bodyIDs[] = {174, 142, 24, 25};
                return new ProductionStub(212, bodyIDs, 4);
            }
        case 169:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ [constant-expression] ]
            {
                int bodyIDs[] = {174, 142, 24, 200, 25};
                return new ProductionStub(212, bodyIDs, 5);
            }
        case 170:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ ] = [initializer]
            {
                int bodyIDs[] = {174, 142, 24, 25, 34, 175};
                return new ProductionStub(212, bodyIDs, 6);
            }
        case 171:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER [ [constant-expression] ] = [initializer]
            {
                int bodyIDs[] = {174, 142, 24, 200, 25, 34, 175};
                return new ProductionStub(212, bodyIDs, 7);
            }
        case 172:  // [single-declaration] ::= [fully-specified-type] IDENTIFIER = [initializer]
            {
                int bodyIDs[] = {174, 142, 34, 175};
                return new ProductionStub(212, bodyIDs, 4);
            }
        case 173:  // [fully-specified-type] ::= [type-specifier]
            {
                int bodyIDs[] = {186};
                return new ProductionStub(174, bodyIDs, 1);
            }
        case 174:  // [fully-specified-type] ::= [type-qualifier] [type-specifier]
            {
                int bodyIDs[] = {204, 186};
                return new ProductionStub(174, bodyIDs, 2);
            }
        case 175:  // [interpolation-qualifier] ::= smooth
            {
                int bodyIDs[] = {57};
                return new ProductionStub(213, bodyIDs, 1);
            }
        case 176:  // [interpolation-qualifier] ::= flat
            {
                int bodyIDs[] = {56};
                return new ProductionStub(213, bodyIDs, 1);
            }
        case 177:  // [interpolation-qualifier] ::= noperspective
            {
                int bodyIDs[] = {58};
                return new ProductionStub(213, bodyIDs, 1);
            }
        case 178:  // [layout-qualifier] ::= layout ( [layout-qualifier-id-list] )
            {
                int bodyIDs[] = {55, 22, 214, 23};
                return new ProductionStub(215, bodyIDs, 4);
            }
        case 179:  // [layout-qualifier-id-list] ::= IDENTIFIER
            {
                int bodyIDs[] = {142};
                return new ProductionStub(214, bodyIDs, 1);
            }
        case 180:  // [layout-qualifier-id-list] ::= IDENTIFIER = INTCONSTANT
            {
                int bodyIDs[] = {142, 34, 145};
                return new ProductionStub(214, bodyIDs, 3);
            }
        case 181:  // [layout-qualifier-id-list] ::= [layout-qualifier-id-list] , IDENTIFIER
            {
                int bodyIDs[] = {214, 31, 142};
                return new ProductionStub(214, bodyIDs, 3);
            }
        case 182:  // [layout-qualifier-id-list] ::= [layout-qualifier-id-list] , IDENTIFIER = INTCONSTANT
            {
                int bodyIDs[] = {214, 31, 142, 34, 145};
                return new ProductionStub(214, bodyIDs, 5);
            }
        case 183:  // [type-qualifier] ::= [storage-qualifier]
            {
                int bodyIDs[] = {216};
                return new ProductionStub(204, bodyIDs, 1);
            }
        case 184:  // [type-qualifier] ::= [layout-qualifier]
            {
                int bodyIDs[] = {215};
                return new ProductionStub(204, bodyIDs, 1);
            }
        case 185:  // [type-qualifier] ::= [layout-qualifier] [storage-qualifier]
            {
                int bodyIDs[] = {215, 216};
                return new ProductionStub(204, bodyIDs, 2);
            }
        case 186:  // [type-qualifier] ::= [interpolation-qualifier] [storage-qualifier]
            {
                int bodyIDs[] = {213, 216};
                return new ProductionStub(204, bodyIDs, 2);
            }
        case 187:  // [type-qualifier] ::= [interpolation-qualifier]
            {
                int bodyIDs[] = {213};
                return new ProductionStub(204, bodyIDs, 1);
            }
        case 188:  // [type-qualifier] ::= invariant [storage-qualifier]
            {
                int bodyIDs[] = {54, 216};
                return new ProductionStub(204, bodyIDs, 2);
            }
        case 189:  // [type-qualifier] ::= invariant [interpolation-qualifier] [storage-qualifier]
            {
                int bodyIDs[] = {54, 213, 216};
                return new ProductionStub(204, bodyIDs, 3);
            }
        case 190:  // [type-qualifier] ::= invariant
            {
                int bodyIDs[] = {54};
                return new ProductionStub(204, bodyIDs, 1);
            }
        case 191:  // [storage-qualifier] ::= const
            {
                int bodyIDs[] = {47};
                return new ProductionStub(216, bodyIDs, 1);
            }
        case 192:  // [storage-qualifier] ::= attribute
            {
                int bodyIDs[] = {46};
                return new ProductionStub(216, bodyIDs, 1);
            }
        case 193:  // [storage-qualifier] ::= varying
            {
                int bodyIDs[] = {52};
                return new ProductionStub(216, bodyIDs, 1);
            }
        case 194:  // [storage-qualifier] ::= centroid varying
            {
                int bodyIDs[] = {53, 52};
                return new ProductionStub(216, bodyIDs, 2);
            }
        case 195:  // [storage-qualifier] ::= in
            {
                int bodyIDs[] = {48};
                return new ProductionStub(216, bodyIDs, 1);
            }
        case 196:  // [storage-qualifier] ::= out
            {
                int bodyIDs[] = {50};
                return new ProductionStub(216, bodyIDs, 1);
            }
        case 197:  // [storage-qualifier] ::= centroid in
            {
                int bodyIDs[] = {53, 48};
                return new ProductionStub(216, bodyIDs, 2);
            }
        case 198:  // [storage-qualifier] ::= centroid out
            {
                int bodyIDs[] = {53, 50};
                return new ProductionStub(216, bodyIDs, 2);
            }
        case 199:  // [storage-qualifier] ::= uniform
            {
                int bodyIDs[] = {51};
                return new ProductionStub(216, bodyIDs, 1);
            }
        case 200:  // [type-specifier] ::= [type-specifier-no-prec]
            {
                int bodyIDs[] = {203};
                return new ProductionStub(186, bodyIDs, 1);
            }
        case 201:  // [type-specifier] ::= [precision-qualifier] [type-specifier-no-prec]
            {
                int bodyIDs[] = {202, 203};
                return new ProductionStub(186, bodyIDs, 2);
            }
        case 202:  // [precision-qualifier] ::= highp
            {
                int bodyIDs[] = {60};
                return new ProductionStub(202, bodyIDs, 1);
            }
        case 203:  // [precision-qualifier] ::= mediump
            {
                int bodyIDs[] = {61};
                return new ProductionStub(202, bodyIDs, 1);
            }
        case 204:  // [precision-qualifier] ::= lowp
            {
                int bodyIDs[] = {59};
                return new ProductionStub(202, bodyIDs, 1);
            }
        case 205:  // [type-specifier-no-prec] ::= [type-specifier-nonarray]
            {
                int bodyIDs[] = {217};
                return new ProductionStub(203, bodyIDs, 1);
            }
        case 206:  // [type-specifier-no-prec] ::= [type-specifier-nonarray] [ ]
            {
                int bodyIDs[] = {217, 24, 25};
                return new ProductionStub(203, bodyIDs, 3);
            }
        case 207:  // [type-specifier-no-prec] ::= [type-specifier-nonarray] [ [constant-expression] ]
            {
                int bodyIDs[] = {217, 24, 200, 25};
                return new ProductionStub(203, bodyIDs, 4);
            }
        case 208:  // [type-specifier-nonarray] ::= void
            {
                int bodyIDs[] = {78};
                return new ProductionStub(217, bodyIDs, 1);
            }
        case 209:  // [type-specifier-nonarray] ::= NATIVE-TYPE
            {
                int bodyIDs[] = {148};
                return new ProductionStub(217, bodyIDs, 1);
            }
        case 210:  // [type-specifier-nonarray] ::= [struct-specifier]
            {
                int bodyIDs[] = {218};
                return new ProductionStub(217, bodyIDs, 1);
            }
        case 211:  // [type-specifier-nonarray] ::= TYPE-NAME
            {
                int bodyIDs[] = {143};
                return new ProductionStub(217, bodyIDs, 1);
            }
        case 212:  // [struct-specifier] ::= struct IDENTIFIER { [struct-declaration-list] }
            {
                int bodyIDs[] = {75, 142, 26, 205, 27};
                return new ProductionStub(218, bodyIDs, 5);
            }
        case 213:  // [struct-specifier] ::= struct { [struct-declaration-list] }
            {
                int bodyIDs[] = {75, 26, 205, 27};
                return new ProductionStub(218, bodyIDs, 4);
            }
        case 214:  // [struct-declaration-list] ::= [struct-declaration]
            {
                int bodyIDs[] = {219};
                return new ProductionStub(205, bodyIDs, 1);
            }
        case 215:  // [struct-declaration-list] ::= [struct-declaration-list] [struct-declaration]
            {
                int bodyIDs[] = {205, 219};
                return new ProductionStub(205, bodyIDs, 2);
            }
        case 216:  // [struct-declaration] ::= [type-specifier] [struct-declarator-list] ;
            {
                int bodyIDs[] = {186, 220, 33};
                return new ProductionStub(219, bodyIDs, 3);
            }
        case 217:  // [struct-declaration] ::= [type-qualifier] [type-specifier] [struct-declarator-list] ;
            {
                int bodyIDs[] = {204, 186, 220, 33};
                return new ProductionStub(219, bodyIDs, 4);
            }
        case 218:  // [struct-declarator-list] ::= [struct-declarator]
            {
                int bodyIDs[] = {221};
                return new ProductionStub(220, bodyIDs, 1);
            }
        case 219:  // [struct-declarator-list] ::= [struct-declarator-list] , [struct-declarator]
            {
                int bodyIDs[] = {220, 31, 221};
                return new ProductionStub(220, bodyIDs, 3);
            }
        case 220:  // [struct-declarator] ::= IDENTIFIER
            {
                int bodyIDs[] = {142};
                return new ProductionStub(221, bodyIDs, 1);
            }
        case 221:  // [struct-declarator] ::= IDENTIFIER [ ]
            {
                int bodyIDs[] = {142, 24, 25};
                return new ProductionStub(221, bodyIDs, 3);
            }
        case 222:  // [struct-declarator] ::= IDENTIFIER [ [constant-expression] ]
            {
                int bodyIDs[] = {142, 24, 200, 25};
                return new ProductionStub(221, bodyIDs, 4);
            }
        case 223:  // [initializer] ::= [assignment-expression]
            {
                int bodyIDs[] = {185};
                return new ProductionStub(175, bodyIDs, 1);
            }
    } // switch (productionID

    return 0;
}   // findProduction()

std::string GLSLParser::symbolID2String(const int symbolID) const {
    switch (symbolID) {
        case 0: return "[$START$]";
        case 1: return "<<";
        case 2: return ">>";
        case 3: return "++";
        case 4: return "--";
        case 5: return "<=";
        case 6: return ">=";
        case 7: return "==";
        case 8: return "!=";
        case 9: return "&&";
        case 10: return "||";
        case 11: return "^^";
        case 12: return "<<=";
        case 13: return ">>=";
        case 14: return "+=";
        case 15: return "/=";
        case 16: return "*=";
        case 17: return "-=";
        case 18: return "%=";
        case 19: return "&=";
        case 20: return "|=";
        case 21: return "^=";
        case 22: return "(";
        case 23: return ")";
        case 24: return "[";
        case 25: return "]";
        case 26: return "{";
        case 27: return "}";
        case 28: return "<";
        case 29: return ">";
        case 30: return ".";
        case 31: return ",";
        case 32: return ":";
        case 33: return ";";
        case 34: return "=";
        case 35: return "!";
        case 36: return "-";
        case 37: return "~";
        case 38: return "+";
        case 39: return "*";
        case 40: return "/";
        case 41: return "%";
        case 42: return "|";
        case 43: return "^";
        case 44: return "&";
        case 45: return "?";
        case 46: return "attribute";
        case 47: return "const";
        case 48: return "in";
        case 49: return "inout";
        case 50: return "out";
        case 51: return "uniform";
        case 52: return "varying";
        case 53: return "centroid";
        case 54: return "invariant";
        case 55: return "layout";
        case 56: return "flat";
        case 57: return "smooth";
        case 58: return "noperspective";
        case 59: return "lowp";
        case 60: return "highp";
        case 61: return "mediump";
        case 62: return "precision";
        case 63: return "do";
        case 64: return "else";
        case 65: return "for";
        case 66: return "if";
        case 67: return "while";
        case 68: return "case";
        case 69: return "default";
        case 70: return "switch";
        case 71: return "break";
        case 72: return "continue";
        case 73: return "discard";
        case 74: return "return";
        case 75: return "struct";
        case 76: return "false";
        case 77: return "true";
        case 78: return "void";
        // added manually
        //
        case GLSLTerminals::ID_BOOL: return "bool";
        case GLSLTerminals::ID_FLOAT: return "float";
        case GLSLTerminals::ID_INT: return "int";
        case GLSLTerminals::ID_UINT: return "uint";
        case GLSLTerminals::ID_BVEC2: return "bvec2";
        case GLSLTerminals::ID_BVEC3: return "bvec3";
        case GLSLTerminals::ID_BVEC4: return "bvec4";
        case GLSLTerminals::ID_IVEC2: return "ivec2";
        case GLSLTerminals::ID_IVEC3: return "ivec3";
        case GLSLTerminals::ID_IVEC4: return "ivec4";
        case GLSLTerminals::ID_UVEC2: return "uvec2";
        case GLSLTerminals::ID_UVEC3: return "uvec3";
        case GLSLTerminals::ID_UVEC4: return "uvec4";
        case GLSLTerminals::ID_VEC2: return "vec2";
        case GLSLTerminals::ID_VEC3: return "vec3";
        case GLSLTerminals::ID_VEC4: return "vec4";
        case GLSLTerminals::ID_MAT2: return "mat2";
        case GLSLTerminals::ID_MAT3: return "mat3";
        case GLSLTerminals::ID_MAT4: return "mat4";
        case GLSLTerminals::ID_MAT2X2: return "mat2x2";
        case GLSLTerminals::ID_MAT2X3: return "mat2x3";
        case GLSLTerminals::ID_MAT2X4: return "mat2x4";
        case GLSLTerminals::ID_MAT3X2: return "mat3x2";
        case GLSLTerminals::ID_MAT3X3: return "mat3x3";
        case GLSLTerminals::ID_MAT3X4: return "mat3x4";
        case GLSLTerminals::ID_MAT4X2: return "mat4x2";
        case GLSLTerminals::ID_MAT4X3: return "mat4x3";
        case GLSLTerminals::ID_MAT4X4: return "mat4x4";
        case GLSLTerminals::ID_SAMPLER1D: return "sampler1D";
        case GLSLTerminals::ID_SAMPLER2D: return "sampler2D";
        case GLSLTerminals::ID_SAMPLER3D: return "sampler3D";
        case GLSLTerminals::ID_SAMPLERCUBE: return "samplerCube";
        case GLSLTerminals::ID_SAMPLER1DSHADOW: return "sampler1DShadow";
        case GLSLTerminals::ID_SAMPLER2DSHADOW: return "sampler2DShadow";
        case GLSLTerminals::ID_SAMPLER1DARRAY: return "sampler1DArray";
        case GLSLTerminals::ID_SAMPLER2DARRAY: return "sampler2DArray";
        case GLSLTerminals::ID_SAMPLER1DARRAYSHADOW: return "sampler1DArrayShadow";
        case GLSLTerminals::ID_SAMPLER2DARRAYSHADOW: return "sampler2DArrayShadow";
        case GLSLTerminals::ID_ISAMPLER1D: return "isampler1D";
        case GLSLTerminals::ID_ISAMPLER2D: return "isampler2D";
        case GLSLTerminals::ID_ISAMPLER3D: return "isampler3D";
        case GLSLTerminals::ID_ISAMPLERCUBE: return "isamplerCube";
        case GLSLTerminals::ID_ISAMPLER1DARRAY: return "isampler1DArray";
        case GLSLTerminals::ID_ISAMPLER2DARRAY: return "isampler2DArray";
        case GLSLTerminals::ID_USAMPLER1D: return "usampler1D";
        case GLSLTerminals::ID_USAMPLER2D: return "usampler2D";
        case GLSLTerminals::ID_USAMPLER3D: return "usampler3D";
        case GLSLTerminals::ID_USAMPLERCUBE: return "usamplerCube";
        case GLSLTerminals::ID_USAMPLER1DARRAY: return "usampler1DArray";
        case GLSLTerminals::ID_USAMPLER2DARRAY: return "usampler2DArray";
        case GLSLTerminals::ID_SAMPLER2DRECT: return "sampler2DRect";
        case GLSLTerminals::ID_SAMPLER2DRECTSHADOW: return "sampler2DRectShadow";
        case GLSLTerminals::ID_ISAMPLER2DRECT: return "isampler2DRect";
        case GLSLTerminals::ID_USAMPLER2DRECT: return "usampler2DRect";
        case GLSLTerminals::ID_SAMPLERBUFFER: return "samplerBuffer";
        case GLSLTerminals::ID_ISAMPLERBUFFER: return "isamplerBuffer";
        case GLSLTerminals::ID_USAMPLEBUFFER: return "usamplerBuffer";
        case GLSLTerminals::ID_SAMPLER2DMS: return "sampler2DMS";
        case GLSLTerminals::ID_ISAMPLER2DMS: return "isampler2DMS";
        case GLSLTerminals::ID_USAMPLER2DMS: return "usampler2DMS";
        case GLSLTerminals::ID_SAMPLER2DMSARRY: return "sampler2DMSArray";
        case GLSLTerminals::ID_ISAMPLER2DMSARRAY: return "isampler2DMSArray";
        case GLSLTerminals::ID_USAMPLER2DMSARRAY: return "usampler2DMSArray";
        //
        // added manually
        case 142: return "IDENTIFIER";
        case 143: return "TYPE-NAME";
        case 144: return "FLOATCONSTANT";
        case 145: return "INTCONSTANT";
        case 146: return "UINTCONSTANT";
        case 147: return "FIELD-SELECTION";
        case 148: return "NATIVE-TYPE";
        case 149: return "LEADING-ANNOTATION";
        case 150: return "TRAILING-ANNOTATION";
        case 151: return "RESERVED-KEYWORD";
        case 152: return "[program]";
        case 153: return "[translation-unit]";
        case 154: return "[external-declaration]";
        case 155: return "[function-definition]";
        case 156: return "[declaration]";
        case 157: return "[enhanced-declaration]";
        case 158: return "[function-prototype]";
        case 159: return "[compound-statement-no-new-scope]";
        case 160: return "[compound-statement]";
        case 161: return "[statement]";
        case 162: return "[simple-statement]";
        case 163: return "[declaration-statement]";
        case 164: return "[expression-statement]";
        case 165: return "[selection-statement]";
        case 166: return "[switch-statement]";
        case 167: return "[case-label]";
        case 168: return "[iteration-statement]";
        case 169: return "[jump-statement]";
        case 170: return "[statement-list]";
        case 171: return "[statement-no-new-scope]";
        case 172: return "[expression]";
        case 173: return "[condition]";
        case 174: return "[fully-specified-type]";
        case 175: return "[initializer]";
        case 176: return "[for-init-statement]";
        case 177: return "[for-rest-statement]";
        case 178: return "[primary-expression]";
        case 179: return "[postfix-expression]";
        case 180: return "[integer-expression]";
        case 181: return "[function-call-or-method]";
        case 182: return "[function-call-generic]";
        case 183: return "[function-call-header]";
        case 184: return "[function-call-header-with-parameters]";
        case 185: return "[assignment-expression]";
        case 186: return "[type-specifier]";
        case 187: return "[unary-expression]";
        case 188: return "[multiplicative-expression]";
        case 189: return "[additive-expression]";
        case 190: return "[shift-expression]";
        case 191: return "[relational-expression]";
        case 192: return "[equality-expression]";
        case 193: return "[and-expression]";
        case 194: return "[exclusive-or-expression]";
        case 195: return "[inclusive-or-expression]";
        case 196: return "[logical-and-expression]";
        case 197: return "[logical-xor-expression]";
        case 198: return "[logical-or-expression]";
        case 199: return "[conditional-expression]";
        case 200: return "[constant-expression]";
        case 201: return "[init-declarator-list]";
        case 202: return "[precision-qualifier]";
        case 203: return "[type-specifier-no-prec]";
        case 204: return "[type-qualifier]";
        case 205: return "[struct-declaration-list]";
        case 206: return "[function-header]";
        case 207: return "[function-header-with-parameters]";
        case 208: return "[parameter-declaration]";
        case 209: return "[parameter-declarator]";
        case 210: return "[parameter-qualifier]";
        case 211: return "[parameter-type-specifier]";
        case 212: return "[single-declaration]";
        case 213: return "[invariant-qualifier]";
        case 214: return "[interpolation-qualifier]";
        case 215: return "[layout-qualifier-id-list]";
        case 216: return "[layout-qualifier]";
        case 217: return "[storage-qualifier]";
        case 218: return "[type-specifier-nonarray]";
        case 219: return "[struct-specifier]";
        case 220: return "[struct-declaration]";
        case 221: return "[struct-declarator-list]";
        case 222: return "[struct-declarator]";
        case 223: return "$END$";
    }  // switch (symbolID
    return "<missing>";
}

}   // namespace glslparser

}   // namespace voreen
