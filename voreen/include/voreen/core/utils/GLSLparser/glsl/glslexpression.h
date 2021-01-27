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

#ifndef VRN_GLSLEXPRESSION_H
#define VRN_GLSLEXPRESSION_H

#include "voreen/core/utils/GLSLparser/glsl/glslparsetreenode.h"
#include "voreen/core/utils/GLSLparser/glsl/glslvalue.h"

namespace voreen {

namespace glslparser {

/**
 * NOTE: do not confuse this class GLSLExpression with the item [expression]
 * in GLSL grammar! This class is only for performing semantic actions on (i.e.
 * evaluating) the parse tree when it is traversed by an appropriate visitor.
 * The parse tree is constructed by the parsers' <code>expandParseTree()</code>
 * method. That method "packs" the grammar (i.e. its production) into a parse
 * tree.
 */
class GLSLExpression : public GLSLNode {
public:
    GLSLExpression()
        : GLSLNode(GLSLTerminals::ID_UNKNOWN),
        token_(0)
    {
    }

    GLSLExpression(const Token& token)
        : GLSLNode(token.getTokenID()),
        token_(token.getCopy())
    {
    }

    virtual ~GLSLExpression() {
        delete token_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_EXPRESSION; }

    Token* getToken() const { return token_; }

protected:
    Token* const token_;
};

// ============================================================================

class GLSLPrimaryExpression : public GLSLExpression {
public:
    GLSLPrimaryExpression(GLSLExpression* const expr)
        : GLSLExpression(Token(GLSLTerminals::ID_LPAREN)),
        expr_(expr)
    {
    }

    GLSLPrimaryExpression(const Token& token)
        : GLSLExpression(token),
        expr_(0)
    {
    }

    virtual ~GLSLPrimaryExpression() {
        delete expr_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_PRIMARY_EXPRESSION; }

    /**
     * @return  Returns only non-null value if the primary expression is
     *          a parenthesis expression, i.e. the ID of getToken() is
     *          ID_LPAREN
     */
    GLSLExpression* getExpression() const { return expr_; }

protected:
    GLSLExpression* const expr_;
};

// ============================================================================

class GLSLPostfixExpression : public GLSLExpression {
public:
    GLSLPostfixExpression(Token* const token, GLSLPostfixExpression* const operand)
        : GLSLExpression(*token),
        operand_(operand),
        field_(0),
        intExpression_(0)
    {
    }

    GLSLPostfixExpression(GLSLExpression* const expr)
        : GLSLExpression(),
        operand_(expr),
        field_(0),
        intExpression_(0)
    {
    }

    GLSLPostfixExpression(GLSLPostfixExpression* const operand, IdentifierToken* const field)
        : GLSLExpression(Token(GLSLTerminals::ID_DOT)),
        operand_(operand),
        field_(dynamic_cast<IdentifierToken* const>(field->getCopy())),
        intExpression_(0)
    {
    }

    GLSLPostfixExpression(GLSLPostfixExpression* const arr, GLSLExpression* const intExpr)
        : GLSLExpression(Token(GLSLTerminals::ID_LBRACKET)),
        operand_(arr),
        field_(0),
        intExpression_(intExpr)
    {
    }

    virtual ~GLSLPostfixExpression() {
        delete operand_;
        delete field_;
        delete intExpression_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_POSTFIX_EXPRESSION; }

    /**
     * Returns only a value != NULL, if the ID of the token returned by getToken()
     * returns ID_DOT.
     */
    IdentifierToken* getField() { return field_; }

    /**
     * Returns only a value != NULL, if the ID of the token returned by getToken()
     * returns ID_LBRACKET.
     */
    GLSLExpression* getIntExpression() { return intExpression_; }

    GLSLExpression* getOperand() { return operand_; }

protected:
    GLSLExpression* const operand_;
    IdentifierToken* const field_;
    GLSLExpression* const intExpression_;
};

// ============================================================================

class GLSLUnaryExpression : public GLSLExpression {
public:
    GLSLUnaryExpression(GLSLPostfixExpression* const expr)
        : GLSLExpression(),
        expression_(expr)
    {
    }

    GLSLUnaryExpression(Token* const token, GLSLUnaryExpression* const expr)
        : GLSLExpression(*token),
        expression_(expr)
    {
    }

    virtual ~GLSLUnaryExpression() {
        delete expression_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_UNARY_EXPRESSION; }

    GLSLExpression* getExpression() { return expression_; }

protected:
    GLSLExpression* const expression_;
};

// ============================================================================

class GLSLBinaryExpression : public GLSLExpression {
public:
    GLSLBinaryExpression(Token* const token, GLSLExpression* const lhs, GLSLExpression* const rhs)
        : GLSLExpression(*token),
        lhs_(lhs),
        rhs_(rhs)
    {
    }

    virtual ~GLSLBinaryExpression() {
        delete lhs_;
        delete rhs_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_BINARY_EXPRESSION; }

    GLSLExpression* getLHS() { return lhs_; }
    GLSLExpression* getRHS() { return rhs_; }

protected:
    GLSLExpression* const lhs_;
    GLSLExpression* const rhs_;
};

// ============================================================================

class GLSLAssignmentExpression : public GLSLExpression {
public:
    GLSLAssignmentExpression(Token* const token, GLSLUnaryExpression* const lvalue,
        GLSLExpression* const rvalue)
        : GLSLExpression(*token),
        lvalue_(lvalue),
        rvalue_(rvalue)
    {
    }

    virtual ~GLSLAssignmentExpression() {
        delete lvalue_;
        delete rvalue_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_ASSIGNMENT_EXPRESSION; }

    GLSLUnaryExpression* getLValue() const { return lvalue_; }
    GLSLExpression* getRValue() const { return rvalue_; }

protected:
    GLSLUnaryExpression* const lvalue_;
    GLSLExpression* const rvalue_;
};

// ============================================================================

class GLSLConditionalExpression : public GLSLExpression {
public:
    GLSLConditionalExpression(GLSLExpression* const condExpr,
        GLSLExpression* const trueExpr, GLSLExpression* const falseExpr)
        : GLSLExpression(Token(GLSLTerminals::ID_QUESTION)),
        condition_(condExpr),
        trueExpression_(trueExpr),
        falseExpression_(falseExpr)
    {
    }

    virtual ~GLSLConditionalExpression() {
        delete condition_;
        delete trueExpression_;
        delete falseExpression_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_CONDITIONAL_EXPRESSION; }

protected:
    GLSLExpression* const condition_;
    GLSLExpression* const trueExpression_;
    GLSLExpression* const falseExpression_;
};

// ============================================================================

class GLSLExpressionList : public GLSLExpression {
public:
    GLSLExpressionList(GLSLExpression* const expr)
        : GLSLExpression(Token(GLSLTerminals::ID_COMMA))
    {
        addExpression(expr);
    }

    virtual ~GLSLExpressionList() {
        for (size_t i = 0; i < expressions_.size(); ++i)
            delete expressions_[i];
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_EXPRESSION_LIST; }

    void addExpression(GLSLExpression* const expr) {
        if (expr != 0)
            expressions_.push_back(expr);
    }

    const std::vector<GLSLExpression*>& getExpressions() { return expressions_; }

protected:
    std::vector<GLSLExpression*> expressions_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
