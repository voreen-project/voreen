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

#ifndef VRN_PPEXPRESSION_H
#define VRN_PPEXPRESSION_H

#include "voreen/core/utils/GLSLparser/preprocessor/ppparsetreenode.h"

namespace voreen {

namespace glslparser {

class Expression : public ParseTreeNode {
public:
    Expression(Token* const token)
        : ParseTreeNode(token->getTokenID()),
        token_(token->getCopy())
    {
    }

    virtual ~Expression() {
        delete token_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_EXPRESSION; }
    Token* getToken() const { return token_; }

    virtual void tokenize(TokenList* const list) const {
        list->addToken(token_->getCopy());
    }

private:
    Token* const token_;
};

// ============================================================================

class ExpressionList : public ParseTreeNode {
public:
    ExpressionList(Expression* const exp)
        : ParseTreeNode(PreprocessorTerminals::ID_UNKNOWN)
    {
        addExpression(exp);
    }

    void addExpression(Expression* const exp) {
        if (exp != 0)
            expressions_.push_back(exp);
    }

    virtual ~ExpressionList() {
        for (size_t i = 0; i < expressions_.size(); ++i)
            delete expressions_[i];
    }

    Expression* operator[](const size_t index) {
        if (index < expressions_.size())
            return expressions_[index];

        throw std::runtime_error("ExpressionList::operator[]: index out of bounds!");
    }

    const Expression* operator[](const size_t index) const {
        if (index < expressions_.size())
            return expressions_[index];

        throw std::runtime_error("ExpressionList::operator[]: index out of bounds!");
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_EXPRESSIONLIST; }
    size_t size() const { return expressions_.size(); }

protected:
    std::vector<Expression*> expressions_;
};

// ============================================================================

class ParenthesisExpression : public Expression {
public:
    ParenthesisExpression(Token* const lparen, Expression* const interior)
        : Expression(lparen),
        interior_(interior)
    {
    }

    virtual ~ParenthesisExpression() {
        delete interior_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_PARENTHESIS; }

    Expression* getInterior() const { return interior_; }

    virtual void tokenize(TokenList* const list) const {
        list->addToken(new Token(PreprocessorTerminals::ID_LPAREN));
        interior_->tokenize(list);
        list->addToken(new Token(PreprocessorTerminals::ID_RPAREN));
    }

protected:
    Expression* const interior_;
};

// ============================================================================

class ArithmeticExpression : public Expression {
public:
    ArithmeticExpression(Token* const token)
        : Expression(token)
    {
    }

    virtual ~ArithmeticExpression() {}

    int getNodeType() const { return PreprocessorNodeTypes::NODE_ARITHMETIC; }
};

// ============================================================================

class IntConstant : public ArithmeticExpression {
public:
    IntConstant(ConstantToken* const token)
        : ArithmeticExpression(token),
        value_(token->convert<int>())
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_INT_CONSTANT; }

    int getValue() const { return value_; }

protected:
    int value_;
};

// ============================================================================

class UnaryExpression : public ArithmeticExpression {
public:
    UnaryExpression(Token* const token, Expression* const expr)
        : ArithmeticExpression(token),
        expr_(expr)
    {
    }

    virtual ~UnaryExpression() {
        delete expr_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_UNARY_EXPRESSION; }

    Expression* getExpression() const { return expr_; }

    virtual void tokenize(TokenList* const list) const {
        Expression::tokenize(list);
        expr_->tokenize(list);
    }

protected:
    Expression* const expr_;
};

// ============================================================================

class BinaryExpression : public ArithmeticExpression {
public:
    BinaryExpression(Token* const token, Expression* const lhs, Expression* const rhs)
        : ArithmeticExpression(token),
        lhs_(lhs),
        rhs_(rhs)
    {
    }

    virtual ~BinaryExpression() {
        delete lhs_;
        delete rhs_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_BINARY_EXPRESSION; }

    Expression* getLeft() const { return lhs_; }
    Expression* getRight() const { return rhs_; }

    virtual void tokenize(TokenList* const list) const {
        lhs_->tokenize(list);
        Expression::tokenize(list);
        rhs_->tokenize(list);
    }

protected:
    Expression* const lhs_;
    Expression* const rhs_;
};

// ============================================================================

class LogicalExpression : public Expression {
public:
    LogicalExpression(Token* const token) : Expression(token) {
    }

    //virtual ~LogicalExpression() {}
    //int getNodeType() const { return PreprocessorNodeTypes::NODE_LOGICAL_EXPRESSION; }
};

// ============================================================================

class LogicalBinaryExpression : public LogicalExpression {
public:
    LogicalBinaryExpression(Token* const token, Expression* const lhs, Expression* const rhs)
        : LogicalExpression(token),
        lhs_(lhs),
        rhs_(rhs)
    {
    }

    virtual ~LogicalBinaryExpression() {
        delete lhs_;
        delete rhs_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_BINARY_LOGICAL; }

    Expression* getLeft() const { return lhs_; }
    Expression* getRight() const { return rhs_; }

    virtual void tokenize(TokenList* const list) const {
        lhs_->tokenize(list);
        Expression::tokenize(list);
        rhs_->tokenize(list);
    }

protected:
    Expression* const lhs_;
    Expression* const rhs_;
};

// ============================================================================

class DefinedOperator : public LogicalExpression {
public:
    DefinedOperator(IdentifierToken* const token)
        : LogicalExpression(token),
        identifier_(token->getValue())
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_DEFINED_OPERATOR; }

    const std::string& getIdentifier() const { return identifier_; }

protected:
    std::string identifier_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
