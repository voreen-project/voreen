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

#ifndef VRN_GLSLSTATEMENT_H
#define VRN_GLSLSTATEMENT_H

#include "voreen/core/utils/GLSLparser/glsl/glslcondition.h"
#include "voreen/core/utils/GLSLparser/glsl/glsldeclaration.h"
#include <vector>

namespace voreen {

namespace glslparser {

class GLSLStatement : public GLSLNode
{
public:
    GLSLStatement(const int symbolID, const bool newScope)
        : GLSLNode(symbolID),
        newScope_(newScope)
    {
    }

    virtual ~GLSLStatement() {}

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_STATEMENT ; }

    bool ownsNewScope() const { return newScope_; }
    void setNewScope(const bool newScope) { newScope_ = newScope; }

private:
    bool newScope_;
};  // class GLSLStatement

// ============================================================================

class GLSLStatementList : public GLSLNode
{
public:
    GLSLStatementList(GLSLStatement* const statement)
        : GLSLNode(GLSLTerminals::ID_COMMA)
    {
        addStatement(statement);
    }

    virtual ~GLSLStatementList() {
        for (size_t i = 0; i < statements_.size(); ++i)
            delete statements_[i];
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_STATEMENT_LIST ; }

    void addStatement(GLSLStatement* const statement) {
        if (statement != 0)
            statements_.push_back(statement);
    }

    const std::vector<GLSLStatement*>& getStatements() const { return statements_; }

protected:
    std::vector<GLSLStatement*> statements_;

};  // class GLSLStatementList

// ============================================================================

class GLSLSimpleStatement : public GLSLStatement
{
public:
    GLSLSimpleStatement(const int symbolID)
        : GLSLStatement(symbolID, false)
    {
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_SIMPLE_STATEMENT ; }

};  // class GLSLSimpleStatement

// ============================================================================

class GLSLCaseLabel : public GLSLSimpleStatement
{
public:
    GLSLCaseLabel()
        : GLSLSimpleStatement(GLSLTerminals::ID_DEFAULT),
        cond_()
    {
    }

    GLSLCaseLabel(GLSLExpression* const cond)
        : GLSLSimpleStatement(GLSLTerminals::ID_CASE),
        cond_(cond)
    {
    }

    virtual ~GLSLCaseLabel() {
        delete cond_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_CASE_LABEL; }

    GLSLExpression* getCondition() { return cond_; }

protected:
    GLSLExpression* const cond_;
};  // class GLSLCaseLabel

// ============================================================================

class GLSLDeclarationStatement : public GLSLSimpleStatement
{
public:
    GLSLDeclarationStatement(GLSLDeclaration* const decl)
        : GLSLSimpleStatement(GLSLTerminals::ID_UNKNOWN),
        decl_(decl)
    {
    }

    virtual ~GLSLDeclarationStatement() {
        delete decl_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_DECLARATION_STATEMENT; }

    GLSLDeclaration* getDeclaration() { return decl_; }

protected:
    GLSLDeclaration* const decl_;
};  // class GLSLDeclarationStatement

// ============================================================================

class GLSLDoWhileStatement : public GLSLSimpleStatement {
public:
    GLSLDoWhileStatement(GLSLExpression* const cond, GLSLStatement* const stmt)
        : GLSLSimpleStatement(GLSLTerminals::ID_DO),
        stmt_(stmt),
        cond_(cond)
    {
    }

    virtual ~GLSLDoWhileStatement() {
        delete stmt_;
        delete cond_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_DO_WHILE_STATEMENT; }

    GLSLStatement* getBody() { return stmt_; }
    GLSLExpression* getCondition() { return cond_; }

protected:
    GLSLStatement* const stmt_;
    GLSLExpression* const cond_;
};  // class GLSLDoWhileStatement

// ============================================================================

class GLSLExpressionStatement : public GLSLSimpleStatement
{
public:
    GLSLExpressionStatement(GLSLExpression* const expr)
        : GLSLSimpleStatement(GLSLTerminals::ID_SEMICOLON),
        expr_(expr)
    {
    }

    virtual ~GLSLExpressionStatement() {
        delete expr_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_EXPRESSION_STATEMENT; }

    GLSLExpression* getExpression() { return expr_; }

protected:
    GLSLExpression* const expr_;
};  // class GLSLExpressionStatement

// ============================================================================

class GLSLForStatement : public GLSLSimpleStatement
{
public:
    GLSLForStatement()
        : GLSLSimpleStatement(GLSLTerminals::ID_FOR),
        init_(0),
        cond_(0),
        iter_(0),
        stmt_(0)
    {
    }

    GLSLForStatement(GLSLCondition* const cond, GLSLExpression* const expr)
        : GLSLSimpleStatement(GLSLTerminals::ID_FOR),
        init_(0),
        cond_(cond),
        iter_(expr),
        stmt_(0)
    {
    }

    virtual ~GLSLForStatement() {
        delete init_;
        delete cond_;
        delete iter_;
        delete stmt_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_FOR_STATEMENT; }

    GLSLSimpleStatement* getInit() { return init_; }
    GLSLCondition* getCondition() { return cond_; }
    GLSLExpression* getIterationExpr() { return iter_; }
    GLSLStatement* getBody() { return stmt_; }

    void setBody(GLSLStatement* const body) { stmt_ = body; }
    void setInitStatement(GLSLSimpleStatement* const init) {init_ = init; }

protected:
    GLSLSimpleStatement* init_;  /** Must be expression statement or declaration statement! */
    GLSLCondition* const cond_;
    GLSLExpression* const iter_;
    GLSLStatement* stmt_;
};  // class GLSLForStatement

// ============================================================================

class GLSLJumpStatement : public GLSLSimpleStatement
{
public:
    GLSLJumpStatement(Token* const token)
        : GLSLSimpleStatement(token->getTokenID()),
        expr_(0)
    {
    }

    GLSLJumpStatement(GLSLExpression* const expr)
        : GLSLSimpleStatement(GLSLTerminals::ID_RETURN),
        expr_(expr)
    {
    }

    virtual ~GLSLJumpStatement() {
        delete expr_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_JUMP_STATEMENT; }

    /**
     * @return  Result is only non-null if the token return by getToken() is ID_RETURN
     */
    GLSLExpression* getExpression() { return expr_; }

protected:
    GLSLExpression* const expr_;

};  // class GLSLJumpStatement

// ============================================================================

class GLSLSelectionStatement : public GLSLSimpleStatement
{
public:
    GLSLSelectionStatement(GLSLExpression* const cond, GLSLStatement* const trueStmt)
        : GLSLSimpleStatement(GLSLTerminals::ID_IF),
        cond_(cond),
        trueStmt_(trueStmt),
        falseStmt_(0)
    {
    }

    GLSLSelectionStatement(GLSLExpression* const cond, GLSLStatement* const trueStmt,
        GLSLStatement*  const falseStmt)
        : GLSLSimpleStatement(GLSLTerminals::ID_IF),
        cond_(cond),
        trueStmt_(trueStmt),
        falseStmt_(falseStmt)
    {
    }

    virtual ~GLSLSelectionStatement() {
        delete cond_;
        delete trueStmt_;
        delete falseStmt_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_SELECTION_STATEMENT; }

    GLSLExpression* getCondition() { return cond_; }
    GLSLStatement* getTrueStatement() { return trueStmt_; }
    GLSLStatement* getFalseStatement() { return falseStmt_; }

protected:
    GLSLExpression* const cond_;
    GLSLStatement* const trueStmt_;
    GLSLStatement* const falseStmt_;
};  // class GLSLSelectionStatement

// ============================================================================

class GLSLSwitchStatement : public GLSLSimpleStatement
{
public:
    GLSLSwitchStatement(GLSLExpression* const expr)
        : GLSLSimpleStatement(GLSLTerminals::ID_SWITCH),
        expr_(expr),
        stmts_(0)
    {
    }

    GLSLSwitchStatement(GLSLExpression* const expr, GLSLStatementList* const stmts)
        : GLSLSimpleStatement(GLSLTerminals::ID_SWITCH),
        expr_(expr),
        stmts_(stmts)
    {
    }

    virtual ~GLSLSwitchStatement() {
        delete expr_;
        delete stmts_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_SWITCH_STATEMENT; }

    GLSLExpression* getExpression() { return expr_; }
    GLSLStatementList* getStatements() { return stmts_; }

protected:
    GLSLExpression* const expr_;
    GLSLStatementList* const stmts_;
};  // class GLSLSwitchStatement

// ============================================================================

class GLSLWhileStatement : public GLSLSimpleStatement
{
public:
    GLSLWhileStatement(GLSLCondition* const cond, GLSLStatement* const stmt)
        : GLSLSimpleStatement(GLSLTerminals::ID_WHILE),
        cond_(cond),
        stmt_(stmt)
    {
    }

    virtual ~GLSLWhileStatement() {
        delete cond_;
        delete stmt_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_WHILE_STATEMENT; }

    GLSLCondition* getCondition() { return cond_; }
    GLSLStatement* getBody() { return stmt_; }

protected:
    GLSLCondition* const cond_;
    GLSLStatement* const stmt_;
};  // class GLSLWhileStatement

// ============================================================================

class GLSLCompoundStatement : public GLSLStatement
{
public:
    GLSLCompoundStatement(GLSLStatementList* const statements, const bool newScope)
        : GLSLStatement(GLSLTerminals::ID_LBRACKET, newScope),
        statements_(statements)
    {
    }

    virtual ~GLSLCompoundStatement() {
        delete statements_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_COMPOUND_STATEMENT; }

    GLSLStatementList* getStatementList() { return statements_; }

protected:
    GLSLStatementList* const statements_;
};  // class GLSLCompoundStatement

}   // namespace glslparser

}   // namespace voreen

#endif  // VNR_GLSLSTATEMENT
