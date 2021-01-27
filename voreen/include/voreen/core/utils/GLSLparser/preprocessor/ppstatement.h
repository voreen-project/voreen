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

#ifndef VRN_PPSTATEMENT_H
#define VRN_PPSTATEMENT_H

#include "voreen/core/utils/GLSLparser/preprocessor/ppexpression.h"

namespace voreen {

namespace glslparser {

class Statement : public ParseTreeNode {
public:
    Statement(const int symbolID)
        : ParseTreeNode(symbolID)
    {
    }

    virtual ~Statement() {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_STATEMENT; }
};

// ============================================================================

class IdentifierList : public Statement {
public:
    IdentifierList()
        : Statement(PreprocessorTerminals::ID_IDENTIFIER)
    {
    }

    IdentifierList(const std::string& identifier)
        : Statement(PreprocessorTerminals::ID_IDENTIFIER)
    {
        addIdentifier(identifier);
    }

    virtual ~IdentifierList() {
    }

    void addIdentifier(const std::string& identifier) {
        if (identifier != "")
            identifiers_.push_back(identifier);
    }

    const std::vector<std::string>& getIdentifiers() const { return identifiers_; }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_IDENTIFIER_LIST; }

protected:
    std::vector<std::string> identifiers_;
};

// ============================================================================

class DefineDirective : public Statement {
public:
    DefineDirective(const std::string& identifier, TokenList* const list,
                    const bool isFunction, IdentifierList* const formals = 0)
        : Statement(PreprocessorTerminals::ID_DEFINE),
        identifier_(identifier),
        formals_(formals),
        body_(list),
        isFunction_(isFunction)
    {
    }

    virtual ~DefineDirective() {
        delete formals_;
        //delete body_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_DEFINE; }

    const std::string& getIdentifier() const { return identifier_; }
    IdentifierList* getFormals() const { return formals_; }
    TokenList* getBody() const { return body_; }
    bool isFunction() const { return isFunction_; }

protected:
    std::string identifier_;
    IdentifierList* const formals_;
    TokenList* const body_;
    bool isFunction_;
};

// ============================================================================

class ErrorDirective : public Statement {
public:
    ErrorDirective(TokenList* const tokens)
        : Statement(PreprocessorTerminals::ID_ERROR),
        tokenList_(tokens)
    {
    }

    virtual ~ErrorDirective() {
        delete tokenList_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_ERROR; }

protected:
    TokenList* const tokenList_;
};

// ============================================================================

class ExtensionDirective : public Statement {
public:
    ExtensionDirective(TokenList* const tokens)
        : Statement(PreprocessorTerminals::ID_EXTENSION),
        tokenList_(tokens)
    {
    }

    virtual ~ExtensionDirective() {
        delete tokenList_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_EXTENSION; }

protected:
    TokenList* const tokenList_;
};

// ============================================================================

class IncludeDirective : public Statement {
public:
    IncludeDirective(const std::string& filename)
        : Statement(PreprocessorTerminals::ID_INCLUDE),
        filename_(filename)
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_INCLUDE; }

    const std::string& getFileName() const { return filename_; }

protected:
    std::string filename_;
};

// ============================================================================

class LineDirective : public Statement {
public:
    LineDirective(const ConstantToken& lineNo)
        : Statement(PreprocessorTerminals::ID_LINE),
        lineNumber_(lineNo),
        fileName_(PreprocessorTerminals::ID_STRING, "")
    {
    }

    LineDirective(const ConstantToken& lineNo, const StringToken& fileName)
        : Statement(PreprocessorTerminals::ID_LINE),
        lineNumber_(lineNo),
        fileName_(fileName)
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_LINE; }

    int getLineNumber() const { return lineNumber_.convert<int>(); }

protected:
    ConstantToken lineNumber_;
    StringToken fileName_;
};

// ============================================================================

class PragmaDirective : public Statement {
public:
    PragmaDirective(TokenList* const tokens)
        : Statement(PreprocessorTerminals::ID_PRAGMA),
        tokenList_(tokens)
    {
    }

    virtual ~PragmaDirective() {
        delete tokenList_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_PRAGMA; }

protected:
    TokenList* const tokenList_;
};

// ============================================================================

class UndefineDirective : public Statement {
public:
    UndefineDirective(const std::string& identifier)
        : Statement(PreprocessorTerminals::ID_UNDEF),
        identifier_(identifier)
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_UNDEFINE; }

    const std::string& getIdentifier() const { return identifier_; }

protected:
    std::string identifier_;
};

// ============================================================================

class VersionDirective : public Statement {
public:
    VersionDirective(const ConstantToken& version)
        : Statement(PreprocessorTerminals::ID_VERSION),
        version_(version),
        profile_(PreprocessorTerminals::ID_IDENTIFIER, "")
    {
    }

    VersionDirective(const ConstantToken& version, const IdentifierToken& profile)
        : Statement(PreprocessorTerminals::ID_VERSION),
        version_(version),
        profile_(profile)
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_VERSION; }

    int getVersion() const { return version_.convert<int>(); }

protected:
    ConstantToken version_;
    IdentifierToken profile_;
};

// ============================================================================

class ConditionalDirective : public Statement {
public:
    ConditionalDirective(const int symbolID, Expression* const condition,
                         ParseTreeNode* const fulfillment)
        : Statement(symbolID),
        condition_(condition),
        true_(fulfillment),
        false_(0)
    {
    }

    ConditionalDirective(const int symbolID, Expression* const condition,
                         ParseTreeNode* const fulfillment, ParseTreeNode* const alternative)
        : Statement(symbolID),
        condition_(condition),
        true_(fulfillment),
        false_(alternative)
    {
    }

    virtual ~ConditionalDirective() = 0;

    Expression* getCondition() const { return condition_; }
    ParseTreeNode* getTrue() const { return true_; }
    ParseTreeNode* getFalse() const { return false_; }

    void setTrue(ParseTreeNode* const exp, const bool deleteExisting = true) {
        if (deleteExisting)
            delete true_;
        true_ = exp;
    }

    void setFalse(ParseTreeNode* const exp, const bool deleteExisting = true) {
        if (deleteExisting)
            delete false_;
        false_ = exp;
    }

protected:
    Expression* const condition_;
    ParseTreeNode* true_;
    ParseTreeNode* false_;
};

// ============================================================================

class IfDirective : public ConditionalDirective {
public:
    IfDirective(Expression* const condition, ParseTreeNode* const fulfillment)
        : ConditionalDirective(PreprocessorTerminals::ID_IF, condition, fulfillment)
    {
    }

    IfDirective(Expression* const condition, ParseTreeNode* const fulfillment,
                ParseTreeNode* const alternative)
        : ConditionalDirective(PreprocessorTerminals::ID_IF, condition, fulfillment, alternative)
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_IF; }
};

// ============================================================================

class IfdefDirective : public ConditionalDirective {
public:
    IfdefDirective(IdentifierToken* const token, ParseTreeNode* const fulfillment)
        : ConditionalDirective(PreprocessorTerminals::ID_IFDEF,
        new DefinedOperator(token), fulfillment)
    {
    }

    IfdefDirective(IdentifierToken* const token, ParseTreeNode* const fulfillment,
        ParseTreeNode* const alternative)
        : ConditionalDirective(PreprocessorTerminals::ID_IF,
        new DefinedOperator(token), fulfillment, alternative)
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_IFDEF; }
};

// ============================================================================

class IfndefDirective : public ConditionalDirective {
public:
    IfndefDirective(IdentifierToken* const token, Statement* const fulfillment)
        : ConditionalDirective(PreprocessorTerminals::ID_IFNDEF,
        new UnaryExpression(new Token(PreprocessorTerminals::ID_OP_NOT),
        new DefinedOperator(token)), 0, fulfillment)
    {
    }

    IfndefDirective(IdentifierToken* const token, Statement* const fulfillment,
        Statement* const alternative)
        : ConditionalDirective(PreprocessorTerminals::ID_IFNDEF,
        new UnaryExpression(new Token(PreprocessorTerminals::ID_OP_NOT),
        new DefinedOperator(token)), alternative, fulfillment)
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_IFNDEF; }
};

}   // namespace glslparser

}   // namespace voreen

#endif
