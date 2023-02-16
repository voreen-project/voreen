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

#ifndef VRN_PPPARSETREENODE_H
#define VRN_PPPARSETREENODE_H

#include "voreen/core/utils/GLSLparser/parsetreenode.h"

#include "voreen/core/utils/GLSLparser/glsl/glsltoken.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppterminals.h"

#include <list>
#include <sstream>

namespace voreen {

namespace glslparser {

class PreprocessorSymbol;

/**
 * A struct to facilitate the runtime-identification of different
 * ParseTreeNode sub-classes for the visitor design pattern.
 */
struct PreprocessorNodeTypes {
    enum {
        NODE_EXPRESSION,
            NODE_ARITHMETIC,
            NODE_BINARY_EXPRESSION,
            NODE_BINARY_LOGICAL,
            NODE_DEFINED_OPERATOR,
            NODE_INT_CONSTANT,
            NODE_LOGICAL_EXPRESSION,
            NODE_PARENTHESIS,
            NODE_UNARY_EXPRESSION,
        NODE_EXPRESSIONLIST,
        NODE_MACRO,
        NODE_STATEMENT,
            NODE_DEFINE,
            NODE_ERROR,
            NODE_EXTENSION,
            NODE_IDENTIFIER_LIST,
            NODE_INCLUDE,
            NODE_LINE,
            NODE_PRAGMA,
            NODE_UNDEFINE,
            NODE_VERSION,
            NODE_CONDITIONAL,
                NODE_IF,
                NODE_IFDEF,
                NODE_IFNDEF,
        NODE_TEXT,
        NODE_TOKEN_LIST
    };
};

// ============================================================================

class TokenList : public ParseTreeNode {
public:
    TokenList()
        : ParseTreeNode(PreprocessorTerminals::ID_TOKEN)
    {
    }

    TokenList(Token* const token)
        : ParseTreeNode(PreprocessorTerminals::ID_TOKEN)
    {
        addToken(token);
    }

    virtual ~TokenList() {
        for (std::list<Token*>::iterator it = tokens_.begin(); it != tokens_.end(); ++it)
            delete *it;
        tokens_.clear();
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_TOKEN_LIST; }
    const std::list<Token*>& getTokens() const { return tokens_; }

    std::list<Token*> getCopy() {
        std::list<Token*> copy;
        for (std::list<Token*>::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it)
            copy.push_back((*it)->getCopy());
        return copy;
    }

    void addToken(Token* const token) {
        if (token != 0)
            tokens_.push_back(token->getCopy());
    }

    std::string toString(const PreprocessorTerminals& terminals) const {
        std::ostringstream oss;
        
        bool isFirst = true; //  FIX Indentation
        for (std::list<Token*>::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it) {
            if (isFirst)
              isFirst = false;
            else
              oss << " ";
            GenericToken<std::string>* const strToken =
                dynamic_cast<GenericToken<std::string>* const>(*it);
            if (strToken != 0)
                oss << strToken->getValue();
            else
                oss << terminals.findLexeme((*it)->getTokenID());
        }
        return oss.str();
    }

protected:
    std::list<Token*> tokens_;
};

// ============================================================================

class TextNode : public ParseTreeNode {
public:
    TextNode(const TextToken& token)
        : ParseTreeNode(token.getTokenID()),
        text_(token.getValue())
    {
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_TEXT; }

    virtual void print(std::ostream& os) const {
        os << text_;
    }

    const std::string& getRawText() const { return text_; }

protected:
    std::string text_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
