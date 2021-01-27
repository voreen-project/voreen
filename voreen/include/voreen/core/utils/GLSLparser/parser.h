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

#ifndef VRN_PARSER_H
#define VRN_PARSER_H

#include "voreen/core/utils/GLSLparser/lexer.h"
#include "voreen/core/utils/GLSLparser/parseraction.h"
#include "voreen/core/utils/GLSLparser/parsetreenode.h"

#include <sstream>
#include <stack>

namespace voreen {

namespace glslparser {

class Parser {
public:
    Parser(Lexer* const lexer, const bool deleteLexer = false);
    virtual ~Parser();

    bool hasLexer() const { return (lexer_ != 0); }

    const std::ostringstream& getLog() const { return log_; }

    /**
     * Parses the given input.
     *
     * @return  Pointer to the root node of the parse tree, which is identified with
     *          that root node. The node may have children.
     */
    ParseTreeNode* parse();

    void setDebugging(const bool debug) { debugging_ = debug; }

protected:
    class ParserSymbol {
    public:
        ParserSymbol(const int symbolID)
            : symbolID_(symbolID),
            token_(0)
        {
        }

        explicit ParserSymbol(Token* const token)
            : symbolID_(token->getTokenID()),
            token_(token)
        {
        }

        ~ParserSymbol() {
            delete token_;
        }

        int getSymbolID() const { return symbolID_; }
        Token* getToken() const { return token_; }
        bool isTerminal() const { return token_ != 0; }


    private:
        ParserSymbol(const ParserSymbol&);
        ParserSymbol& operator=(const ParserSymbol&);

    private:
        const int symbolID_;
        Token* const token_;
    };

    friend bool operator==(const Parser::ParserSymbol& lhs, const Parser::ParserSymbol& rhs);
    friend bool operator!=(const Parser::ParserSymbol& lhs, const Parser::ParserSymbol& rhs);

protected:
    /**
     * Action-Table for this parser.
     *
     * The caller has to delete the returned object using C++ delete-operator.
     */
    virtual ParserAction* action(const int state, const int symbol) const = 0;

    /**
     * Nodes-stack operations to perform when reducing the production related
     * to the given ID. This actually builds the ParseTreeNodes which hold the
     * semantics of the parsed input.
     */
    virtual void expandParseTree(const int productionID,
        const std::vector<Parser::ParserSymbol*>& reductionBody) = 0;

    /**
     * Returns a ProductionStub which is consists of the ID of the head
     * symbol and the IDs of the symbol from the body. The production
     * stub is used to check integrity while reducing the states-stack.
     *
     * This method is called by <code>Parser::reduce()</code> and calls
     * <code>expandParseTree()</code> itself.
     */
    virtual ProductionStub* findProduction(const int productionID) const = 0;

    /**
     * Goto-Table for this parser.
     * Determines the transition for the parser's state according to the
     * given state and the symbol (the symbol ID to be more precise).
     */
    virtual int gotoState(const int state, const int symbol) const = 0;

    /**
     * Returns the corresponding lexeme for the symbol with the given ID
     * if such a symbol exists. Otherwise the returned string is empty.
     */
    virtual std::string symbolID2String(const int symbolID) const = 0;

    /**
     * Pops the top node from the nodes-stack, casts it to the given template
     * parameter and returns it, if the statck is not empty. If the nodes-stack
     * is empty or the topmost node cannot be casted to the desired type, the
     * method returns NULL.
     */
    template<class NODE>
    NODE* popNode() {
        if (nodes_.empty())
            return 0;

        NODE* n = dynamic_cast<NODE*>(nodes_.top());

        // Prevents memory leaks in case of dynamic_cast failures.
        // The topmost pointer would not be accessible anymore, because it would
        // be poped, but not returned...
        //
        if (n == 0)
            delete nodes_.top();

        nodes_.pop();
        return n;
    }

    ParseTreeNode* popNode();

    /**
     * Pushes the given node onto the nodes-stack, if it is not NULL.
     */
    void pushNode(ParseTreeNode* const node);

    /**
     * While reducing a production to its head, this method is called in case that
     * the considered symbol's ID does not match the required one from the body of
     * the prodcution. If the symbol within the body of the production is a wildcard
     * symbol, the reduction will continue instead of raising a fatal error.
     *
     * This is useful in case of the GLSL preprocessor where there productions
     * containing a terminal symbol called TOKEN which is actually a wildcard for all
     * possible tokens which might have been emitted by the driving lexical analysator
     * (lexer).
     */
    virtual bool isProxySymbol(const int /*symbolID*/, const int /*originalID*/) const {
        return false;
    }

    /**
     * Called within <code>parse()</code> to provide the next token for parsing.
     * The default implementation takes the next token from the lexer. If the
     * returned token is NULL, the parser will stop parsing after it has performed
     * all possible reductions left.
     */
    virtual Token* nextToken();

private:
    Parser(const Parser&);
    Parser& operator=(const Parser&);

    /**
     * Performs deletion of stack-objects (e.g. nodes-stack) and resets the lexer.
     * Also called by the dtor.
     */
    void cleanup();

    void logError(ParserActionError* const error, Token* const token);

    /**
     * Dumps the current parser state (states-stack, symbols-stack and the stack
     * containing all ParserTreeNode nodes which haven been created so far) for
     * debug purposes to the log.
     */
    void logState();

    /**
     * Performs a reduction of the parser stack depending on the given ParserActionReduce
     * object. If the reduction is successful, the method returns true, otherwise false
     * will be returned.
     *
     * @param   r   ParserAction hold a ProductionStub object to perform the stack reduction.
     * @return  true if the reduction was performed successfully, false otherwise.
     */
    bool reduce(ParserActionReduce* const r);

private:
    bool debugging_;

    Lexer* const lexer_;
    const bool deleteLexer_;

    std::stack<int> states_;
    std::stack<ParserSymbol*> symbols_;
    std::stack<ParseTreeNode*> nodes_;

    std::ostringstream log_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
