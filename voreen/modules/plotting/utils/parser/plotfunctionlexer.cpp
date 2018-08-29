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

#include "plotfunctionlexer.h"

#include "plotfunctionterminals.h"
#include "plotfunctiontoken.h"

#include <list>

namespace voreen {

namespace glslparser {

// There are only two keywords: "true" and "false", so this static member is really small.
//
std::map<std::string, Word> PlotFunctionLexer::keywords_ = PlotFunctionTerminals().getKeywords();

PlotFunctionLexer::PlotFunctionLexer(std::istream* const is, const std::map<std::string, Word>& keywords) 
    : Lexer(is, keywords)
    //, readingString_(false)
{
    PlotFunctionLexer::keywords_ = keywords;
}

PlotFunctionLexer::PlotFunctionLexer(const std::string& fileName, const std::map<std::string, Word>& keywords) 
    : Lexer(fileName, keywords)
    //, readingString_(false)
{
    PlotFunctionLexer::keywords_ = keywords;
}


PlotFunctionLexer::~PlotFunctionLexer() {
}

Token* PlotFunctionLexer::scan() {
    Token* token = nextToken();
    if (token)
        return token;
    else
        return new Token(-1);
}

// private methods
//
Token* PlotFunctionLexer::nextToken() {
    Token* token = 0;
    int tokenLine = getCurrentLineNumber();
    int tokenStart = getCurrentCharNumberInLine();

    try {
        do {
            // Skip whitespace including newlines
            //
            do {
                readChar();
            } while (isWhitespace(peek_));

            switch (peek_) {
                case '.':
                    if (isDigit(nextChar()))
                        token = scanNumber();
                    else
                        token = new Token(PlotFunctionTerminals::ID_DOT);
                    break;  // case '.'

                case '+':
                    token = new OperatorToken(PlotFunctionTerminals::ID_PLUS,peek_);
                    break;  // case '+'

                case '-':
                    token = new OperatorToken(PlotFunctionTerminals::ID_DASH,peek_);
                    break;  // case '-'

                case '/':
                    token = new OperatorToken(PlotFunctionTerminals::ID_SLASH,peek_);
                    break;  // case '/'

                case '*':
                    token = new OperatorToken(PlotFunctionTerminals::ID_STAR,peek_);
                    break;  // case '*'

                case '[':
                    token = new BracketToken(PlotFunctionTerminals::ID_LBRACKET,peek_);
                    break;  // case '['

                case ']':
                    token = new BracketToken(PlotFunctionTerminals::ID_RBRACKET,peek_);
                    break;  // case ']'

                case '(':
                    token = new BracketToken(PlotFunctionTerminals::ID_LPAREN,peek_);
                    break;  // case '('

                case ')':
                    token = new BracketToken(PlotFunctionTerminals::ID_RPAREN,peek_);
                    break;  // case ')'

                case '^':
                    token = new OperatorToken(PlotFunctionTerminals::ID_CARET,peek_);
                    break;  // case '^'

                case ':':
                    token = new IsolatorToken(PlotFunctionTerminals::ID_COLON,peek_);
                    break;  // case ':'

                case ';':
                    token = new IsolatorToken(PlotFunctionTerminals::ID_SEMICOLON,peek_);
                    break;  // case ';'

                case ',':
                    token = new IsolatorToken(PlotFunctionTerminals::ID_COMMA,peek_);
                    break;  // case ','

                case '|':
                    token = new IsolatorToken(PlotFunctionTerminals::ID_VERTICAL_BAR,peek_);
                    break;  // case '|'

            }   // switch (peek_)

            if (token == 0) {
                if (isAlpha(peek_))
                    token = scanIdentifier();
                else if (isDigit(peek_) || (peek_ == '.'))
                    token = scanNumber();
            }
        } while (token == 0);
    } catch (std::runtime_error&) {
        if (token != 0) {
            token->setLineNumber(tokenLine);
            token->setCharNumber(tokenStart);
        }
        else if (tokenStart == 0)
            token = new IdentifierToken(PlotFunctionTerminals::ID_EMPTY,""); // EMPTY Expression

        return token;
    }

    if (token != 0) {
        token->setLineNumber(tokenLine);
        token->setCharNumber(tokenStart);
    }
    return token;
}

IdentifierToken* PlotFunctionLexer::scanIdentifier() {
    std::string id;
    int tokenID = PlotFunctionTerminals::ID_VARIABLE;
    if (! isAlpha(peek_))
        return 0;

   while (true) {
        id += peek_;
        if ((! isAlpha(nextChar())) && (! isDigit(nextChar())))
            break;
        readChar();
    }

    if (id.empty())
        return 0;
    else if (id.size() > 1) {
        bool allowed = false;
        // verified the allowed function like log, sin,...
        std::vector<std::string> list = PlotFunctionTerminals::getAllowedFunctionList();
        int index = -1;
        for (size_t i = 0; i < list.size(); ++i) {
            if (id == list[i]) {
                index = static_cast<int>(i);
                allowed = true;
                break;
            }
        }
        if (allowed)
            tokenID = PlotFunctionTerminals::ID_FUNCTION;
        else
            return 0;
        return new FunctionToken(tokenID,id,index);
    }

    return new VariablesToken(tokenID, id);
}

ConstantToken* PlotFunctionLexer::scanNumber() {
    ConstantToken* t = Lexer::scanNumber(PlotFunctionTerminals::ID_UNKNOWN);

    if (t != 0) {
        switch (t->getType()) {
            case ConstantToken::TYPE_FLOAT:
                t->replaceTokenID(PlotFunctionTerminals::ID_FLOATCONST);
                break;

            case ConstantToken::TYPE_UINT:
            case ConstantToken::TYPE_INT:
                t->replaceTokenID(PlotFunctionTerminals::ID_INTCONST);
                break;
        }
    }

    return t;
}

}   // namespace glslparser

}   // namespace voreen
