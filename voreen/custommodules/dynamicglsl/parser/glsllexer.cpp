/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "glsllexer.h"

#include "voreen/core/utils/GLSLparser/annotations/annotationlexer.h"
#include "voreen/core/utils/GLSLparser/glsl/glslterminals.h"

namespace voreen {

namespace glslparser {

GLSLLexer::GLSLLexer(std::istream* const is, const std::map<std::string, Word>& keywords) 
    : Lexer(is, keywords),
    state_(STATE_NORMAL),
    previousTokenID_(-1)
{
    initPredefinedStructs();
}

GLSLLexer::GLSLLexer(const std::string& fileName, const std::map<std::string, Word>& keywords) 
    : Lexer(fileName, keywords),
    state_(STATE_NORMAL),
    previousTokenID_(-1)
{
    initPredefinedStructs();
}

GLSLLexer::~GLSLLexer() {
}

Token* GLSLLexer::scan() {
    Token* next = 0;
    try {
        do {
            next = nextToken();
        } while(next == 0);
    } catch (std::runtime_error& ) {
    }
    return next;
}

// private methods
//

void GLSLLexer::initPredefinedStructs() {
    std::set<std::string>& structs = GLSLTerminals::getPredefinedStructs();
    for (std::set<std::string>::const_iterator it = structs.begin(); it != structs.end(); ++it)
        typeNames_.insert( std::make_pair(*it, Word(GLSLTerminals::ID_TYPE_NAME, *it)) );
}

Token* GLSLLexer::nextToken() {
    switch (state_) {
        case STATE_INSIDE_ANNOTATION_LINE:
            {
                std::ostringstream oss;
                do {
                    readChar();
                    if (isNewline(peek_))
                        break;
                    else
                        oss << peek_;
                } while (true);
                state_ = STATE_NORMAL;
                return scanAnnotation(oss.str(), (previousTokenID_ != GLSLTerminals::ID_SEMICOLON));
            }
            break;

        case STATE_INSIDE_ANNOTATION_BLOCK:
            {
                std::ostringstream oss;
                do {
                    readChar();
                    if ((peek_ == '*') && (nextChar() == '/'))
                        break;
                    else
                        oss << peek_;
                } while (true);
                readChar();
                state_ = STATE_NORMAL;

                // Annotations within block comments can only be leading ones.
                //
                return scanAnnotation(oss.str(), true);
            }
            break;

        case STATE_INSIDE_LINE_COMMENT:
            do {
                readChar();
            } while (! isNewline(peek_));
            state_ = STATE_NORMAL;
            return 0;

        case STATE_INSIDE_BLOCK_COMMENT:
            do {
                readChar();
            } while ((peek_ != '*') || (nextChar() != '/'));
            readChar();
            state_ = STATE_NORMAL;
            // no break here

        // Detect newline characters and omit whitespaces
        //
        case STATE_NORMAL:
            do {
                readChar();

                // In case that the previously returned token was a semicolon, this information
                // is erased in order to not identify any possibly following annotation as
                // a trailing one, for it does not occure in the same line as the semicolon
                // anymore!
                //
                if ((previousTokenID_ == GLSLTerminals::ID_SEMICOLON) && (isNewline(peek_)))
                    previousTokenID_ = -1;
            } while ((isWhitespace(peek_)) || (peek_ == 0));
            break;
    }

    int tokenLine = getCurrentLineNumber();
    int tokenStart = getCurrentCharNumberInLine();
    Token* token = 0;
    switch (peek_) {
        case '.':
            if (isDigit(nextChar()))
                token = scanNumber();
            else
                token = new Token(GLSLTerminals::ID_DOT);
            break;  // case '.'

        case '+':
            if (readCharConditioned('+'))
               token = new Token(GLSLTerminals::ID_INC_OP);
            else if (readCharConditioned('='))
                token = new Token(GLSLTerminals::ID_ADD_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_PLUS);
            break;  // case '+'

        case '-':
            if (readCharConditioned('-') == true)
               token = new Token(GLSLTerminals::ID_DEC_OP);
            else if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_SUB_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_DASH);
            break;  // case '-'

        case '/':
            if (readCharConditioned('/'))
                if (readCharConditioned('$'))
                    state_ = STATE_INSIDE_ANNOTATION_LINE;
                else
                    state_ = STATE_INSIDE_LINE_COMMENT;
            else if (readCharConditioned('*'))
                if (readCharConditioned('$'))
                    state_ = STATE_INSIDE_ANNOTATION_BLOCK;
                else
                    state_ = STATE_INSIDE_BLOCK_COMMENT;
            else if (readCharConditioned('='))
                token = new Token(GLSLTerminals::ID_DIV_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_SLASH);
            break;  // case '/'

        case '*':
            if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_MUL_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_STAR);
            break;  // case '*'

        case '%':
            if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_MOD_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_PERCENT);
            break;  // case '%'

        case '<':
            if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_LE_OP);
            else if (readCharConditioned('<') == true) {
                if (readCharConditioned('=') == true)
                    token = new Token(GLSLTerminals::ID_LSHIFT_ASSIGN);
                else
                    token = new Token(GLSLTerminals::ID_LSHIFT_OP);
            } else
                token = new Token(GLSLTerminals::ID_LANGLE);
            break;  // case '<'

        case '>':
            if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_GE_OP);
            else if (readCharConditioned('>') == true) {
                if (readCharConditioned('=') == true)
                    token = new Token(GLSLTerminals::ID_RSHIFT_ASSIGN);
                else
                    token = new Token(GLSLTerminals::ID_RSHIFT_OP);
            } else
                token = new Token(GLSLTerminals::ID_RANGLE);
            break;  // case '>'

        case '[':
            token = new Token(GLSLTerminals::ID_LBRACKET);
            break;

        case ']':
            token = new Token(GLSLTerminals::ID_RBRACKET);
            break;

        case '(':
            token = new Token(GLSLTerminals::ID_LPAREN);
            break;

        case ')':
            token = new Token(GLSLTerminals::ID_RPAREN);
            break;

        case '{':
            token = new Token(GLSLTerminals::ID_LBRACE);
            break;

        case '}':
            token = new Token(GLSLTerminals::ID_RBRACE);
            break;

        case '^':
            if (readCharConditioned('^') == true)
                token = new Token(GLSLTerminals::ID_XOR_OP);
            else if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_XOR_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_CARET);
            break;  // case '^'

        case '|':
            if (readCharConditioned('|') == true)
                token = new Token(GLSLTerminals::ID_OR_OP);
            else if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_OR_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_VERTICAL_BAR);
            break;  // case '|'

        case '&':
            if (readCharConditioned('&') == true)
                token = new Token(GLSLTerminals::ID_AND_OP);
            else if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_AND_ASSIGN);
            else
                token = new Token(GLSLTerminals::ID_AMPERSAND);
            break;  // case '&'

        case '~':
            token = new Token(GLSLTerminals::ID_TILDE);
            break;

        case '=':
            if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_EQ_OP);
            else
                token = new Token(GLSLTerminals::ID_EQUAL);
            break;  // case '='

        case '!':
            if (readCharConditioned('=') == true)
                token = new Token(GLSLTerminals::ID_NE_OP);
            else
                token = new Token(GLSLTerminals::ID_BANG);
            break;  // case '!'

        case ':':
            token = new Token(GLSLTerminals::ID_COLON);
            break;

        case ';':
            token = new Token(GLSLTerminals::ID_SEMICOLON);
            break;

        case ',':
            token = new Token(GLSLTerminals::ID_COMMA);
            break;

        case '?':
            token = new Token(GLSLTerminals::ID_QUESTION);
            break;

        default:
            break;
    }   // switch (peek_)

    if ((token == 0) && (state_ == STATE_NORMAL)) {
        if (isDigit(peek_))
            token = scanNumber();
        else if (isAlpha(peek_)) {
            IdentifierToken* const identifier = scanIdentifier(GLSLTerminals::ID_IDENTIFIER);

            // When having read an identifier, first determine whether it is a reserved
            // keyword or a previously declared type name which have special token IDs.
            // In this case, do not return an identifier token but replace it with a
            // copy of the more specialized version instead.
            //
            if (identifier != 0) {
                KeywordMap::const_iterator itKeyword = keywords_.find(identifier->getValue());

                if (itKeyword != keywords_.end()) {
                    delete identifier;
                    token = itKeyword->second.getCopy();
                } else {
                    TypeNameMap::const_iterator itTypeName = typeNames_.find(identifier->getValue());

                    if (itTypeName != typeNames_.end()) {
                        delete identifier;
                        token = itTypeName->second.getCopy();
                    } else
                        token = identifier;
                }
            }   // if (identifier != 0)
        }   // else if (isAlpha(peek_

        if ((token == 0) && (state_ == STATE_NORMAL))
            token = new Word(GLSLTerminals::ID_UNKNOWN, std::string(&peek_, 1));
    }

    // Remember the ID of the returning token in order to adjust following
    // IDENTIFIER token after having previously returned a DOT or STRUCT token.
    //
    if (token != 0) {
        if (token->getTokenID() == GLSLTerminals::ID_IDENTIFIER) {
            switch (previousTokenID_) {
                case GLSLTerminals::ID_DOT:
                    token->replaceTokenID(GLSLTerminals::ID_FIELDSELECTION);
                    break;

                // A struct has been declared, so insert the new identifer into
                // the map of type names
                //
                case GLSLTerminals::ID_STRUCT:
                    if (IdentifierToken* const id = dynamic_cast<IdentifierToken*>(token)) {
                        typeNames_.insert( std::make_pair(id->getValue(),
                            Word(GLSLTerminals::ID_TYPE_NAME, id->getValue())) );
                    }
                    break;
            }
        }

        token->setLineNumber(tokenLine);
        token->setCharNumber(tokenStart);
        previousTokenID_ = token->getTokenID();
    }

    return token;
}

// private methods
//

AnnotationToken* GLSLLexer::scanAnnotation(const std::string& annotation, const bool isLeading) {
    if (annotation.empty())
        return 0;

    const int tokenID = ((isLeading) ? GLSLTerminals::ID_LEADING_ANNOTATION :
        GLSLTerminals::ID_TRAILING_ANNOTATION);

    std::istringstream is(annotation);
    AnnotationLexer anLexer(&is);
    AnnotationToken* const anToken = dynamic_cast<AnnotationToken* const>(anLexer.scan());

    if (anToken != 0)
        anToken->replaceTokenID(tokenID);

    return anToken;
}

ConstantToken* GLSLLexer::scanNumber() {
    ConstantToken* t = Lexer::scanNumber(GLSLTerminals::ID_UNKNOWN);

    if (t != 0) {
        switch (t->getType()) {
            case ConstantToken::TYPE_FLOAT:
                t->replaceTokenID(GLSLTerminals::ID_FLOATCONST);
                break;

            case ConstantToken::TYPE_INT:
                t->replaceTokenID(GLSLTerminals::ID_INTCONST);
                break;

            case ConstantToken::TYPE_UINT:
                t->replaceTokenID(GLSLTerminals::ID_UINTCONST);
                break;
        }
    }

    return t;
}

}   // namespace glslparser

}   // namespace voreen
