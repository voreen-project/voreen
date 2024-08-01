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

#include "voreen/core/utils/GLSLparser/annotations/annotationlexer.h"

#include "voreen/core/utils/GLSLparser/annotations/annotationterminals.h"
#include "voreen/core/utils/GLSLparser/glsl/glsltoken.h"

#include <list>

namespace voreen {

namespace glslparser {

// There are only two keywords: "true" and "false", so this static member is really small.
//
std::map<std::string, Word> AnnotationLexer::keywords_ = AnnotationTerminals().getKeywords();

AnnotationLexer::AnnotationLexer(std::istream* const is) 
    : Lexer(is, AnnotationLexer::keywords_),
    readingString_(false)
{
}

AnnotationLexer::~AnnotationLexer() {
}

Token* AnnotationLexer::scan() {
    std::list<Token*> parts;

    for (Token* token = nextToken(); token != 0; token = nextToken())
        parts.push_back(token);

    return new AnnotationToken(-1, parts);
}

// private methods
//

Token* AnnotationLexer::nextToken() {
    Token* token = 0;
    int tokenLine = getCurrentLineNumber();
    int tokenStart = getCurrentCharNumberInLine();

    try {
        do {
            if (readingString_) {
                readChar();
                if (peek_ == '"') {
                    readingString_ = false;
                    return new Token(AnnotationTerminals::ID_QUOTE);
                } else {
                    std::string str;
                    do {
                        str += peek_;
                        if (nextChar() != '"')
                            readChar();
                        else
                            break;
                    } while (true);
                    return new StringToken(AnnotationTerminals::ID_STRING, str);
                }
            }

            // Skip whitespace including newlines
            //
            do {
                readChar();
            } while (isWhitespace(peek_));

            switch (peek_) {
                case '@':
                    token = new Token(AnnotationTerminals::ID_AT);
                    break;

                case '=':
                    token = new Token(AnnotationTerminals::ID_EQUAL);
                    break;

                case '[':
                    token = new Token(AnnotationTerminals::ID_LBRACKET);
                    break;

                case ']':
                    token = new Token(AnnotationTerminals::ID_RBRACKET);
                    break;

                case ',':
                    token = new Token(AnnotationTerminals::ID_COMMA);
                    break;

                case '\"':
                    token = new Token(AnnotationTerminals::ID_QUOTE);
                    readingString_ = true;
                    break;

                case '\'':
                    token = new Token(AnnotationTerminals::ID_APOSTROPHE);
                    break;
            }   // switch (peek_)

            if (token == 0) {
                if (isAlpha(peek_))
                    token = scanIdentifier(AnnotationTerminals::ID_IDENTIFIER);
                else if (isDigit(peek_) || (peek_ == '.'))
                    token = scanNumber();
            }
        } while (token == 0);
    } catch (std::runtime_error&) {
        if (token != 0) {
            token->setLineNumber(tokenLine);
            token->setCharNumber(tokenStart);
        }

        return token;
    }

    if (token != 0) {
        token->setLineNumber(tokenLine);
        token->setCharNumber(tokenStart);
    }
    return token;
}

ConstantToken* AnnotationLexer::scanNumber() {
    ConstantToken* t = Lexer::scanNumber(AnnotationTerminals::ID_UNKNOWN);

    if (t != 0) {
        switch (t->getType()) {
            case ConstantToken::TYPE_FLOAT:
                t->replaceTokenID(AnnotationTerminals::ID_FLOATCONST);
                break;

            case ConstantToken::TYPE_UINT:
            case ConstantToken::TYPE_INT:
                t->replaceTokenID(AnnotationTerminals::ID_INTCONST);
                break;
        }
    }

    return t;
}

}   // namespace glslparser

}   // namespace voreen
