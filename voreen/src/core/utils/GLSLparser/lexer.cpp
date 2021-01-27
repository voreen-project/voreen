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

#include "voreen/core/utils/GLSLparser/lexer.h"

#include <fstream>
#include <new>

namespace voreen {

namespace glslparser {

Lexer::Lexer(std::istream* const is, const std::map<std::string, Word>& keywords) 
    : stream_(is),
    fileName_("[stream]"),
    peek_(0),
    keywords_(keywords),
    numLine_(1),
    numCharInLine_(0),
    freeStream_(false)
{
    if (stream_->fail() || stream_->bad())
        throw std::bad_alloc();

    stream_->seekg(0, std::istream::beg);
}

Lexer::Lexer(const std::string& fileName, const std::map<std::string, Word>& keywords) 
    : stream_(new std::ifstream(fileName.c_str(), std::ios_base::in | std::ios_base::binary)),
    fileName_(fileName),
    peek_(0),
    keywords_(keywords),
    numLine_(1),
    numCharInLine_(0),
    freeStream_(true)
{
    if (stream_->fail() || stream_->bad())
        throw std::bad_alloc();

    stream_->seekg(0, std::istream::beg);
}

Lexer::~Lexer() {
    if (freeStream_ == true) {
        std::ifstream* ifs = dynamic_cast<std::ifstream*>(stream_);
        if (ifs != 0)
            ifs->close();
        delete stream_;
        stream_ = 0;
    }
}

bool Lexer::isAlpha(const char ref) {
    if ((ref >= 'A') && (ref <= 'Z'))
        return true;
    else if ((ref >= 'a') && (ref <= 'z'))
        return true;
    else if (ref == '_')
        return true;
    return false;
}

bool Lexer::isDigit(const char ref) {
    switch (ref) {
        case '0': case '1':
        case '2': case '3':
        case '4': case '5':
        case '6': case '7':
        case '8': case '9':
            return true;
        default:
            break;
    }
    return false;
}

bool Lexer::isOctalDigit(const char ref) {
    switch (ref) {
        case '0': case '1': case '2': case '3':
        case '4': case '5': case '6': case '7':
            return true;
        default:
            break;
    }
    return false;
}

bool Lexer::isHexDigit(const char ref) {
    if (isDigit(ref) == true)
        return true;
    switch (ref) {
        case 'a': case 'A':
        case 'b': case 'B':
        case 'c': case 'C':
        case 'd': case 'D':
        case 'e': case 'E':
        case 'f': case 'F':
            return true;
        default:
            break;
    }
    return false;
}

bool Lexer::isNewline(const char ref) {
    switch (ref) {
        case '\n':
        case '\r':
            return true;
    }

    return false;
}

bool Lexer::isWhitespace(const char ref) {
    switch (ref) {
        case ' ':   // space
        case '\t':  // horizontal tab
        case '\v':  // vertical tab
        case '\f':  // form feed
            return true;
        default:
            break;
    }
    return Lexer::isNewline(ref);
}

void Lexer::reset() {
    stream_->seekg(0, std::ios_base::beg);
    numLine_ = 1;
    numCharInLine_ = 0;
    peek_ = 0;
}

// protected methods
//

char Lexer::readChar() {
    stream_->get(peek_);
    ++numCharInLine_;

    if (stream_->fail()) {
        peek_ = 0;
        throw std::runtime_error("Lexer::readChar(): No more characters to read from the stream!");
    }

    if (peek_ == '\n') {
        ++numLine_;
        numCharInLine_ = 1;
    }

    return peek_;
}

char Lexer::nextChar() {
    char next = stream_->peek();
    if (stream_->fail())
        next = 0;

    return next;
}

bool Lexer::readCharConditioned(const char ref) {
    char old = peek_;
    int numLine = numLine_;
    int numChar = numCharInLine_;

    if (readChar() != ref) {
        // restore previous state in case the condition is not true
        //
        stream_->seekg(-1, std::ios_base::cur);

        numLine_ = numLine;
        numCharInLine_ = numChar;
        peek_ = old;
        return false;
    }

    return true;
}

bool Lexer::readCharConditioned(bool (*pred)(const char)) {
    if (pred == 0)
        return false;

    char old = peek_;
    int numLine = numLine_;
    int numChar = numCharInLine_;

    if ((*pred)(readChar()) == false) {
        // restore previous state in case the condition is not true
        //
        stream_->seekg(-1, std::ios_base::cur);

        peek_ = old;
        numLine_ = numLine;
        numCharInLine_ = numChar;
        return false;
    }

    return true;
}

IdentifierToken* Lexer::scanIdentifier(const int tokenID) {
    std::string id;
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

    return new IdentifierToken(tokenID, id);
}

ConstantToken* Lexer::scanNumber(const int tokenID) {
    std::string numberString("");
    ConstantToken::Base base = ConstantToken::BASE_10;
    ConstantToken::ConstantType type = ConstantToken::TYPE_INT;
    bool hasExponent = false;
    bool (*pred)(const char ref) = &isDigit;

    size_t pos = 0;
    for (bool atEnd = false; (! atEnd); ++pos) {
        if ((pred != 0) && ((*pred)(peek_) == true)) {
            numberString += peek_;

            if ((pos == 0) && (peek_ == '0')) {
                if ((nextChar() == 'x') || (nextChar() == 'X')) {
                    numberString += readChar();
                    ++pos;
                    pred = &isHexDigit;
                    base = ConstantToken::BASE_16;
                } else if (nextChar() == '.')
                    type = ConstantToken::TYPE_FLOAT;
                else {
                    pred = &isOctalDigit;
                    base = ConstantToken::BASE_8;
                }
            }

            switch (nextChar()) {
                case 'u':
                case 'U':
                    atEnd = true;
                    if (type == ConstantToken::TYPE_INT) {
                        numberString += readChar();
                        ++pos;
                        type = ConstantToken::TYPE_UINT;
                    }
                    break;

                case 'f':
                case 'F':
                    atEnd = true;
                    if (type == ConstantToken::TYPE_FLOAT) {
                        numberString += readChar();
                        ++pos;
                    }
                    break;

                case '.':
                    if (base == ConstantToken::BASE_10) {
                        type = ConstantToken::TYPE_FLOAT;
                        numberString += readChar();
                        ++pos;
                    }
                    break;

                case 'e':
                case 'E':
                    if (type == ConstantToken::TYPE_FLOAT) {
                        hasExponent = true;
                        numberString += readChar();
                        ++pos;
                    }
                    break;

                default:
                    if (! (*pred)(nextChar()))
                        atEnd = true;
                    break;
            }   // switch
        }
        else if ((peek_ == 'f') && (type == ConstantToken::TYPE_FLOAT))
            numberString += peek_;
        else if ((peek_ == '.') && (base == ConstantToken::BASE_10))
        {
            // Accept the dot '.' only, if it is not the first character or
            // if a digit is following to avoid invalid ".f" strings.
            //
            if ((pos > 0) || (isDigit(nextChar()) == true)) {
                type = ConstantToken::TYPE_FLOAT;
                numberString += peek_;
            }
        }
        else if (((peek_ == '+') || (peek_ == '-')) && (hasExponent == true)) {
            numberString += peek_;
        } else
            break;

        if (! atEnd)
            readChar();
    }   // for (pos

    if (numberString.empty() == true)
        return 0;

    return new ConstantToken(tokenID, numberString, type, base);
}

}   // namespace glslparser

}   // namespace voreen
