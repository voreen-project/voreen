/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_LEXER_H
#define VRN_LEXER_H

#include "voreen/core/utils/GLSLparser/token.h"

#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace voreen {

namespace glslparser {

class Lexer {
public:
    Lexer(std::istream* const is, const std::map<std::string, Word>& keywords);

    Lexer(const std::string& fileName, const std::map<std::string, Word>& keywords);

    virtual ~Lexer();

    int getCurrentLineNumber() const { return static_cast<int>(numLine_); }
    int getCurrentCharNumberInLine() const { return static_cast<int>(numCharInLine_); }

    /**
     * Decides whether the given character is one of those
     * [A-Za-z_]
     */
    static bool isAlpha(const char ref);

    /**
     * Decides whether the given character is one of those
     * [0-9]
     */
    static bool isDigit(const char ref);

    static bool isHexDigit(const char ref);

    static bool isOctalDigit(const char ref);

    /**
     * Decides whether the given character is one of those
     * [\n\r] (line feed or carriage return)
     */
    static bool isNewline(const char ref);

    /**
     * Decides whether the given character is one of those
     * [ \t\v\f\r\n] (space, horizontal tab, vertical tab,
     * form feed, carriage return or line feed)
     */
    static bool isWhitespace(const char ref);

    /**
     * Resets the lexer to its initial state to read the input stream from
     * the beginning.
     */
    virtual void reset();

    /**
     * Intended to implemented as follows:
     *
     * Scans and tokenizes the input stream. The next possible token is returned
     * and the stream is advanced so that the entire stream can be processed by
     * subsequent calls to this method. If no further token can be generated
     * (e.g. for the input stream has reached its end), NULL shall be returned so
     * that the caller can stop calling scan().
     */
    virtual Token* scan() = 0;

protected:
    const std::map<std::string, Word>& getKeywords() const { return keywords_; }

    /**
     * Reads the next character from the stream, stores it in peek_
     * and returns it as well.
     */
    char readChar();

    char nextChar();

    /**
     * Reads the next character from the stream, compares it to the
     * character given in ref and sets the result stored in peek_, if the
     * result of the comparison is true.
     * Otherwise, if the result of the comparison is false, the method will
     * restore the previous state.
     */
    bool readCharConditioned(const char ref);
    bool readCharConditioned(bool (*pred)(const char));

    IdentifierToken* scanIdentifier(const int tokenID);
    ConstantToken* scanNumber(const int tokenID);

protected:
    std::istream* stream_;
    std::string fileName_;
    char peek_;

    typedef std::map<std::string, Word> KeywordMap;
    KeywordMap keywords_;

private:
    int numLine_;
    int numCharInLine_;
    const bool freeStream_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
