/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/utils/GLSLparser/preprocessor/pplexer.h"

#include "voreen/core/utils/GLSLparser/preprocessor/ppterminals.h"

namespace voreen {

namespace glslparser {

PreprocessorLexer::PreprocessorLexer(std::istream* const is, const std::map<std::string, Word>& keywords) 
    : Lexer(is, keywords),
    scanState_(STATE_START),
    currentDirective_()
{
}

PreprocessorLexer::PreprocessorLexer(const std::string& fileName,
                                     const std::map<std::string, Word>& keywords) 
    : Lexer(fileName, keywords),
    scanState_(STATE_START),
    currentDirective_()
{
}

PreprocessorLexer::~PreprocessorLexer() {
    for (; (currentDirective_.empty() == false); currentDirective_.pop_back())
        delete currentDirective_.back();
}

Token* PreprocessorLexer::scan() {
    if (scanState_ == STATE_END)
        return 0;

    // TODO: add correct line number and char position to returned token

    Token* token = 0;
    std::streampos stringStart = stream_->tellg();

    try {
        while (token == 0) {
            switch (scanState_) {
                case STATE_START:
                    scanState_ = STATE_SKIPPING_WHITESPACE;
                    // no break here

                case STATE_SKIPPING_WHITESPACE:
                    do {
                        readChar();
                    } while (isWhitespace(peek_));
                    break;

                case STATE_INSIDE_BLOCK_COMMENT:
                    for (int i = 0; i < 2; ) {
                        readChar();

                        switch (peek_) {
                            case '*':
                                if (i == 0)
                                    i = 1;
                                break;
                            case '/':
                                if (i == 1)
                                    i = 2;
                                break;
                            default:
                                i = 0;
                        }
                    }
                    scanState_ = STATE_START;
                    break;

                case STATE_INSIDE_LINE_COMMENT:
                    do {
                        readChar();
                    } while (! isNewline(peek_));
                    scanState_ = STATE_START;
                    break;

                case STATE_READING_DIRECTIVE:
                    // If the tokens which are held from scanning the last directive
                    // are empty, then we need to scan a new directive and keep it as
                    // the current one until all tokens from it are processed
                    //
                    if (currentDirective_.empty()) {
                        // read the directive
                        //
                        std::string directive;
                        try {
                            do {
                                if (peek_ != 0)
                                    directive += peek_;

                                readChar();

                                // if a backslash preceeds a newline, the directive does not end
                                // yet, so remove any further newlines (and whitespace)
                                //
                                if ((peek_ == '\\') && (isNewline(stream_->peek()) == true)) {
                                    do {
                                        readChar();
                                    } while (isWhitespace(peek_));
                                    directive += ' ';
                                }
                            } while (isNewline(peek_) == false);

                            currentDirective_ = scanDirective(directive);
                            currentDirective_.push_back(new Token(PreprocessorTerminals::ID_NEWLINE));

                        // Catch runtime errors and create a token from the parts of the directive
                        // which have already been read. This happens for instance, if a directive
                        // is the last occurance in the current stream.
                        //
                        } catch (std::runtime_error&) {
                            if (! directive.empty()) {
                                currentDirective_ = scanDirective(directive);
                                currentDirective_.push_back(new Token(PreprocessorTerminals::ID_NEWLINE));
                            } else {
                                // There is nothing left within the stream so prepare to
                                // reset the scanner and tell it to not emit any further
                                // text token by adjusting the start position for scanning.
                                //
                                stream_->seekg(0, std::ios_base::end);
                                stringStart = stream_->tellg();
                                scanState_ = STATE_START;
                            }
                        }
                    } // if (currentDirective_.empty()

                    // If the last scanned directive is not finished yet, return the first
                    // token from it until it is empty. Afterwards reset the scanner.
                    //
                    if (! currentDirective_.empty()) {
                        token = currentDirective_.front();
                        currentDirective_.pop_front();

                        // If the last token from the directive has been taken,
                        // start new scanning
                        //
                        if (currentDirective_.empty())
                            scanState_ = (stream_->eof()) ? STATE_END : STATE_START;
                    }
                    break;

                default:
                    readChar();
                    break;
            }   // switch (scanState)

            switch (peek_) {
                case '/':
                    if (readCharConditioned('/') == true)
                        scanState_ = STATE_INSIDE_LINE_COMMENT;
                    else if (readCharConditioned('*') == true)
                        scanState_ = STATE_INSIDE_BLOCK_COMMENT;
                    break;

                case '#':
                    // create a token containing the text read up to here
                    token = createTextToken(stringStart, stream_->tellg());
                    scanState_ = STATE_READING_DIRECTIVE;
                    break;

                case '\\':
                    if (readCharConditioned(&isNewline) == true)
                        scanState_ = STATE_SKIPPING_WHITESPACE;
                    break;
            }   // switch (peek_)
        }   // while (token == 0)

        // In case we created a runtime error by reading up to the end of the stream,
        // we have to create a TextToken for the remaining characters.
        // In addition, the scanner ist put to state STATE_END in order to indicate
        // the end of its work.
        //
    } catch (std::runtime_error&) {
        if (stream_->eof() == true) {
            scanState_ = STATE_END;
            stream_->clear();
            stream_->seekg(0, std::ios_base::end);
            if (token == 0)
                token = createTextToken(stringStart, (stream_->tellg() + std::streampos(1)));
        }
        return token;
    }

    return token;
}

// protected methods
//

std::deque<Token*> PreprocessorLexer::scanDirective(const std::string& directive) const {
    std::deque<Token*> tokenList;

    size_t stringStart = std::string::npos;
    for (size_t pos = 0; pos < directive.length(); ++pos) {
        const char top = directive[pos];
        const char next = ((pos < (directive.length() - 1)) ? directive[pos + 1] : '\0');

        if (stringStart != std::string::npos) {
            switch (top) {
                case '\\':
                    if (next == '"')
                        ++pos;
                    break;

                case '"':
                    tokenList.push_back( new StringToken(PreprocessorTerminals::ID_STRING,
                            directive.substr(stringStart, (pos - stringStart))) );
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_QUOTE));
                    stringStart = std::string::npos;
                    break;

                default:
                    break;
            } // switch (top)

            continue;
        }   // if (readingString

        // omit any kind of whitespace within the directive
        //
        if (isWhitespace(top) == true)
            continue;

        bool consumed = true;

        // at first determine characters which correspond to preprocessor operators
        //
        switch (top) {
            case '#':
                if (next == '#') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_FFENCE));
                    ++pos;
                } else if ((next == 0) || (isWhitespace(next)))
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_FENCE));
                else
                    consumed = false;
                break;

            case '(':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_LPAREN));
                break;

            case ')':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_RPAREN));
                break;

            case ',':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_COMMA));
                break;

            case '+':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_PLUS));
                break;

            case '-':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_DASH));
                break;

            case '~':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_COMPLEMENT));
                break;

            case '!':
                if (next == '=') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_NEQ));
                    ++pos;
                } else
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_NOT));
                break;

            case '*':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_MUL));
                break;

            case '/':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_DIV));
                break;

            case '%':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_MOD));
                break;

            case '<':
                if (next == '<') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_LSHIFT));
                    ++pos;
                } else if (next == '=') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_LEQ));
                    ++pos;
                } else
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_LESS));
                break;

            case '>':
                if (next == '>') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_OP_RSHIFT));
                    ++pos;
                } else if (next == '=') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_GEQ));
                    ++pos;
                } else
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_GREATER));
                break;

            case '=':
                if (next == '=') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_EQ));
                    ++pos;
                } else
                    consumed = false;
                break;

            case '&':
                if (next == '&') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_LOGICAL_AND));
                    ++pos;
                } else
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_BIT_AND));
                break;

            case '|':
                if (next == '|') {
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_LOGICAL_OR));
                    ++pos;
                } else
                    tokenList.push_back(new Token(PreprocessorTerminals::ID_BIT_OR));
                break;

            case '^':
                tokenList.push_back(new Token(PreprocessorTerminals::ID_BIT_XOR));
                break;

            case '"':
                stringStart = pos + 1;  // indicates the start of a string which will be read.
                tokenList.push_back(new Token(PreprocessorTerminals::ID_QUOTE));
                break;

            default:
                consumed = false;
                break;
        }   // switch (top)

        if (consumed == true)
            continue;

        if ((isAlpha(top) == true) || (top == '#')) {
            std::string identifier(&top, 1);
            for (++pos; pos < directive.length(); ++pos) {
                const char top = directive[pos];
                if ((isAlpha(top) == true) || (isDigit(top) == true))
                    identifier += top;
                else {
                    --pos;
                    break;
                }
            }

            KeywordMap::const_iterator it = keywords_.find(identifier);
            if (it != keywords_.end())
                tokenList.push_back(new Word(it->second));
            else
                tokenList.push_back(new IdentifierToken(PreprocessorTerminals::ID_IDENTIFIER, identifier));
        } else if ((isDigit(top) == true) || ((top == '.') && (isDigit(next) == true))) {
            size_t end = 0;
            ConstantToken* ct = scanNumber(directive.substr(pos), end);
            if (ct != 0)
                tokenList.push_back(ct);
            pos += end;
        } else
            tokenList.push_back(new TextToken(PreprocessorTerminals::ID_TEXT, std::string(&top, 1)));
    }   // for (pos < length)

    return tokenList;
}

ConstantToken* PreprocessorLexer::scanNumber(const std::string& input, size_t& end) const
{
    std::string numberString("");
    ConstantToken::Base base = ConstantToken::BASE_10;
    ConstantToken::ConstantType type = ConstantToken::TYPE_INT;
    bool hasExponent = false;
    bool (*pred)(const char ref) = &isDigit;

    size_t pos = 0;
    for (bool atEnd = false; ((! atEnd) && (pos < input.length())); ++pos) {
        const char top = input[pos];
        const char next = ((pos < (input.length() - 1)) ? input[pos + 1] : '\0');

        if ((pred != 0) && ((*pred)(top) == true)) {
            numberString += top;
            if ((pos == 0) && (top == '0')) {
                if ((next == 'x') || (next == 'X')) {
                    numberString += next;
                    ++pos;
                    pred = &isHexDigit;
                    base = ConstantToken::BASE_16;
                } else if (next == '.')
                    type = ConstantToken::TYPE_FLOAT;
                else {
                    pred = &isOctalDigit;
                    base = ConstantToken::BASE_8;
                }
            }

            switch (next) {
                case 'u':
                case 'U':
                    atEnd = true;
                    if (type == ConstantToken::TYPE_INT) {
                        numberString += next;
                        ++pos;
                        type = ConstantToken::TYPE_UINT;
                    }
                    break;

                case 'f':
                case 'F':
                    atEnd = true;
                    if (type == ConstantToken::TYPE_FLOAT) {
                        numberString += next;
                        ++pos;
                    }
                    break;

                case '.':
                    if (base == ConstantToken::BASE_10) {
                        type = ConstantToken::TYPE_FLOAT;
                        numberString += next;
                        ++pos;
                    }
                    break;

                case 'e':
                case 'E':
                    if (type == ConstantToken::TYPE_FLOAT) {
                        hasExponent = true;
                        numberString += next;
                        ++pos;
                    }
                    break;

                default:
                    break;
            }
        }  else if ((top == 'f') && (type == ConstantToken::TYPE_FLOAT)) {
            numberString += top;
        } else if ((top == '.') && (base == ConstantToken::BASE_10)) {
            // Accept the dot '.' only, if it is not the first character or
            // if a digit is following to avoid invalid ".f" strings.
            //
            if ((pos > 0) || (isDigit(next) == true)) {
                type = ConstantToken::TYPE_FLOAT;
                numberString += top;
            }
        } else if (((top == '+') || (top == '-')) && (hasExponent == true)) {
            numberString += top;
        } else
            break;
    }   // for (pos

    end = (pos > 0) ? --pos : 0;
    if (numberString.empty() == true)
        return 0;

    return new ConstantToken(PreprocessorTerminals::ID_CONSTANT, numberString, type, base);
}

Token* PreprocessorLexer::createTextToken(const std::streampos& start,
                                          const std::streampos& end)
{
    // copy the number of characters to the buffer
    //
    std::streamsize count = (end - std::streampos(1) - start);
    char* buffer = 0;
    if (count > 0) {
        buffer = new char[static_cast<size_t>(count) + 1];
        stream_->seekg(start, std::ios_base::beg);
        stream_->get(buffer, count + 1, char(255)); // char(255) prevents get() to stop at newlines.
        stream_->seekg(end, std::ios_base::beg);
        buffer[count] = 0;
    } else
        return 0;

    if (buffer != 0) {
        // Return a special 'newline' token, if the text only contains newlines.
        // Those tokens are required to terminate preprocessor directives
        //
        bool newlinesOnly = true;
        for (int i = 0; ((newlinesOnly) && (i < static_cast<int>(count))); ++i) {
            if (buffer[i] != 0)
                newlinesOnly = isNewline(buffer[i]);
        }

        Token* token = 0;
        if (! newlinesOnly)
            token = new TextToken(PreprocessorTerminals::ID_TEXT, std::string(buffer, static_cast<size_t>(count)));

        delete [] buffer;
        return token;
    }

    return new TextToken(PreprocessorTerminals::ID_TEXT, "");
}

}   // namespace glslparser

}   // namespace voreen
