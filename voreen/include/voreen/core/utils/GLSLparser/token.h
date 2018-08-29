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

#ifndef VRN_TOKEN_H
#define VRN_TOKEN_H

#include <string>
#include <sstream>

namespace voreen {

namespace glslparser {

class Token {
public:
    Token(const int tokenID);
    virtual ~Token() {}

    bool operator==(const Token& rhs) const;
    bool operator!=(const Token& rhs) const;

    virtual Token* getCopy() const { return new Token(*this); }

    int getTokenID() const { return tokenID_; }
    int getLineNumber() const { return numLine_; }
    int getCharNumber() const { return numCharInLine_; }

    void replaceTokenID(const int newID) { tokenID_ = newID; }
    void setLineNumber(const int numLine) { numLine_ = numLine; }
    void setCharNumber(const int numCharInLine) { numCharInLine_ = numCharInLine; }

    virtual std::string toString() const;

private:
    int tokenID_;
    int numLine_;
    int numCharInLine_;
};  // class Token

// ============================================================================

template<typename T>
class GenericToken : public Token {
public:
    GenericToken(const int tokenID, const T value)
        : Token(tokenID), value_(value)
    {
    }

    GenericToken(const GenericToken<T>& rhs)
        : Token(rhs.getTokenID()), value_(rhs.value_)
    {
    }

    virtual ~GenericToken() {}

    virtual Token* getCopy() const { return new GenericToken<T>(*this); }

    const T& getValue() const { return value_; }
    virtual std::string toString() const;

protected:
    T value_;
};   // class GenericToken

template<typename T>
std::string GenericToken<T>::toString() const {
    std::ostringstream oss;
    oss << "<" << value_ << ", ID = " << getTokenID() << ">";
    return oss.str();
}

// ============================================================================

typedef GenericToken<std::string> Word;
typedef GenericToken<std::string> IdentifierToken;
typedef GenericToken<std::string> StringToken;



// ============================================================================

class ConstantToken : public GenericToken<std::string> {
public:
    enum ConstantType { TYPE_INT, TYPE_UINT, TYPE_FLOAT };
    enum Base { BASE_8, BASE_10, BASE_16};

public:
    ConstantToken(const int tokenID, const std::string& strValue,
        const ConstantToken::ConstantType type, const ConstantToken::Base base = BASE_10)
        : GenericToken<std::string>(tokenID, strValue),
        type_(type),
        base_(base)
    {
    }

    virtual Token* getCopy() const { return new ConstantToken(*this); }

    template<typename T>
    T convert() const {
        T c = T(0);
        std::istringstream iss(value_);
        iss >> c;
        return c;
    }

    Base getBase() const { return base_; }
    ConstantType getType() const { return type_; }

private:
    ConstantType type_;
    Base base_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
