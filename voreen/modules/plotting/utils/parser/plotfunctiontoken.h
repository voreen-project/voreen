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

#ifndef VRN_PLOTFUNCTIONTOKEN_H
#define VRN_PLOTFUNCTIONTOKEN_H


#include "voreen/core/utils/GLSLparser/token.h"

namespace voreen {

namespace glslparser {

//----------  FunctionToken ---------------------------------------------------

/**
 * \brief represents the function as a token
 */
class FunctionToken : public IdentifierToken {
public:

    FunctionToken(const int tokenID, const std::string& value, const int index);

    /// Destructor
    virtual ~FunctionToken() {};
    // create a new Token with the same properties
    virtual Token* getCopy() const;
    // return the indexvalue
    int getIndex();


private:
    int index_;

};

//----------  VariablesToken ---------------------------------------------------

/**
 * \brief represents the variable as a token
 */
class VariablesToken : public IdentifierToken {
public:

    VariablesToken(const int tokenID, const std::string& value);

    /// Destructor
    virtual ~VariablesToken() {};
    // create a new Token with the same properties
    virtual Token* getCopy() const;

};

//----------  OperatorToken ---------------------------------------------------
/**
 * \brief represents the operator as a token
 */
class OperatorToken : public GenericToken<char> {
public:
    OperatorToken(const int tokenID, const char value, const int parameter = -1);

    /// Destructor
    virtual ~OperatorToken() {};
    // create a new Token with the same properties
    virtual Token* getCopy() const;
    // set the parametercount of this operator (e.g. + can be 2[a+b] or 1 [+1]
    void setParameter(int parameter);
    // gives back the parametercount
    int getParameter();


private:
    int parameter_;
};

//----------  Further Tokens ---------------------------------------------------
/**
 * \brief represents the Isolator as a token (",", ";","|")
 */
typedef GenericToken<char> IsolatorToken;

/**
 * \brief represents the bracket as a token
 */
typedef GenericToken<char> BracketToken;

}// namespace glslparser
}// namespace voreen


#endif // VRN_EXPRESSION_H
