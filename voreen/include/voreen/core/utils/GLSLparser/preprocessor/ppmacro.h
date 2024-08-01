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

#ifndef VRN_PPMACRO_H
#define VRN_PPMACRO_H

#include "voreen/core/utils/GLSLparser/preprocessor/ppexpression.h"

namespace voreen {

namespace glslparser {

class Macro : public Expression {
public:
    Macro(IdentifierToken* const token, const bool isFunction)
        : Expression(token),
        identifier_(token->getValue()),
        isFunction_(isFunction),
        parameters_(0)
    {
    }

    Macro(IdentifierToken* const token, const bool isFunction, ExpressionList* const params)
        : Expression(token),
        identifier_(token->getValue()),
        isFunction_(isFunction),
        parameters_(params)
    {
    }

    virtual ~Macro() {
        delete parameters_;
    }

    int getNodeType() const { return PreprocessorNodeTypes::NODE_MACRO; }

    const std::string& getIdentifier() const { return identifier_; }
    bool isFunction() const { return isFunction_; }
    ExpressionList* getParameters() const { return parameters_; }
    size_t getNumParameters() const { return parameters_->size(); }

    virtual void tokenize(TokenList* const list) const;

    std::list<Token*> expandMacro(PreprocessorSymbol* const symbol);

protected:
    std::string identifier_;
    bool isFunction_;
    ExpressionList* const parameters_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
