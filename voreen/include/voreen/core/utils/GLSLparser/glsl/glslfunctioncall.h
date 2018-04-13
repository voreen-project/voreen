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

#ifndef VRN_GLSLFUNCTIONCALL_H
#define VRN_GLSLFUNCTIONCALL_H

#include "voreen/core/utils/GLSLparser/glsl/glsltypespecifier.h"

namespace voreen {

namespace glslparser {

class GLSLFunctionCall : public GLSLExpression {
public:
    GLSLFunctionCall(IdentifierToken* const name)
        : GLSLExpression(*name),
        name_(dynamic_cast<IdentifierToken* const>(name->getCopy())),
        typeSpecifier_(0),
        postfixExpression_(0)
    {
    }

    GLSLFunctionCall(GLSLTypeSpecifier* const typeSpec)
        : GLSLExpression(Token(typeSpec->getSymbolID())),
        name_(0),
        typeSpecifier_(typeSpec),
        postfixExpression_(0)
    {
    }

    virtual ~GLSLFunctionCall() {
        delete name_;
        delete typeSpecifier_;
        delete postfixExpression_;

        for (size_t i = 0; i < parameters_.size(); ++i)
            delete parameters_[i];
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_FUNCTION_CALL; }

    void addParameter(GLSLExpression* const param) {
        if (param != 0)
            parameters_.push_back(param);
    }

    void setPostfixExpression(GLSLExpression* const postfixExpr) {
        postfixExpression_ = postfixExpr;
    }

protected:
    IdentifierToken* const name_;   // either one of them is used at a time only
    GLSLTypeSpecifier* const typeSpecifier_;
    GLSLExpression* postfixExpression_;
    std::vector<GLSLExpression*> parameters_;
};

}   // namespace glslparser

}   // namespace

#endif
