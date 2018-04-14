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

#ifndef VRN_GLSLVARIABLE_H
#define VRN_GLSLVARIABLE_H

#include "voreen/core/utils/GLSLparser/glsl/glslexpression.h"

namespace voreen {

namespace glslparser {

class GLSLVariable : public GLSLNode {
public:
    GLSLVariable(IdentifierToken* const token, const bool isArray,
        GLSLConditionalExpression* const numArrayElements,
        GLSLExpression* const initializer = 0)
        : GLSLNode(token->getTokenID()),
        name_(token->getValue()),
        isArray_(isArray),
        numArrayElements_(numArrayElements),
        initializer_(initializer)
    {
    }

    ~GLSLVariable() {
        delete numArrayElements_;
        delete initializer_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_VARIABLE; }

    GLSLExpression* getInitializer() const { return initializer_; }
    const std::string& getName() const { return name_; }
    GLSLConditionalExpression* getNumArrayElements() const { return numArrayElements_; }
    bool isArray() const { return isArray_; }

public:
    std::string name_;
    bool isArray_;
    GLSLConditionalExpression* const numArrayElements_;
    GLSLExpression* const initializer_;
};

}   // namespace glslparser

}   // namespace

#endif
