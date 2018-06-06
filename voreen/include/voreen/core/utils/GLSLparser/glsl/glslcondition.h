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

#ifndef VRN_GLSLCONDITION_H
#define VRN_GLSLCONDITION_H

#include "voreen/core/utils/GLSLparser/glsl/glsltypespecifier.h"

namespace voreen {

namespace glslparser {

class GLSLCondition : public GLSLNode
{
public:
    GLSLCondition(GLSLExpression* const expr)
        : GLSLNode(GLSLTerminals::ID_UNKNOWN),
        expr_(expr),
        spec_(0),
        identifier_(0)
    {
    }

    GLSLCondition(GLSLTypeSpecifier* const spec, IdentifierToken* const identifier,
        GLSLExpression* const init)
        : GLSLNode(identifier->getTokenID()),
        expr_(init),
        spec_(spec),
        identifier_(identifier)
    {
    }

    virtual ~GLSLCondition() {
        delete expr_;
        delete spec_;
        delete identifier_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_CONDITION; }

    GLSLExpression* getExpression() { return expr_; }
    IdentifierToken* getIdentifier() { return identifier_; }
    GLSLTypeSpecifier* getTypeSpecifier() { return spec_; }

protected:
    GLSLExpression* const expr_;
    GLSLTypeSpecifier* const spec_;
    IdentifierToken* const identifier_;
};  // GLSLCondition

}   // namespace glslparser

}   // namespace voreen

#endif  // VRN_GLSLCONDITION_H

