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

#ifndef VRN_GLSLFUNCTIONPROTOTYPE_H
#define VRN_GLSLFUNCTIONPROTOTYPE_H

#include "voreen/core/utils/GLSLparser/glsl/glslparameter.h"

namespace voreen
{
namespace glslparser
{

class GLSLFunctionPrototype : public GLSLNode
{
public:
    /**
     * @param   name    Must not be NULL.
     * @param   typeSpecifier   Must not be NULL.
     */
    GLSLFunctionPrototype(IdentifierToken* const name, GLSLTypeSpecifier* const typeSpecifier)
        : GLSLNode(name->getTokenID()),
        name_(dynamic_cast<IdentifierToken* const>(name->getCopy())),
        typeSpecifier_(typeSpecifier)
    {
    }

    virtual ~GLSLFunctionPrototype() {
        delete name_;
        delete typeSpecifier_;

        for (size_t i = 0; i < params_.size(); ++i)
            delete params_[i];
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_FUNCTION_PROTOTYPE; }

    void addParameter(GLSLParameter* const param) {
        if (param != 0)
            params_.push_back(param);
    }

    IdentifierToken* getName() { return name_; }
    GLSLTypeSpecifier* getSpecifier() { return typeSpecifier_; }
    const std::vector<GLSLParameter*>& getParameters() const { return params_; }

protected:
    IdentifierToken* const name_;
    GLSLTypeSpecifier* const typeSpecifier_;
    std::vector<GLSLParameter*> params_;
};  // class GLSLFunctionPrototype

}   // namespace glslparser
}   // namespace voreen

#endif  // VRN_GLSLFUNCTIONPROTOTYPE_H
