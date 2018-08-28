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

#ifndef VRN_GLSLPARAMETER_H
#define VRN_GLSLPARAMETER_H

#include "voreen/core/utils/GLSLparser/glsl/glsltypespecifier.h"
#include "voreen/core/utils/GLSLparser/glsl/glslvariable.h"

namespace voreen {

namespace glslparser {

class GLSLParameter : public GLSLNode {
public:
    GLSLParameter(GLSLTypeSpecifier* const typeSpecifier,
                  GLSLParameterQualifier* const paramQualifier,
                  const bool isConst, GLSLVariable* const paramName)
        : GLSLNode(GLSLTerminals::ID_UNKNOWN),
        typeSpecifier_(typeSpecifier),
        paramQualifier_(paramQualifier),
        isConst_(isConst),
        name_(paramName)
    {
    }

    virtual ~GLSLParameter() {
        delete typeSpecifier_;
        delete paramQualifier_;
        delete name_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_PARAMETER; }

    void setParameterQulifier(GLSLParameterQualifier* const paramQualifier) {
        paramQualifier_ = paramQualifier;
    }

    void setIsConst(const bool isConst) { isConst_ = isConst; }

    GLSLTypeSpecifier* getSpecifier() { return typeSpecifier_; }
    GLSLParameterQualifier* getQualifier() { return paramQualifier_; }
    bool isConst() const { return isConst_; }
    GLSLVariable* getName() { return name_; }

protected:
    GLSLTypeSpecifier* typeSpecifier_;
    GLSLParameterQualifier* paramQualifier_;
    bool isConst_;
    GLSLVariable* name_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
