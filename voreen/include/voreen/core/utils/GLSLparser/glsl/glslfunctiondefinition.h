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

#ifndef VRN_GLSLFUNCTIONDEFINITION_H
#define VRN_GLSLFUNCTIONDEFINITION_H

#include "voreen/core/utils/GLSLparser/glsl/glslexternaldeclaration.h"
#include "voreen/core/utils/GLSLparser/glsl/glslfunctionprototype.h"
#include "voreen/core/utils/GLSLparser/glsl/glslstatement.h"

namespace voreen {

namespace glslparser {

class GLSLFunctionDefinition : public GLSLExternalDeclaration
{
public:
    /**
     * @param   funcProto   Must not be null.
     * @param   statements  Must not be null and newScope must be true.
     */
    GLSLFunctionDefinition(GLSLFunctionPrototype* const funcProto,
        GLSLCompoundStatement* const statements)
        : GLSLExternalDeclaration(),
        function_(funcProto),
        statements_(statements)
    {
    }

    virtual ~GLSLFunctionDefinition() {
        delete function_;
        delete statements_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_FUNCTION_DEFINITION; }

    GLSLFunctionPrototype* getFunctionPrototype() { return function_; }
    GLSLCompoundStatement* getStatements() { return statements_; }

protected:
    GLSLFunctionPrototype* const function_;
    GLSLCompoundStatement* const statements_;
};  // class GLSLFunctionDefinition

}   // namespace glslparser

}   // namespace voreen

#endif  // VRN_GLSLFUNCTIONDEFINITION_H

