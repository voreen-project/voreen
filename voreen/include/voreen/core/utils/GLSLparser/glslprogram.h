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

#ifndef VRN_GLSLPROGRAM_H
#define VRN_GLSLPROGRAM_H

#include "voreen/core/utils/GLSLparser/glsl/glslsymbol.h"

#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace voreen {

namespace glslparser {

class GLSLProgram {
public:
    GLSLProgram(std::istream* is);
    GLSLProgram(const std::string& fileName);

    ~GLSLProgram();

    void freeDeclarations();

    const std::ostringstream& getLog() const { return log_; }
    const std::string& getShaderHeader() const { return shaderHeader_; }
    const std::vector<GLSLVariableSymbol*>& getUniformDeclarations() const { return uniformDecls_; }
    const std::vector<GLSLVariableSymbol*>& getOutDeclarations() const { return outDecls_; }
    const std::set<unsigned int>& getReferencedGlFragData() const { return referencedGlFragData_; }

    void setShaderHeader(const std::string& shaderHeader) { shaderHeader_ = shaderHeader; }

    bool parse();

private:
    std::string fileName_;
    std::istream* const is_;
    std::string shaderHeader_;
    std::vector<GLSLVariableSymbol*> uniformDecls_;
    std::vector<GLSLVariableSymbol*> outDecls_;
    std::set<unsigned int> referencedGlFragData_;
    std::ostringstream log_;

    static const std::string loggerCat_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
