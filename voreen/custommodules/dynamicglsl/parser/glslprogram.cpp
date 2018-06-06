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

#include "glslprogram.h"

#include "glslparser.h"

#include "voreen/core/utils/GLSLparser/glsl/glslvisitor.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppparser.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppvisitor.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

namespace glslparser {

const std::string GLSLProgram::loggerCat_("voreen.dynamicglsl.GLSLProgram");

GLSLProgram::GLSLProgram(std::istream* is)
    : fileName_(""),
    is_(is),
    shaderHeader_(""),
    log_(std::ios_base::out | std::ios_base::binary)
{
}

GLSLProgram::GLSLProgram(const std::string& fileName)
    : fileName_(fileName),
    is_(0),
    shaderHeader_(""),
    log_(std::ios_base::out | std::ios_base::binary)
{
}

GLSLProgram::~GLSLProgram() {
    freeDeclarations();
}

void GLSLProgram::freeDeclarations() {
    for (size_t i = 0; i < uniformDecls_.size(); ++i)
        delete uniformDecls_[i];

    uniformDecls_.clear();
}

// private methods
//

bool GLSLProgram::parse() {
    freeDeclarations();

    PreprocessorVisitor preprocessor;
    preprocessor.setShaderHeader(shaderHeader_);
    std::ostringstream& program = preprocessor.translate(is_, VoreenApplication::app()->getBasePath()+"/glsl");

    log_ << "parsing GLSL program in file '" << fileName_ << "'...\n\n";
    log_ << "Preprocessor:\n================\n";
    log_ << preprocessor.getLog().str() << "\n\n";

    /*
    std::ofstream pre("translated.txt", std::ios_base::binary | std::ios_base::out);
    pre << program.str();
    pre.close();
    */

    std::istream prog(program.rdbuf());
    GLSLParser glslParser(&prog);
    //glslParser.setDebugging(true);
    ParseTreeNode* root = glslParser.parse();
    log_ << "GLSL parser:\n================\n";
    log_ << glslParser.getLog().str() << "\n\n";

    bool res = false;
    if (root != 0) {
        res = true;
        GLSLVisitor glsl;
        glsl.visitAll(root);

        uniformDecls_ = glsl.getUniforms(false);
        outDecls_ = glsl.getOuts(false);
        referencedGlFragData_ = glsl.getReferencedGlFragData();

        if (glsl.getNumWarnings() > 0)
            std::cout << glsl.getLog().str() << std::endl;

    } else
        LINFO(glslParser.getLog().str());
    delete root;

    return res;
}

}   // namespace glslparser

}   // namespace voreen
