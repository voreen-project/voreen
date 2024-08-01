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

#include "voreen/core/utils/GLSLparser/glsl/glslsymbol.h"

namespace voreen {

namespace glslparser {

GLSLSymbol::GLSLSymbol(const std::string& identifier, const InternalType type, const int elementCount,
                       const bool isDeclared)
    : Symbol(identifier),
    internalType_(type),
    precision_(PRECQ_NONE),
    elementCount_(elementCount),
    isDeclared_(isDeclared)
{
}

GLSLSymbol::~GLSLSymbol() {
    for (size_t i = 0; i < annotations_.size(); ++i)
        delete annotations_[i];
}

const std::string GLSLSymbol::getAnnotationValueString(std::string name) const {
    for (unsigned int i=0; i<annotations_.size(); i++)
        if (annotations_[i]->getName() == name)
            return annotations_[i]->getValueString();
    return "";
}

const GLSLAnnotation* GLSLSymbol::getAnnotation(std::string name) const {
    for (unsigned int i=0; i<annotations_.size(); i++)
        if (annotations_[i]->getName() == name)
            return annotations_[i];
    return 0;
}

// ============================================================================

GLSLVariableSymbol::GLSLVariableSymbol(const std::string& identifier, const InternalType type,
        const int elementCount, const bool isDeclared)
        : GLSLSymbol(identifier, type, elementCount, isDeclared),
        storage_(SQ_NONE),
        interpolation_(IQ_NONE),
        hasInvariantQualifier_(false),
        isArray_(false),
        numArrayElements_(0)
{
}

}   // namespace glslparser

}   // namespace voreen
