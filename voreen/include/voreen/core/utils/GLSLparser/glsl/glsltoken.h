/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_GLSLTOKEN_H
#define VRN_GLSLTOKEN_H

#include "voreen/core/utils/GLSLparser/token.h"

#include <list>

namespace voreen {

namespace glslparser {

typedef GenericToken<std::string> TextToken;

// ============================================================================

class AnnotationToken : public Token {
public:
    AnnotationToken(const int symbolID, const std::list<Token*>& parts)
        : Token(symbolID),
        parts_(parts)
    {
    }

    AnnotationToken(const AnnotationToken& other)
        : Token(other.getTokenID())
    {
        for (std::list<Token*>::const_iterator it = other.parts_.begin();
            it != other.parts_.end(); ++it)
        {
            parts_.push_back((*it)->getCopy());
        }
    }

    virtual Token* getCopy() const { return new AnnotationToken(*this); }

    const std::list<Token*>& getParts() const { return parts_; }

protected:
    std::list<Token*> parts_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
