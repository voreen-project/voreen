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

#ifndef VRN_ANNOTATION_TERMINALS_H
#define VRN_ANNOTATION_TERMINALS_H

#include "voreen/core/utils/GLSLparser/terminalsmap.h"

namespace voreen {

namespace glslparser {

class AnnotationTerminals : public TerminalsMap {
public:
    enum SymbolIDs {
        ID_START, ID_LBRACKET, ID_RBRACKET, ID_EQUAL, ID_AT,
        ID_COMMA, ID_QUOTE, ID_APOSTROPHE,

        ID_TRUE, ID_FALSE, ID_IDENTIFIER, ID_INTCONST,
        ID_FLOATCONST, ID_STRING, ID_CHARACTER, ID_UNKNOWN,
        ID_END, ID_EMPTY
    };

public:
    AnnotationTerminals();

};

}   // namespace glslparser

}   // namespace voreen

#endif
