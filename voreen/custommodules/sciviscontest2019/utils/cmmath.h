
/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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

#ifndef VRN_CMMATH_H
#define VRN_CMMATH_H

#include "tgt/matrix.h"
#include "tgt/bounds.h"
namespace voreen{
static tgt::vec3 CMtransform(tgt::mat4 M, tgt::vec3 v){
    tgt::vec4 w = M*tgt::vec4(v, 1.0f);
    return w.xyz();
}

static tgt::mat4 CMinvert(tgt::mat4 M){
    tgt::mat4 O;
    if (!M.invert(O))
        tgtAssert(false, "Could not invert matrix");
    return O;
}

static tgt::Bounds CMtransformBounds(tgt::mat4 M, tgt::Bounds b){
    return tgt::Bounds(CMtransform(M, b.getLLF()), CMtransform(M, b.getURB()));
}
}

#endif