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

#include "voreen/core/datastructures/roi/roicube.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

using tgt::vec3;
using tgt::vec4;
using tgt::mat4;

const std::string ROICube::loggerCat_ = "voreen.ROICube";

ROICube::ROICube() : ROINormalizedGeometry()
{
}

ROICube::ROICube(Grid grid, tgt::vec3 center, tgt::vec3 dimensions) : ROINormalizedGeometry(grid, center, dimensions)
{
}

bool ROICube::inROINormalized(tgt::vec3 p) const {
    if(tgt::hand(tgt::greaterThanEqual(p, vec3(-1.0f))) && tgt::hand(tgt::lessThanEqual(p, vec3(1.0f))))
        return true;
    else
        return false;
}

Geometry* ROICube::generateNormalizedMesh() const {
    return TriangleMeshGeometrySimple::createCube(vec3(-1.0f), vec3(1.0f));
}

Geometry* ROICube::generateNormalizedMesh(tgt::plane pl) const {
    TriangleMeshGeometrySimple* planeFace = createQuad(pl);

    tgt::plane xm(-1.0f, 0.0f, 0.0f, 1.0f);
    tgt::plane xp(1.0f, 0.0f, 0.0f, 1.0f);
    tgt::plane ym(0.0f, -1.0f, 0.0f, 1.0f);
    tgt::plane yp(0.0f, 1.0f, 0.0f, 1.0f);
    tgt::plane zm(0.0f, 0.0f, -1.0f, 1.0f);
    tgt::plane zp(0.0f, 0.0f, 1.0f, 1.0f);

    planeFace->clip(xm);
    planeFace->clip(xp);
    planeFace->clip(ym);
    planeFace->clip(yp);
    planeFace->clip(zm);
    planeFace->clip(zp);

    return planeFace;
}

} // namespace
