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

#include "voreen/core/datastructures/roi/roicylinder.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

using tgt::vec2;
using tgt::vec3;
using tgt::vec4;
using tgt::mat4;

const std::string ROICylinder::loggerCat_ = "voreen.ROICylinder";

ROICylinder::ROICylinder() : ROINormalizedGeometry()
{
}

ROICylinder::ROICylinder(Grid grid, tgt::vec3 center, tgt::vec3 dimensions, int primaryDir) : ROINormalizedGeometry(grid, center, dimensions)
{
    if((primaryDir < 0) || (primaryDir > 2)) {
        primaryDir = 0;
        LERROR("Invalid primary direction!");
    }
    setPrimaryDir(primaryDir);
}

ROICylinder::ROICylinder(Grid grid, tgt::vec3 center, float radius, float length, int primaryDir)  : ROINormalizedGeometry(grid, center, vec3(radius*2.0f)) {
    if((primaryDir < 0) || (primaryDir > 2)) {
        primaryDir = 0;
        LERROR("Invalid primary direction!");
    }
    setPrimaryDir(primaryDir);

    vec3 d = vec3(radius*2.0f);
    d[primaryDir] = length;
    setDimensions(d);
}

bool ROICylinder::inROINormalized(tgt::vec3 p) const {
    switch(getPrimaryDir()) {
        case 0: return ((std::abs(p.x) < 1.0f) && (distance(p.yz(), vec2(0.0f)) <= 1.0f));
        case 1: return ((std::abs(p.y) < 1.0f) && (distance(vec2(p.x, p.z), vec2(0.0f)) <= 1.0f));
        case 2: return ((std::abs(p.z) < 1.0f) && (distance(p.xy(), vec2(0.0f)) <= 1.0f));
    }
    LWARNING("Primary direction invalid!");
    return false;
}

Geometry* ROICylinder::generateNormalizedMesh() const {
    TriangleMeshGeometrySimple* geometry = new TriangleMeshGeometrySimple();
    vec3 e1(0.0f);
    vec3 e2(0.0f);
    vec3 e3(0.0f);

    if(getPrimaryDir() == 0) {
        e1 = vec3(0.0f, 1.0f, 0.0f);
        e2 = vec3(0.0f, 0.0f, 1.0f);
        e3 = vec3(1.0f, 0.0f, 0.0f);
    }
    else if(getPrimaryDir() == 1) {
        e1 = vec3(1.0f, 0.0f, 0.0f);
        e2 = vec3(0.0f, 0.0f, 1.0f);
        e3 = vec3(0.0f, 1.0f, 0.0f);
    }
    else if(getPrimaryDir() == 2) {
        e1 = vec3(1.0f, 0.0f, 0.0f);
        e2 = vec3(0.0f, 1.0f, 0.0f);
        e3 = vec3(0.0f, 0.0f, 1.0f);
    }

    int aSteps = 20;
    int lSteps = std::max(1, tgt::ifloor(getDimensions()[getPrimaryDir()] / getGrid().getSpacing()[getPrimaryDir()]));
    for(int i=0; i<aSteps; i++) {
        for(int j=0; j<lSteps; j++) {
            float angleA1 = ((float)i / (float)aSteps) * 2.0f * tgt::PIf;
            float angleA2 = ((float)(i+1) / (float)aSteps) * 2.0f * tgt::PIf;
            float l1 = (2.0f * (float)j / (float) lSteps) - 1.0f;
            float l2 = (2.0f * (float)(j+1) / (float) lSteps) - 1.0f;

            geometry->addQuad(VertexBase(e1*cosf(angleA1) + e2*sinf(angleA1) + e3*l1),
                              VertexBase(e1*cosf(angleA1) + e2*sinf(angleA1) + e3*l2),
                              VertexBase(e1*cosf(angleA2) + e2*sinf(angleA2) + e3*l2),
                              VertexBase(e1*cosf(angleA2) + e2*sinf(angleA2) + e3*l1));
        }
    }

    return geometry;
}

Geometry* ROICylinder::generateNormalizedMesh(tgt::plane pl) const {
    // We do this by clipping a quad
    TriangleMeshGeometrySimple* planeFace = createQuad(pl);

    // First against the caps:
    if(getPrimaryDir() == 0) {
        tgt::plane xm(-1.0f, 0.0f, 0.0f, 1.0f);
        tgt::plane xp(1.0f, 0.0f, 0.0f, 1.0f);
        planeFace->clip(xm);
        planeFace->clip(xp);
    }
    else if(getPrimaryDir() == 1) {
        tgt::plane ym(0.0f, -1.0f, 0.0f, 1.0f);
        tgt::plane yp(0.0f, 1.0f, 0.0f, 1.0f);
        planeFace->clip(ym);
        planeFace->clip(yp);
    }
    else if(getPrimaryDir() == 2) {
        tgt::plane zm(0.0f, 0.0f, -1.0f, 1.0f);
        tgt::plane zp(0.0f, 0.0f, 1.0f, 1.0f);
        planeFace->clip(zm);
        planeFace->clip(zp);
    }

    // Then against 20 tangential planes that describe the cylinder:
    vec3 e1(0.0f);
    vec3 e2(0.0f);
    vec3 e3(0.0f);

    if(getPrimaryDir() == 0) {
        e1 = vec3(0.0f, 1.0f, 0.0f);
        e2 = vec3(0.0f, 0.0f, 1.0f);
        e3 = vec3(1.0f, 0.0f, 0.0f);
    }
    else if(getPrimaryDir() == 1) {
        e1 = vec3(1.0f, 0.0f, 0.0f);
        e2 = vec3(0.0f, 0.0f, 1.0f);
        e3 = vec3(0.0f, 1.0f, 0.0f);
    }
    else if(getPrimaryDir() == 2) {
        e1 = vec3(1.0f, 0.0f, 0.0f);
        e2 = vec3(0.0f, 1.0f, 0.0f);
        e3 = vec3(0.0f, 0.0f, 1.0f);
    }

    int steps = 20;
    for(int i=0; i<steps; i++) {
        float angle = ((float)i / (float)steps) * 2.0f * tgt::PIf;
        //face.addVertex(VertexGeometry( + e3*l1, vec3(0.0f), getColor()));

        tgt::plane pl(e1*cosf(angle) + e2*sinf(angle), 1.0f);
        planeFace->clip(pl);
    }

    return planeFace;
}

} // namespace
