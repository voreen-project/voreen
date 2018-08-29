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

#include "voreen/core/datastructures/roi/roisphere.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

using tgt::vec3;
using tgt::mat4;

const std::string ROISphere::loggerCat_ = "voreen.ROISphere";

ROISphere::ROISphere() : ROINormalizedGeometry()
{
}

ROISphere::ROISphere(Grid grid, tgt::vec3 center, tgt::vec3 dimensions) : ROINormalizedGeometry(grid, center, dimensions)
{
}

ROISphere::ROISphere(Grid grid, tgt::vec3 center, float radius)  : ROINormalizedGeometry(grid, center, vec3(radius*2.0f)) {
}

bool ROISphere::inROINormalized(tgt::vec3 p) const {
    if(distance(p, vec3(0.0f)) <= 1.0f)
        return true;
    else
        return false;
}

Geometry* ROISphere::generateNormalizedMesh() const {
    TriangleMeshGeometrySimple* geometry = new TriangleMeshGeometrySimple();

    int steps = 20;
    for(int i=0; i<steps; i++) {
        for(int j=0; j<steps; j++) {
            float angleA1 = ((float)i / (float)steps) * 2.0f * tgt::PIf;
            float angleA2 = ((float)(i+1) / (float)steps) * 2.0f * tgt::PIf;
            float angleB1 = ((float)j / (float)steps) * tgt::PIf;
            float angleB2 = ((float)(j+1) / (float)steps) * tgt::PIf;

            geometry->addQuad(VertexBase(vec3(sin(angleB1)*cos(angleA1), sin(angleB1)*sin(angleA1), cosf(angleB1))),
                              VertexBase(vec3(sin(angleB2)*cos(angleA1), sin(angleB2)*sin(angleA1), cosf(angleB2))),
                              VertexBase(vec3(sin(angleB2)*cos(angleA2), sin(angleB2)*sin(angleA2), cosf(angleB2))),
                              VertexBase(vec3(sin(angleB1)*cos(angleA2), sin(angleB1)*sin(angleA2), cosf(angleB1))));
        }
    }

    return geometry;
}

Geometry* ROISphere::generateNormalizedMesh(tgt::plane pl) const {
    TriangleMeshGeometrySimple* mlg = new TriangleMeshGeometrySimple();

    vec3 temp(1.0f, 1.0f, 0.0f);
    if(dot(temp, pl.n) > 0.9f)
        temp = vec3(0.0f, 1.0f, 1.0f);

    if(fabs(pl.d) > 1.0f)
        return mlg;

    vec3 t1 = normalize(cross(pl.n, temp));
    vec3 t2 = normalize(cross(pl.n, t1));
    float circleRadius = sqrt(1.0f - pl.d*pl.d);
    t1 *= circleRadius;
    t2 *= circleRadius;
    vec3 base = -(pl.d * pl.n);

    int steps = 40;
    std::vector<VertexBase> vertices;
    for(int i=0; i<steps; i++) {
        float angle = ((float)i / (float)steps) * 2.0f * tgt::PIf;
        vertices.push_back(VertexBase(base + sinf(angle)*t1 + cosf(angle)*t2));
    }
    mlg->triangulate(vertices);

    return mlg;
}
} // namespace
