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

//#version 430 // is set in application as a header to the shader

layout(triangles) in;
layout(points, max_vertices = 3) out;

//position in model coordinates
in VertexData {
    vec4 position;
    vec4 color;
} vin[];

//clip plane (in model coordinates)
uniform vec4 plane;

// write out the edge as a single point to guarantee subsequent vertices for each edge
layout(stream = 0) out edgeData {
    vec4 position1;
    vec4 color1;
    vec4 position2;
    vec4 color2;
} eOut;

void main() {
    // distances to the plane
    precise vec3 dist = vec3(dot(plane.xyz, gl_in[0].gl_Position.xyz) - plane.w,
                     dot(plane.xyz, gl_in[1].gl_Position.xyz) - plane.w,
                     dot(plane.xyz, gl_in[2].gl_Position.xyz) - plane.w);

    if (dist == vec3(0.0)) {
        // the whole triangle lies in the plane => output all three edges
        for (int i = 0; i < 3; i++) {
            eOut.position1 = vin[i].position;
            eOut.color1 = vin[i].color;
            eOut.position2 = vin[(i+1)%3].position;
            eOut.color2 = vin[(i+1)%3].color;
            EmitStreamVertex(0);
            EndStreamPrimitive(0);
        }
        return;
    }

    if (all(greaterThan(dist,vec3(0.0)))) {
        // Case 1 (none clipped)
        return;
    }

    if (!any(greaterThanEqual(dist, vec3(0.0)))) {
        // Case 2 (all clipped)
        return;
    }

    // there are two special cases that have to be handled when only either v0 or v1 touches the plane
    for (int i = 0; i < 2; i++) {
        if ((dist[i] == 0.0) && (dist[i+1] != 0.0) && (dist[(i+2)%3] != 0.0)) {
            // write out the single touching point
            eOut.position1 = vin[i].position;
            eOut.color1 = vin[i].color;
            eOut.position2 = vin[i].position;
            eOut.color2 = vin[i].color;
            EmitStreamVertex(0);
            EndStreamPrimitive(0);
            return;
        }
    }

    // There are either 1 or 2 vertices above the clipping plane.
    ivec3 indices = ivec3(0,1,2);
    bvec3 above = greaterThanEqual(dist, vec3(0.0));
    bool nextIsAbove;

    // Find the CCW - most vertex above the plane.
    if (above[1] && !above[0]) {
        // Cycle once CCW .
        nextIsAbove = above[2];
        indices = ivec3(1,2,0);
        dist = dist.yzx ;
    } else if (above[2] && !above[1]) {
        // Cycle once CW.
        nextIsAbove = above[0];
        indices = ivec3(2,0,1);
        dist = dist.zxy ;
    } else
        nextIsAbove = above[1];

    //values for clipped edge vertices
    precise vec4 ePos1;
    precise vec4 ePos2;
    precise vec4 eColor1;
    precise vec4 eColor2;

    // We always need to clip v2 - v0.

    // prevent problems with division by zero when all vertices are above the plane
    float interpolationFactor2 = ((dist[0] == 0.0) || ((dist[0] - dist[2]) == 0.0)) ? 0.0 : dist[0] / (dist[0] - dist[2]);
    ePos2 = mix(vin[indices[0]].position, vin[indices[2]].position, interpolationFactor2);

    eColor2 = mix(vin[indices[0]].color, vin[indices[2]].color, interpolationFactor2);

    if (nextIsAbove) {
        // geometric case 3

        // prevent problems with division by zero when all vertices are above the plane
        float interpolationFactor1 = ((dist[1] == 0.0) || ((dist[1] - dist[2]) == 0.0)) ? 0.0 : dist[1] / (dist[1] - dist[2]);
        ePos1 = mix(vin[indices[1]].position, vin[indices[2]].position, interpolationFactor1);

        eColor1 = mix(vin[indices[1]].color, vin[indices[2]].color, interpolationFactor1);

    } else {
        // geometric case 4

        // prevent problems with division by zero when all vertices are above the plane
        float interpolationFactor1 = ((dist[0] == 0.0) || ((dist[0] - dist[1]) == 0.0)) ? 0.0 : dist[0] / (dist[0] - dist[1]);
        ePos1 = mix(vin[indices[0]].position, vin[indices[1]].position, interpolationFactor1);

        eColor1 = mix(vin[indices[0]].color, vin[indices[1]].color, interpolationFactor1);
    }

    // write out the clipping edge
    // correct orientation of edge is p2 -> p1
    eOut.position1 = ePos2;
    eOut.color1 = eColor2;
    eOut.position2 = ePos1;
    eOut.color2 = eColor1;
    EmitStreamVertex(0);
    EndStreamPrimitive(0);
}
