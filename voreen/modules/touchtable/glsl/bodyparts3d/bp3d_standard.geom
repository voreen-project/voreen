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

layout(triangles) in;
layout(triangle_strip, max_vertices=4) out;

in VertexData {
  vec4 color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
  vec3 normal;
#endif
} vertexData[];

struct FragStruct {
  vec4 posEye;
  vec4 color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
  vec3 normal;
#endif
};

out FragStruct fragData;

uniform mat4 MVP;
uniform vec4 plane_;
uniform bool enableClipping_;

void main() {

    bvec3 vertPlaneStates;
    vec3 distances;

    if(!enableClipping_)
        vertPlaneStates = bvec3(false);
    else {
        distances = vec3(dot(gl_in[0].gl_Position.xyz, plane_.xyz),
                         dot(gl_in[1].gl_Position.xyz, plane_.xyz),
                         dot(gl_in[2].gl_Position.xyz, plane_.xyz));

        vertPlaneStates = bvec3(distances.x > plane_.w, distances.y > plane_.w, distances.z > plane_.w);
    }

    if(all(vertPlaneStates))
        return;

    if(all(not(vertPlaneStates))) {
        fragData.posEye = gl_in[0].gl_Position;
        fragData.color = vertexData[0].color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
        fragData.normal = vertexData[0].normal;
#endif
        gl_Position = MVP * fragData.posEye; EmitVertex();
        fragData.posEye = gl_in[1].gl_Position;
        fragData.color = vertexData[1].color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
        fragData.normal = vertexData[1].normal;
#endif
        gl_Position = MVP * fragData.posEye; EmitVertex();
        fragData.posEye = gl_in[2].gl_Position;
        fragData.color = vertexData[2].color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
        fragData.normal = vertexData[2].normal;
#endif
        gl_Position = MVP * fragData.posEye; EmitVertex();
        EndPrimitive();
        return;
    }

    distances -= plane_.w;
    float lastDistance = distances.x;

    FragStruct[4] outVertices;
    int outIndex = 0;
    float epsilon = 0.0001;

    // Process face edges...
    for (int i = 0; i < 3; ++i) {
        float dist = distances[(i + 1) % 3];

        // Keep both vertices?
        if (lastDistance <= 0 && dist <= 0) {
            // If processing the first edge, insert first vertex...
            if (i == 0) {
                outVertices[outIndex].posEye = gl_in[i].gl_Position;
                outVertices[outIndex].color = vertexData[i].color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
                outVertices[outIndex].normal = vertexData[i].normal;
#endif
                outIndex++;
            }

            // If NOT processing the last edge, insert second vertex...
            if (i < 2) {
                outVertices[outIndex].posEye = gl_in[i + 1].gl_Position;
                outVertices[outIndex].color = vertexData[i + 1].color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
                outVertices[outIndex].normal = vertexData[i + 1].normal;
#endif
                outIndex++;
            }
        }
        // Discard first vertex, but keep second vertex?
        else if (lastDistance > 0.0 && dist <= 0.0) {
            // If NOT clipplane intersection vertex and second vertex are equal, insert clipplane intersection vertex...
            float fac = lastDistance / (lastDistance - dist);
            vec3 intersectionPos = mix(gl_in[i].gl_Position.xyz, gl_in[(i + 1) % 3].gl_Position.xyz, fac);
            if (distance(gl_in[(i + 1) % 3].gl_Position.xyz, intersectionPos) > epsilon) {
                outVertices[outIndex].posEye = vec4(intersectionPos, 1.0);
                outVertices[outIndex].color  = mix(vertexData[i].color, vertexData[(i + 1) % 3].color, fac);
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
                outVertices[outIndex].normal = mix(vertexData[i].normal, vertexData[(i + 1) % 3].normal, fac);
#endif
                outIndex++;
            }

            // If NOT processing the last edge, insert second vertex...
            if (i < 2) {
                outVertices[outIndex].posEye = gl_in[i + 1].gl_Position;
                outVertices[outIndex].color = vertexData[i + 1].color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
                outVertices[outIndex].normal = vertexData[i + 1].normal;
#endif
                outIndex++;
            }
        }
        // Keep first vertex, but discard second vertex?
        else if (lastDistance <= 0.0 && dist > 0.0) {
            // If processing the first edge, insert first vertex...
            if (i == 0) {
                outVertices[outIndex].posEye = gl_in[i].gl_Position;
                outVertices[outIndex].color = vertexData[i].color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
                outVertices[outIndex].normal = vertexData[i].normal;
#endif
                outIndex++;
            }

            // If NOT clipplane intersection vertex and first vertex are equal, insert clipplane intersection vertex...
            float fac = lastDistance / (lastDistance - dist);
            vec3 intersectionPos = mix(gl_in[i].gl_Position.xyz, gl_in[(i + 1) % 3].gl_Position.xyz, fac);
            if (distance(gl_in[i].gl_Position.xyz, intersectionPos) > epsilon) {
                outVertices[outIndex].posEye = vec4(intersectionPos, 1.0);
                outVertices[outIndex].color  = mix(vertexData[i].color, vertexData[(i + 1) % 3].color, fac);
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
                outVertices[outIndex].normal = mix(vertexData[i].normal, vertexData[(i + 1) % 3].normal, fac);
#endif
                outIndex++;
            }
        }

        lastDistance = dist;
    }

    // Create triangles from output vertices:
    if(outIndex == 3) {
        for(int i = 0; i < 3; i++) {
            fragData = outVertices[i];
            gl_Position = MVP * fragData.posEye; EmitVertex();
        }
        EndPrimitive();
        return;
    } else if(outIndex == 4) {
        for(int i = 0; i < 4; i++) {
            int index = i;
            if(i == 2)
                index = 3;
            if(i == 3)
                index = 2;
            fragData = outVertices[index];
            gl_Position = MVP * fragData.posEye; EmitVertex();
        }
        EndPrimitive();
        return;
    }
}

