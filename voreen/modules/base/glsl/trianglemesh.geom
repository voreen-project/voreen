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

#version 330

layout(triangles) in;
layout(triangle_strip, max_vertices=4) out;

in VertexData {
  vec4 color;
#if defined(USE_NORMAL)
  vec3 normal;
#endif
#if defined(USE_TEXCOORD)
  vec2  texCoord;
  ivec2 texIndex;
#endif
} vertexData[];

struct FragStructHelper {
  vec4 posEye;
  vec4 color;
#if defined(USE_NORMAL)
  vec3 normal;
#endif
#if defined(USE_TEXCOORD)
  vec2  texCoord;
  ivec2 texIndex;
#endif
};

out FragStruct {
  vec4 posEye;
  vec4 color;
#if defined(USE_NORMAL)
  vec3 normal;
#endif
#if defined(USE_TEXCOORD)
  vec2  texCoord;
  flat ivec2 texIndex;
#endif
} fragData;

uniform mat4 viewMatrix_;
uniform mat4 projectionMatrix_;
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
        fragData.posEye = viewMatrix_ * gl_in[0].gl_Position;
        fragData.color = vertexData[0].color;
#if defined(USE_NORMAL)
        fragData.normal = vertexData[0].normal;
#endif
#if defined(USE_TEXCOORD)
        fragData.texCoord = vertexData[0].texCoord;
        fragData.texIndex = vertexData[0].texIndex;
#endif
        gl_Position = projectionMatrix_ * fragData.posEye; EmitVertex();
        fragData.posEye = viewMatrix_ * gl_in[1].gl_Position;
        fragData.color = vertexData[1].color;
#if defined(USE_NORMAL)
        fragData.normal = vertexData[1].normal;
#endif
#if defined(USE_TEXCOORD)
        fragData.texCoord = vertexData[1].texCoord;
        fragData.texIndex = vertexData[1].texIndex;
#endif
        gl_Position = projectionMatrix_ * fragData.posEye; EmitVertex();
        fragData.posEye = viewMatrix_ * gl_in[2].gl_Position;
        fragData.color = vertexData[2].color;
#if defined(USE_NORMAL)
        fragData.normal = vertexData[2].normal;
#endif
#if defined(USE_TEXCOORD)
        fragData.texCoord = vertexData[2].texCoord;
        fragData.texIndex = vertexData[2].texIndex;
#endif
        gl_Position = projectionMatrix_ * fragData.posEye; EmitVertex();
        EndPrimitive();
        return;
    }

    distances -= plane_.w;
    float lastDistance = distances.x;

    FragStructHelper[4] outVertices;
    int outIndex = 0;
    const float epsilon = 0.0001;

    // Process face edges...
    for (int i = 0; i < 3; ++i) {
        float dist = distances[(i + 1) % 3];

        // Keep both vertices?
        if (lastDistance <= 0 && dist <= 0) {
            // If processing the first edge, insert first vertex...
            if (i == 0) {
                outVertices[outIndex].posEye = viewMatrix_ * gl_in[i].gl_Position;
                outVertices[outIndex].color = vertexData[i].color;
#if defined(USE_NORMAL)
                outVertices[outIndex].normal = vertexData[i].normal;
#endif
#if defined(USE_TEXCOORD)
                outVertices[outIndex].texCoord = vertexData[i].texCoord;
                outVertices[outIndex].texIndex = vertexData[i].texIndex;
#endif
                outIndex++;
            }

            // If NOT processing the last edge, insert second vertex...
            if (i < 2) {
                outVertices[outIndex].posEye = viewMatrix_ * gl_in[i + 1].gl_Position;
                outVertices[outIndex].color = vertexData[i + 1].color;
#if defined(USE_NORMAL)
                outVertices[outIndex].normal = vertexData[i + 1].normal;
#endif
#if defined(USE_TEXCOORD)
                outVertices[outIndex].texCoord = vertexData[i + 1].texCoord;
                outVertices[outIndex].texIndex = vertexData[i + 1].texIndex;
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
                outVertices[outIndex].posEye = viewMatrix_ * vec4(intersectionPos, 1.0);
                outVertices[outIndex].color  = mix(vertexData[i].color, vertexData[(i + 1) % 3].color, fac);
#if defined(USE_NORMAL)
                outVertices[outIndex].normal = mix(vertexData[i].normal, vertexData[(i + 1) % 3].normal, fac);
#endif
#if defined(USE_TEXCOORD)
                // mix only tex coords, as texIndex contains texture layer and function integers which are used with bit arithmetic later
                outVertices[outIndex].texCoord = mix(vertexData[i].texCoord.xy, vertexData[(i + 1) % 3].texCoord.xy, fac);
                outVertices[outIndex].texIndex = vertexData[i].texIndex;
#endif
                outIndex++;
            }

            // If NOT processing the last edge, insert second vertex...
            if (i < 2) {
                outVertices[outIndex].posEye = viewMatrix_ * gl_in[i + 1].gl_Position;
                outVertices[outIndex].color = vertexData[i + 1].color;
#if defined(USE_NORMAL)
                outVertices[outIndex].normal = vertexData[i + 1].normal;
#endif
#if defined(USE_TEXCOORD)
                outVertices[outIndex].texCoord = vertexData[i + 1].texCoord;
                outVertices[outIndex].texIndex = vertexData[i + 1].texIndex;
#endif
                outIndex++;
            }
        }
        // Keep first vertex, but discard second vertex?
        else if (lastDistance <= 0.0 && dist > 0.0) {
            // If processing the first edge, insert first vertex...
            if (i == 0) {
                outVertices[outIndex].posEye = viewMatrix_ * gl_in[i].gl_Position;
                outVertices[outIndex].color = vertexData[i].color;
#if defined(USE_NORMAL)
                outVertices[outIndex].normal = vertexData[i].normal;
#endif
#if defined(USE_TEXCOORD)
                outVertices[outIndex].texCoord = vertexData[i].texCoord;
                outVertices[outIndex].texIndex = vertexData[i].texIndex;
#endif
                outIndex++;
            }

            // If NOT clipplane intersection vertex and first vertex are equal, insert clipplane intersection vertex...
            float fac = lastDistance / (lastDistance - dist);
            vec3 intersectionPos = mix(gl_in[i].gl_Position.xyz, gl_in[(i + 1) % 3].gl_Position.xyz, fac);
            if (distance(gl_in[i].gl_Position.xyz, intersectionPos) > epsilon) {
                outVertices[outIndex].posEye = viewMatrix_ * vec4(intersectionPos, 1.0);
                outVertices[outIndex].color  = mix(vertexData[i].color, vertexData[(i + 1) % 3].color, fac);
#if defined(USE_NORMAL)
                outVertices[outIndex].normal = mix(vertexData[i].normal, vertexData[(i + 1) % 3].normal, fac);
#endif
#if defined(USE_TEXCOORD)
                // mix only tex coords, as texIndex contains texture layer and function integers which are used with bit arithmetic later
                outVertices[outIndex].texCoord = mix(vertexData[i].texCoord.xy, vertexData[(i + 1) % 3].texCoord.xy, fac);
                outVertices[outIndex].texIndex = vertexData[i].texIndex;
#endif
                outIndex++;
            }
        }

        lastDistance = dist;
    }

    // Create triangles from output vertices:
    if(outIndex == 3) {
        for(int i = 0; i < 3; i++) {
            fragData.posEye = outVertices[i].posEye;
            fragData.color  = outVertices[i].color;
            #if defined(USE_NORMAL)
            fragData.normal = outVertices[i].normal;
            #endif
            #if defined(USE_TEXCOORD)
            fragData.texCoord = outVertices[i].texCoord;
            fragData.texIndex = outVertices[i].texIndex;
            #endif
            gl_Position = projectionMatrix_ * fragData.posEye; EmitVertex();
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
            fragData.posEye = outVertices[index].posEye;
            fragData.color  = outVertices[index].color;
            #if defined(USE_NORMAL)
            fragData.normal = outVertices[index].normal;
            #endif
            #if defined(USE_TEXCOORD)
            fragData.texCoord = outVertices[index].texCoord;
            fragData.texIndex = outVertices[index].texIndex;
            #endif
            gl_Position = projectionMatrix_ * fragData.posEye; EmitVertex();
        }
        EndPrimitive();
        return;
    }
}

