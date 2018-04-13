/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#version 330
#extension GL_ARB_arrays_of_arrays : require
#extension GL_ARB_shading_language_420pack : require

layout(points) in;
layout(triangle_strip, max_vertices = 15) out;


in uvec3[] voxelPos;
out vec3 pos;
out vec3 normal;

uniform ivec3 volumeDimensions_;
uniform float isoValue_;
uniform sampler3D colorTex_;

#define uint16_t int
#define int8_t int
#define uint8_t int
#define size_t int

// UBOs
layout(std140) uniform EDGE_INDICES_BLOCK{
    //uint EDGE_INDICES[256];
    uvec4 EDGE_INDICES[64];
};

int getEdgeIndex(int i) {
    return int(EDGE_INDICES[i/4][i%4]);
}
layout(std140) uniform TRIANGLE_VERTEX_INDICES_BLOCK{
    //int TRIANGLE_VERTEX_INDICES[256][16];
    ivec4 TRIANGLE_VERTEX_INDICES[256][4];
};
int getTriangleVertexIndex(int r, int c) {
    return TRIANGLE_VERTEX_INDICES[r][c/4][c%4];
}


const uvec3 VERTEX_OFFSETS[] = {
    uvec3(0, 0, 0),
    uvec3(1, 0, 0),
    uvec3(1, 1, 0),
    uvec3(0, 1, 0),
    uvec3(0, 0, 1),
    uvec3(1, 0, 1),
    uvec3(1, 1, 1),
    uvec3(0, 1, 1),
};

vec3 texPos(in vec3 vp) {
    return (vp+vec3(0.5,0.5,0.5))/volumeDimensions_;
}
float textureLookup(in vec3 voxelcoord) {
    return texture(colorTex_, texPos(voxelcoord)).x;
}
vec3 normalAt(in vec3 p) {
    float val_xp = textureLookup(p + vec3( 1, 0, 0));
    float val_xm = textureLookup(p + vec3(-1, 0, 0));
    float val_yp = textureLookup(p + vec3( 0, 1, 0));
    float val_ym = textureLookup(p + vec3( 0,-1, 0));
    float val_zp = textureLookup(p + vec3( 0, 0, 1));
    float val_zm = textureLookup(p + vec3( 0, 0,-1));


    // Central differences seem to provide better results in edge cases
    float dv_dx = val_xm - val_xp;
    float dv_dy = val_ym - val_yp;
    float dv_dz = val_zm - val_zp;
    vec3 n = vec3(dv_dx, dv_dy, dv_dz);
    if(n == vec3(0,0,0)) {
        // Central differences failed. Try forward differences
        float central_val = textureLookup(p); //Should be about equal to isovalue, but we better get the exact value
        n.x = central_val - val_xp;
        n.y = central_val - val_yp;
        n.z = central_val - val_zp;
        if(n == vec3(0,0,0)) {
            // Forward differences failed, too. Try backward differences
            n.x = val_xm - central_val;
            n.y = val_ym - central_val;
            n.z = val_zm - central_val;
            if(n == vec3(0,0,0)) {
                // Well, none of the above worked. Just set the normal to (1,0,0)
                n = vec3(1,0,0);
            }
        }
    }
    return n;
}
void prepareVertex(in uvec3 p1, in uvec3 p2, out vec3 pos, out vec3 normal) {
    float v_p1 = textureLookup(vec3(p1));
    float v_p2 = textureLookup(vec3(p2));
    float a = (isoValue_ - v_p2) / (v_p1 - v_p2);
    //float a = 0.5f;
    pos = a*vec3(p1) + (1-a)*vec3(p2);
    normal = normalAt(pos);
    if(normal != vec3(0,0,0)) {
        normal = normalize(normal);
    } else {
        //LWARNINGC("voreen.surfaceextraction.surfaceectractor", "normal = 0, setting to (1,0,0)");
        normal = vec3(1,0,0);
    }
}
void main() {
    uvec3 p = voxelPos[0];


    int index = 0;
    if(textureLookup(p + VERTEX_OFFSETS[0]) < isoValue_) index |= (1 << 0);
    if(textureLookup(p + VERTEX_OFFSETS[1]) < isoValue_) index |= (1 << 1);
    if(textureLookup(p + VERTEX_OFFSETS[2]) < isoValue_) index |= (1 << 2);
    if(textureLookup(p + VERTEX_OFFSETS[3]) < isoValue_) index |= (1 << 3);
    if(textureLookup(p + VERTEX_OFFSETS[4]) < isoValue_) index |= (1 << 4);
    if(textureLookup(p + VERTEX_OFFSETS[5]) < isoValue_) index |= (1 << 5);
    if(textureLookup(p + VERTEX_OFFSETS[6]) < isoValue_) index |= (1 << 6);
    if(textureLookup(p + VERTEX_OFFSETS[7]) < isoValue_) index |= (1 << 7);

    if(index == 0 || index == 0xff) {
        return;
    }

    vec3 positions[12];
    vec3 normals[12];
    const int edgeIndex = getEdgeIndex(index);

    if((edgeIndex & (1 << 0x0)) != 0) prepareVertex(p + VERTEX_OFFSETS[0], p + VERTEX_OFFSETS[1], positions[0x0], normals[0x0]);
    if((edgeIndex & (1 << 0x1)) != 0) prepareVertex(p + VERTEX_OFFSETS[1], p + VERTEX_OFFSETS[2], positions[0x1], normals[0x1]);
    if((edgeIndex & (1 << 0x2)) != 0) prepareVertex(p + VERTEX_OFFSETS[2], p + VERTEX_OFFSETS[3], positions[0x2], normals[0x2]);
    if((edgeIndex & (1 << 0x3)) != 0) prepareVertex(p + VERTEX_OFFSETS[3], p + VERTEX_OFFSETS[0], positions[0x3], normals[0x3]);
    if((edgeIndex & (1 << 0x4)) != 0) prepareVertex(p + VERTEX_OFFSETS[4], p + VERTEX_OFFSETS[5], positions[0x4], normals[0x4]);
    if((edgeIndex & (1 << 0x5)) != 0) prepareVertex(p + VERTEX_OFFSETS[5], p + VERTEX_OFFSETS[6], positions[0x5], normals[0x5]);
    if((edgeIndex & (1 << 0x6)) != 0) prepareVertex(p + VERTEX_OFFSETS[6], p + VERTEX_OFFSETS[7], positions[0x6], normals[0x6]);
    if((edgeIndex & (1 << 0x7)) != 0) prepareVertex(p + VERTEX_OFFSETS[7], p + VERTEX_OFFSETS[4], positions[0x7], normals[0x7]);
    if((edgeIndex & (1 << 0x8)) != 0) prepareVertex(p + VERTEX_OFFSETS[0], p + VERTEX_OFFSETS[4], positions[0x8], normals[0x8]);
    if((edgeIndex & (1 << 0x9)) != 0) prepareVertex(p + VERTEX_OFFSETS[1], p + VERTEX_OFFSETS[5], positions[0x9], normals[0x9]);
    if((edgeIndex & (1 << 0xA)) != 0) prepareVertex(p + VERTEX_OFFSETS[2], p + VERTEX_OFFSETS[6], positions[0xA], normals[0xA]);
    if((edgeIndex & (1 << 0xB)) != 0) prepareVertex(p + VERTEX_OFFSETS[3], p + VERTEX_OFFSETS[7], positions[0xB], normals[0xB]);

#define vertices(i) getTriangleVertexIndex(index, i)

    for(int vertexNum = 0; vertices(vertexNum) != -1; vertexNum += 3) {
        int8_t v0 = vertices(vertexNum    );
        pos = positions[v0];
        normal = normals[v0];
        EmitVertex();

        int8_t v1 = vertices(vertexNum + 1);
        pos = positions[v1];
        normal = normals[v1];
        EmitVertex();

        int8_t v2 = vertices(vertexNum + 2);
        pos = positions[v2];
        normal = normals[v2];
        EmitVertex();

        EndPrimitive();
    }
}
