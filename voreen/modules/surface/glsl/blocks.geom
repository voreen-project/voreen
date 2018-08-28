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

#version 330

layout(points) in;
// Reasoning for max_vertices = 24.
// A single cube would have 6 sides with 4 vertices each => 24
// Although in general we extract only 3 sides for each voxel,
// the voxel at (0, 0, 0) might have 6.
// (See main, max. 6 calls to tryAddFace.)
layout(triangle_strip, max_vertices = 24) out;


in ivec3[] voxelPos;
out vec3 pos;
out vec3 normal;

uniform ivec3 volumeDimensions_;
uniform float isoValue_;
uniform sampler3D colorTex_;

vec3 texPos(in vec3 vp) {
    return (vp+vec3(0.5,0.5,0.5))/volumeDimensions_;
}

float textureLookup(in vec3 voxelcoord) {
    return texture(colorTex_, texPos(voxelcoord)).x;
}

void tryAddFace(ivec3 p, ivec3 dp, uint centerSmaller) {
    uint thisSmaller = ((textureLookup(p + dp) < isoValue_)) ? 1u : 0u;
    ivec3 p2 = ivec3(p + dp);
    if((centerSmaller ^ thisSmaller) == 1u || centerSmaller == 0u && any(equal(p2, volumeDimensions_))) {
        vec3 lineDelta = 0.5 * vec3(dp); // Central position between voxels
        vec3 faceDelta1 = lineDelta.yzx; // \|/
        vec3 faceDelta2 = lineDelta.zxy; // Are orthogonal to each other and will span the face
        vec3 n = dp * (1 - 2*float(centerSmaller)); // Invert normal if centerSmaller

        pos = vec3(p) + lineDelta + faceDelta1 + faceDelta2;
        normal = n;
        EmitVertex();

        pos = vec3(p) + lineDelta - faceDelta1 + faceDelta2;
        normal = n;
        EmitVertex();

        pos = vec3(p) + lineDelta + faceDelta1 - faceDelta2;
        normal = n;
        EmitVertex();

        pos = vec3(p) + lineDelta - faceDelta1 - faceDelta2;
        normal = n;
        EmitVertex();

        EndPrimitive();
    }
}
// We assume a low value outside of the volume => smaller than iso value
#define OUTSIDE_VOLUME_SMALLER 1u

void main() {
    ivec3 p = voxelPos[0];

    int index = 0;
    uint centerSmaller = textureLookup(p) < isoValue_ ? 1u : 0u;
    // Possibly add all 3 forward facing faces.
    // Faces facing backwards (to origin) will (in general, see below) be added by other voxels.
    tryAddFace(p, ivec3(1,0,0), centerSmaller);
    tryAddFace(p, ivec3(0,1,0), centerSmaller);
    tryAddFace(p, ivec3(0,0,1), centerSmaller);

    // Special handling for faces at the 3 outside faces of the volume nearest to origin
    if(p.x == 0) {
        tryAddFace(p-ivec3(1,0,0), ivec3(1,0,0), OUTSIDE_VOLUME_SMALLER);
    }
    if(p.y == 0) {
        tryAddFace(p-ivec3(0,1,0), ivec3(0,1,0), OUTSIDE_VOLUME_SMALLER);
    }
    if(p.z == 0) {
        tryAddFace(p-ivec3(0,0,1), ivec3(0,0,1), OUTSIDE_VOLUME_SMALLER);
    }
}
