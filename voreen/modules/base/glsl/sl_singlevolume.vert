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

// Adapted from "Real-Time Volume Graphics"

#include "modules/mod_sampler3d.frag"

layout(location=0) in vec2 vertexPosition;

// uniforms needed for clipping code
uniform int frontIdx_;
uniform vec3 vecView_;
uniform float dPlaneStart_;
uniform float dPlaneIncr_;
uniform int[64] nSequence_;
uniform vec3[8] vecVertices_;
uniform int[24] v1_;
uniform int[24] v2_;

uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;    // texture lookup parameters for volume_

out vec4 textureCoordinate[2];

void main() {
    vec3 position = vec3(0.0);
    float dPlane = dPlaneStart_ + vertexPosition.y * dPlaneIncr_;

    for (int e=0; e<4; e++) {
        int vidx1 = nSequence_[int(frontIdx_*8 + v1_[int(vertexPosition.x*4 + e)])];
        int vidx2 = nSequence_[int(frontIdx_*8 + v2_[int(vertexPosition.x*4 + e)])];
        vec3 vecV1 = vecVertices_[vidx1];
        vec3 vecV2 = vecVertices_[vidx2];

        vec3 vecDir = vecV2-vecV1;
        float denom = dot(vecDir,vecView_);

        if (denom == 0.0)
            continue;

        float lambda = (dPlane - dot(vecV1, vecView_)) / denom;
        if (lambda >= 0.0 && lambda <= 1.0) {
            position = vecV1 + lambda * vecDir;
            break;
        }
    }

    vec4 worldPos = volumeStruct_.physicalToWorldMatrix_ * vec4(position, 1.0);
    // save texture coordinate
    textureCoordinate[0] = volumeStruct_.worldToTextureMatrix_ * worldPos;
    // save physical CS coordinate
    textureCoordinate[1] = vec4(position, 1.0);
    // write clip coordinates
    gl_Position = modelViewProjectionMatrixStack_ * vec4(worldPos.xyz, 1.0);
}
