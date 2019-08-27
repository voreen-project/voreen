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

uniform ivec3 dimensions_;
uniform ivec2 projectionRange_;
uniform mat3 realToProjectedMat_;
uniform mat3 projectedToRealMat_;

uniform sampler3D volumeTex_;
uniform sampler2D labelTex_;

in vec4 frag_texcoord;

out vec4 color;

void main() {
    ivec3 projectionDim = ivec3(realToProjectedMat_ * vec3(dimensions_));
    vec2 fragCoord = frag_texcoord.xy;

    float projectionValue = -1e10;
    for(int d=projectionRange_.x; d<=projectionRange_.y; ++d) {
        vec3 coord = projectedToRealMat_ * vec3(fragCoord, float(d)/float(projectionDim.z));
        float texValue = texture(volumeTex_, coord).x;
        projectionValue = max(texValue, projectionValue);
    }
    int labelValue = int(texture(labelTex_, fragCoord).x*255);
    vec3 col = vec3(projectionValue);
    vec3 mixColor;
    float mixFactor = 0.5;
    if(labelValue == UNLABELED) {
        mixColor = col;
        mixFactor = 0.0;
    } else if(labelValue == FOREGROUND) {
        mixColor = vec3(0.0, 1.0, 0.0);
    } else if (labelValue == BACKGROUND) {
        mixColor = vec3(1.0, 0.0, 0.0);
    } else if (labelValue == SUGGESTED_FOREGROUND) {
        mixColor = vec3(0.0, 0.0, 1.0);
        mixFactor = 0.3;
    } else {
        mixColor = vec3(1.0, 0.0, 1.0);
        mixFactor = 1.0;
    }
    col = mix(col, mixColor, mixFactor);
    color = vec4(col, 1.0);
}
