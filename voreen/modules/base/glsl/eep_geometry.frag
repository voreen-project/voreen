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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_depth.frag"

uniform sampler2D entryParams_;
uniform sampler2D entryParamsDepth_;
uniform TextureParameters entryInfo_;
uniform sampler2D exitParams_;
uniform sampler2D exitParamsDepth_;
uniform TextureParameters exitInfo_;
uniform sampler2D geometryTex_;
uniform sampler2D geometryTexDepth_;
uniform TextureParameters geomInfo_;

uniform bool entry_;
uniform float near_;
uniform float far_;
uniform mat4 viewMatrixInverse_;
uniform mat4 projectionMatrixInverse_;
uniform bool useFloatTarget_;
uniform mat4 worldToTexture_;

void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;

    vec4 entryCol = textureLookup2Dnormalized(entryParams_, entryInfo_, p);
    float entryDepth = textureLookup2Dnormalized(entryParamsDepth_, entryInfo_, p).x;

    vec4 exitCol = textureLookup2Dnormalized(exitParams_, exitInfo_, p);
    float exitDepth = textureLookup2Dnormalized(exitParamsDepth_, exitInfo_, p).x;

    float geomDepth = textureLookup2Dnormalized(geometryTexDepth_, geomInfo_, p).x;

    if (entry_) {
        if (geomDepth <= entryDepth) {
            //FragData0 = exitCol;
            FragData0 = vec4(0.0);
            gl_FragDepth = geomDepth;
        } else {
            FragData0 = entryCol;
            gl_FragDepth = entryDepth;
        }
    } else {
        if (geomDepth <= entryDepth) {
            FragData0 = vec4(0.0);
            //FragData0 = exitCol;
            gl_FragDepth = exitDepth;
        } else if (geomDepth <= exitDepth) {

            // use NDCs from -1 to 1
            vec4 wGeom = vec4(2.0*p - 1.0, 2.0*geomDepth - 1.0, 1.0);

            // reverse perspective division by w (which is equal to the camera-space z)
            float origZG = 2.0*far_*near_ / ((far_ + near_) - wGeom.z*(far_ - near_));
            wGeom *= origZG;

            // unproject and reverse camera-transform
            vec4 result = vec4((viewMatrixInverse_ * (projectionMatrixInverse_ * wGeom)).xyz, 1.0);

            // For single-volume-raycasting we want standard (Krüger-Westermann) entry-exit points:
            if(!useFloatTarget_)
                result.xyz = (worldToTexture_ * result).xyz;
            // Otherwise we just use the world coordinates. (World-entry-exit-points)

            FragData0 = vec4(result.xyz, 1.0);

            gl_FragDepth = geomDepth;

        } else {
            FragData0 = exitCol;
            gl_FragDepth = exitDepth;
        }
    }
}
