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

#include "modules/mod_sampler2d.frag"

uniform sampler2D entryTex_;
uniform TextureParameters entryParameters_;
uniform sampler2D entryTexDepth_;
uniform sampler2D exitTex_;
uniform TextureParameters exitParameters_;

uniform vec3 volumeCubeSizeRCP_;
uniform vec3 volumeOffset_;

/**
 * Calculates missing entry points (clipped by near plane) by back-projecting the
 * fragment to the proxy-geometry's object coordinates.
 **/
void main() {

    vec2 fragCoord = gl_FragCoord.xy;
    fragCoord *= screenDimRCP_;
    vec4 exitCol = textureLookup2Dnormalized(exitTex_, exitParameters_, fragCoord);
    vec4 entryCol = textureLookup2Dnormalized(entryTex_, entryParameters_, fragCoord);

    if (exitCol.a != 0.0) {
        // For some reason the near plane clipping can cause alpha values that lie between 0.0
        // and 1.0, so do a conservative test here.
        if (entryCol.a >= 0.5) {
            // entry point is okay
            FragData0 = entryCol;
            gl_FragDepth = textureLookup2Dnormalized(entryTexDepth_, entryParameters_, fragCoord).x;
        }
        else {
            // alpha value of entry point less than 1.0 because of near plane clipping: fill holes

            // viewport coordinates -> NDC
            vec4 pos = vec4((gl_FragCoord.xy * screenDimRCP_) * 2.0 - vec2(1.0), -1.0, 1.0);

            // NDC -> object coordinates
            pos = gl_ModelViewMatrixInverse * gl_ProjectionMatrixInverse * pos;
            pos /= pos.w;

            // object coordinates -> texture coordinates
            vec3 eep = ((pos.xyz - volumeOffset_)* volumeCubeSizeRCP_) + 0.5;
            // write new entry point
            FragData0 = vec4(eep, 1.0);
            gl_FragDepth = 0.0;
        }
    }
    else {
        // exit point is 0.0, ignore
        discard;
    }
}
