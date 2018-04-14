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

#include "modules/mod_sampler2d.frag"

uniform sampler2D firstBack_;      // ray entry points
//uniform sampler2D firstBackDepth_; // depth values of entry points
uniform TextureParameters firstBackParameters_;

uniform sampler2D firstFront_;       // ray exit points
uniform sampler2D firstFrontDepth_; // depth values of entry points
uniform TextureParameters firstFrontParameters_;

uniform mat4 viewMatrixInverse_;
uniform mat4 projectionMatrixInverse_;
uniform bool useTextureCoordinates_;
uniform mat4 worldToTexture_;

uniform vec3 llf_;
uniform vec3 urb_;

void main() {
    vec2 p = gl_FragCoord.xy;
    p *= screenDimRCP_;

    vec4 firstBackPos = textureLookup2Dnormalized(firstBack_, firstBackParameters_, p);
    //float firstBackDepth = textureLookup2Dnormalized(firstBackDepth_, firstBackParameters_, p).z;

    vec4 firstFrontPos = textureLookup2Dnormalized(firstFront_, firstFrontParameters_, p);
    float firstFrontDepth = textureLookup2Dnormalized(firstFrontDepth_, firstFrontParameters_, p).x;

    vec4 fragColor;
    float fragDepth;

    if (firstBackPos.a == 0.0) { // we have no exit pos
        fragColor = vec4(0.0);
        fragDepth = 1.0;
    } else {
        // use NDCs from -1 to 1
        vec4 wGeom = vec4(2.0*p - 1.0, -1.0, 1.0);

        // unproject and reverse camera-transform
        vec4 tmpResult = viewMatrixInverse_ * (projectionMatrixInverse_ * wGeom);
        vec4 result = tmpResult / tmpResult.w;

        // For single-volume-raycasting we want standard (Krüger-Westermann) entry-exit points:
        if(useTextureCoordinates_)
            result.xyz = (worldToTexture_ * result).xyz;

        if(any(greaterThan(result.xyz, urb_)) || any(lessThan(result.xyz, llf_))) {
            fragColor = firstFrontPos;
            fragDepth = firstFrontDepth;
        } else {
            fragColor = result;
            fragDepth = 0.0;
        }
    }

    FragData0 = fragColor;
    gl_FragDepth = fragDepth;
}
