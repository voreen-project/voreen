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
#include "modules/mod_sampler3d.frag"

#include "modules/mod_depth.frag"

uniform sampler2D firstHitPoints_;             // first hit points
uniform TextureParameters firstHitParameters_;

uniform sampler2D entryPoints_;                // ray entry points
uniform sampler2D entryPointsDepth_;           // ray entry points depth
uniform TextureParameters entryParameters_;
uniform sampler2D exitPoints_;                 // ray exit points
uniform sampler2D exitPointsDepth_;            // ray exit points depth
uniform TextureParameters exitParameters_;

uniform VolumeParameters segmentationParameters_;
uniform sampler3D segmentation_;  // segmented dataset

uniform float penetrationDepth_;


vec4 fillIDBuffer(in vec4 entry, in vec4 firstHitPoint, in vec4 exit) {

    vec4 result;
    float seg = textureLookup3D(segmentation_, segmentationParameters_, firstHitPoint.rgb).a;
    result.rgb = firstHitPoint.rgb;
    // if not hit any segment, penetrate volume until a segment is hit
    // or penetration-depth is reached
    if (seg == 0.0) {
        vec3 direction = normalize(exit.rgb-entry.rgb);
        float stepIncr = 0.005;
        float t_add = 0.0;
        vec3 pos;

        while ((seg == 0.0) && (t_add < penetrationDepth_)) {
            t_add += stepIncr;
            pos = firstHitPoint.rgb + t_add*direction;
            seg = textureLookup3DUnnormalized(segmentation_, segmentationParameters_, pos).a;
        }
        result.rgb = pos;
    }
    if (seg > 0.0)
        result.a = seg;
    else
        result.a = 1.0;

    return result;
}

void main() {

    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec4 firstHitPos = textureLookup2Dnormalized(firstHitPoints_, firstHitParameters_, p);
    vec4 entryPos = textureLookup2Dnormalized(entryPoints_, entryParameters_, p);
    vec4 exitPos = textureLookup2Dnormalized(exitPoints_,exitParameters_,  p);
    float entryPointDepth = textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p).z;
    float exitPointDepth = textureLookup2Dnormalized(exitPointsDepth_, exitParameters_, p).z;

    vec4 result;
    if (firstHitPos.a != 0.0) {
           result = fillIDBuffer(entryPos, firstHitPos, exitPos);
           float t = length(result.rgb-entryPos.rgb)/length(exitPos.rgb-entryPos.rgb);
           gl_FragDepth = calculateDepthValue(t, entryPointDepth, exitPointDepth);
    }
    else {
        result=vec4(0.0);
        gl_FragDepth = 1.0;
    }

    FragData0 = result;
}
