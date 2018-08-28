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

uniform sampler2D entryPoints_;      // ray entry points
uniform sampler2D entryPointsDepth_; // depth values of entry points
uniform TextureParameters entryParameters_;
uniform sampler2D exitPoints_;          // ray exit points
uniform TextureParameters exitParameters_;
uniform sampler2D jitterTexture_;    // 8 bit random values
uniform TextureParameters jitterParameters_;

uniform float stepLength_;                  // raycasting step length

float permute(float i) {
    return mod((62.0*i*i + i), 961.0); // permutation polynomial; 62=2*31; 961=31*31
}

/**
 * Jitter entry points in ray direction
 */
void main() {
    vec2 p = gl_FragCoord.xy;
    p *= screenDimRCP_;
    vec4 frontPos = textureLookup2Dnormalized(entryPoints_, entryParameters_, p);
    vec4 backPos = textureLookup2Dnormalized(exitPoints_, exitParameters_, p);

    vec4 fragColor;
    float fragDepth;

    if ( (backPos.a == 0.0) || (frontPos.a == 0.0) ) {
        fragColor = frontPos;
        fragDepth = 1.0;
    } else {
        float rayLength = length(backPos.rgb - frontPos.rgb);
        // do not jitter very short rays
        if (rayLength <= stepLength_) {
            fragColor = frontPos;
            fragDepth = textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p ).x;
        }
        else {
            vec3 dir = (backPos.rgb - frontPos.rgb)/rayLength;
            float x = gl_FragCoord.x;
            float y = gl_FragCoord.y;
            float jitterValue = float(permute(x + permute(y))) / 961.0;
            vec3 frontPosNew = frontPos.rgb + (jitterValue*stepLength_)*dir;
            // save jitter value in alpha channel (for inverting jittering later if necessary)
            fragColor = vec4(frontPosNew, 1.0 - jitterValue*stepLength_);
            fragDepth = textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p ).x;
        }
    }

    FragData0 = fragColor;
    gl_FragDepth = fragDepth;
}
