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

uniform sampler2D colorTex0_;
uniform sampler2D depthTex0_;
uniform TextureParameters textureParameters0_;
uniform sampler2D colorTex1_;
uniform sampler2D depthTex1_;
uniform TextureParameters textureParameters1_;

uniform int compositingMode_;
uniform float weightingFactor_;
uniform float weightFirst_;
uniform float weightSecond_;
uniform bool addDepth_;

void main() {
    vec2 fragCoord = gl_FragCoord.xy;

    // fetch input textures
    vec4 color0 = textureLookup2Dscreen(colorTex0_, textureParameters0_, fragCoord);
    float depth0 = textureLookup2Dscreen(depthTex0_, textureParameters0_, fragCoord).x;
    vec4 color1 = textureLookup2Dscreen(colorTex1_, textureParameters1_, fragCoord);
    float depth1 = textureLookup2Dscreen(depthTex1_, textureParameters1_, fragCoord).x;

    vec4 fragColor = vec4(0.0);
    float fragDepth = 1.0;

    #if defined(MODE_DEPTH_TEST)
        if (depth0 < depth1) {
            fragColor = color0;
            fragDepth = depth0;
        }
        else {
            fragColor = color1;
            fragDepth = depth1;
        }
    #elif defined(MODE_ALPHA_COMPOSITING)
        if (depth1 <= depth0) {
            fragColor.rgb = color1.rgb + color0.rgb * color0.a * (1.0 - color1.a);
            fragColor.a = color1.a + color0.a * (1.0 - color1.a);
        }
        else {
            fragColor.rgb = color0.rgb + color1.rgb * color1.a * (1.0 - color0.a);
            fragColor.a = color0.a + color1.a * (1.0 - color0.a);
        }
        fragDepth = min(depth0, depth1);
    #elif defined(MODE_ALPHA_BLENDING)
        fragColor.rgb = color0.a * color0.rgb + (1.0 - color0.a) * color1.rgb;
        fragColor.a = color0.a + color1.a * (1.0 - color0.a);
        if (color0.a > 0)
            fragDepth = depth0;
        else
            fragDepth = depth1;
    #elif defined(MODE_ALPHA_BLENDING_B_OVER_A)
        fragColor.rgb = color1.a * color1.rgb + (1.0 - color1.a) * color0.rgb;
        fragColor.a = color1.a + color0.a * (1.0 - color1.a);
        if (color1.a > 0)
            fragDepth = depth1;
        else
            fragDepth = depth0;
    #elif defined(MODE_WEIGHTED_AVERAGE)
        if (weightingFactor_ == 0.0)
            fragDepth = depth1;
        else if (weightingFactor_ == 1.0)
            fragDepth = depth0;
        else
            fragDepth = min(depth0, depth1);
        fragColor = weightingFactor_*color0 + (1.0-weightingFactor_)*color1;
    #elif defined(MODE_TAKE_FIRST)
        fragColor = color0;
        fragDepth = depth0;
    #elif defined(MODE_TAKE_SECOND)
        fragColor = color1;
        fragDepth = depth1;
    #elif defined(MODE_TAKE_SECOND_IF_READY)
        fragColor = color0;
        fragDepth = depth0;
    #elif defined(MODE_FIRST_HAS_PRIORITY)
        if (color0.a > 0.0) {
            fragColor = color0;
            fragDepth = depth0;
        }
        else {
            fragColor = color1;
            fragDepth = depth1;
        }
    #elif defined(MODE_SECOND_HAS_PRIORITY)
        // second has priority
        if (color1.a > 0.0) {
            fragColor = color1;
            fragDepth = depth1;
        }
        else {
            fragColor = color0;
            fragDepth = depth0;
        }
    #elif defined(MODE_MAXIMUM_ALPHA)
        fragDepth = min(depth0, depth1);
        fragColor.rgb = color0.rgb*color0.a + color1.rgb*color1.a;
        fragColor.a = max(color0.a, color1.a);
    #elif defined(MODE_DIFFERENCE)
        fragDepth = abs(depth0 - depth1);
        vec3 diff = vec3(1)-abs(color0.rgb-color1.rgb);
        fragColor.rgba = vec4(diff, max(color0.a,color1.a));
    #elif defined(MODE_DIFFERENCE_ALPHA_MUL)
        fragDepth = min(depth0, depth1);
        vec3 diff = abs(color0.rgb*color0.a-color1.rgb*color1.a);
        fragColor.rgba = vec4(diff, 1.0f);
    #elif defined(MODE_ADD)
        fragColor = weightFirst_*color0 + weightSecond_*color1;
        if (!addDepth_)
            fragDepth = depth0;
        else
            fragDepth = weightFirst_*depth0 + weightSecond_*depth1;
    #endif

    FragData0 = fragColor;
    gl_FragDepth = fragDepth;
}
