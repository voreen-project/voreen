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

uniform sampler2D colorTex0_;
uniform sampler2D depthTex0_;
uniform TextureParameters textureParameters0_;
uniform sampler2D colorTex1_;
uniform sampler2D depthTex1_;

void main() {
    vec2 fragCoord = gl_FragCoord.xy;

    // fetch input textures
    vec4 color0 = textureLookup2Dscreen(colorTex0_, textureParameters0_, fragCoord);
    float depth0 = textureLookup2Dscreen(depthTex0_, textureParameters0_, fragCoord).x;
    vec4 color1 = texture2D(colorTex1_,gl_TexCoord[0].xy);
    float depth1 = texture2D(depthTex1_, gl_TexCoord[0].xy).x;

    vec4 fragColor = vec4(0.0);
    float fragDepth = 1.0;


    fragColor.rgb = color1.a * color1.rgb + (1.0 - color1.a) * color0.rgb;
    fragColor.a = color1.a + color0.a * (1.0 - color1.a);
    if (color1.a > 0)
        fragDepth = depth1;
    else
        fragDepth = depth0;

    FragData0 = fragColor;
    gl_FragDepth = fragDepth;
}
