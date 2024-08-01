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

in vec2 frag_texcoord;
uniform bool useTexcoords = false;

uniform sampler2D colorTex_;
uniform TextureParameters texParams_;

#ifndef NO_DEPTH_TEX
uniform sampler2D depthTex_;
#endif

void main() {
    vec2 fragCoord;
    if (useTexcoords) {
        fragCoord = frag_texcoord;
    } else {
        fragCoord = gl_FragCoord.xy * screenDimRCP_;
    }
    vec4 fragColor = textureLookup2Dnormalized(colorTex_, texParams_, fragCoord);
#ifdef LUMINANCE_TEXTURE
    FragData0 = vec4(fragColor.rgb, fragColor.r > 0.0 ? 1.0 : 0.0);
#else
    FragData0 = fragColor;
#endif

#ifndef NO_DEPTH_TEX
    gl_FragDepth = textureLookup2Dnormalized(depthTex_, texParams_, fragCoord).x;
 #endif
}
