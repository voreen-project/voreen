/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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
#include "modules/mod_normdepth.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters texParams_;

uniform sampler1D chromadepthTex_;
uniform float colorDepthFactor_;
uniform int colorMode_;

/***
 * The main method.
 ***/
void main() {
    vec4 fragCoord = gl_FragCoord;

    // read and normalize depth value
    float depth = textureLookup2Dscreen(depthTex_, texParams_, fragCoord.xy).x;
    float depthNorm = normDepth(depth) * (1.0 / colorDepthFactor_);

    vec4 result = vec4(0.0);

    if (colorMode_ == 0) {
        // light-dark (replace)
        if (depth != 1.0)
            result = vec4(vec3(1.0 - depthNorm), 1.0);
    } else if (colorMode_ == 1) {
        // light-dark (modulate)
        result = vec4(vec3(1.0-depthNorm), 1.0)*textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy);
    } else if (colorMode_ == 2) {
        // chromadepth
      #if defined(GLSL_VERSION_130)
        result = texture(chromadepthTex_, depthNorm)*textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy);
      #else
        result = texture1D(chromadepthTex_, depthNorm)*textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy);
      #endif
    } else if (colorMode_ == 3) {
        // pseudochromadepth
        result = vec4(1.0-depthNorm, 0.0, depthNorm, 1.0)*textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy);
    }

    FragData0 = result;
    gl_FragDepth = depth;
}
