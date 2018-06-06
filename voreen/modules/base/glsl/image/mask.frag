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

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters texParams_;
uniform sampler2D colorTexMask_;
uniform TextureParameters texParamsMask_;

uniform vec4 maskColor_;

void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec4 color = textureLookup2Dnormalized(colorTex_, texParams_, p);
    float depth = textureLookup2Dnormalized(depthTex_, texParams_, p).x;
    vec4 maskValue = textureLookup2Dnormalized(colorTexMask_, texParamsMask_, p);
    if (length(maskValue - maskColor_) < 0.01) {
        FragData0 = vec4(0.0);
        gl_FragDepth = 1.0;
    }
    else {
        FragData0 = color;
        gl_FragDepth = depth;
    }
}
