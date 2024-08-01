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

// textures storing the overlay rendering
uniform sampler2D colorTexMe_;
uniform sampler2D depthTexMe_;
uniform TextureParameters textureParametersMe_;

// textures of incoming rendering (inport)
uniform sampler2D colorTexIn_;
uniform sampler2D depthTexIn_;
uniform TextureParameters textureParametersIn_;

void main() {
    vec2 p = gl_FragCoord.xy;
    vec4 colorMe = textureLookup2Dnormalized(colorTexMe_, textureParametersMe_, p*screenDimRCP_);
    float depthMe = textureLookup2Dnormalized(depthTexMe_, textureParametersMe_, p*screenDimRCP_).x;
    vec4 colorIn = textureLookup2Dnormalized(colorTexIn_, textureParametersIn_, p*screenDimRCP_);
    float depthIn = textureLookup2Dnormalized(depthTexIn_, textureParametersIn_, p*screenDimRCP_).x;

    vec4 fragColor;
    float fragDepth;

    if (colorMe.a > 0.0) {
        fragColor = colorMe;
        fragDepth = depthMe;
    }
    else {
        fragColor = colorIn;
        fragDepth = depthIn;
    }

    FragData0 = fragColor;
    gl_FragDepth = fragDepth;
}
