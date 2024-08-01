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

//include shader libraries (shader modules)

//defines and functions for 2D textures
#include "modules/mod_sampler2d.frag"
//defines and functions for filtering
#include "modules/mod_filtering.frag"

//uniforms for the color and depth textures
uniform sampler2D colorTex_;
uniform sampler2D depthTex_;

//struct defined in mod_sampler2d.frag containing all needed parameters
uniform TextureParameters textureParameters_;

//the saturation (user-defined by a property)
uniform float saturation_;

/**
 * Main function of the shader. It takes the color texture passed as a uniform
 * and modifies the saturation.
 */
void main() {
    // look up input color and depth value (see mod_sampler2d.frag)
    vec4 color = textureLookup2Dscreen(colorTex_, textureParameters_, gl_FragCoord.xy);
    float depth = textureLookup2Dscreen(depthTex_, textureParameters_, gl_FragCoord.xy).z;

    // compute gray value (see mod_filtering.frag) and pass-through depth value
    FragData0 = rgbToGrayScaleSaturated(color, saturation_);
    gl_FragDepth = depth;
}
