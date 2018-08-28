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
#include "modules/mod_filtering.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters textureParameters_;

uniform float saturation_;

void main() {
    // determine normalized fragment position
    vec2 p = gl_FragCoord.xy * screenDimRCP_;

    // lookup input color and depth value (see mod_sampler2d.frag)
    vec4 color = textureLookup2Dnormalized(colorTex_, textureParameters_, p);
    float depth = textureLookup2Dnormalized(depthTex_, textureParameters_, p).x;

    // compute gray value (see mod_filtering.frag) and pass-through depth value
    FragData0 = rgbToGrayScaleSaturated(color, saturation_);
    gl_FragDepth = depth;
}
