/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

uniform float opacity_;
uniform sampler2D puckTex_;
uniform TextureParameters texParams_;
//uniform sampler2D symbTex_;
uniform sampler2D innerTex_;

uniform vec4 colorMod_;
uniform vec2 ll_;
uniform vec2 ur_;

uniform float greyOutFactor_;

void main() {

    vec2 fragCoord = clamp((gl_FragCoord.xy - ll_) / (ur_ - ll_), 0.0, 1.0);

    vec4 puckColor = textureLookup2Dnormalized(puckTex_, texParams_, fragCoord);
    vec4 innerColor = textureLookup2Dnormalized(innerTex_, texParams_, fragCoord);

    vec4 backgroundColor = vec4(colorMod_.rgb, colorMod_.a * innerColor.a);
    backgroundColor = backgroundColor * backgroundColor.a + (1-backgroundColor.a) * puckColor;

    backgroundColor.a = backgroundColor.a * opacity_;

    vec4 color = mix(backgroundColor, vec4(vec3(0.3), backgroundColor.a), greyOutFactor_);

    FragData0 = color;

}
