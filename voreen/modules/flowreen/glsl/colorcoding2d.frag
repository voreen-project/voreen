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

#include "mod_colorcoding.frag"

uniform sampler2D vectorField2D_;

float getFlowMagnitude(const vec2 r) {
    vec2 v = texture(vectorField2D_, r).xy;
    if (v != vec2(0.0)) {
        v = (v * (maxValue_ - minValue_)) + minValue_;
        return length(v);
    }
    return 0.0;
}

vec4 getColorFromFlowDirection(const vec2 r) {
    vec2 v = texture(vectorField2D_, r).xy;
    if (v != vec2(0.0))
        v = normalize((v * (maxValue_ - minValue_)) + minValue_);
    return vec4(abs(v.x), 0.0, abs(v.y), 1.0);
}

void main() {
#if COLOR_MODE == COLOR_MODE_DIRECTION
    FragData0 = getColorFromFlowDirection(gl_TexCoord[0].st);
#elif COLOR_MODE == COLOR_MODE_MONOCHROME
    FragData0 = clamp(color_, vec4(0.0), vec4(1.0));
#else
    // look up the flow at the position and calculate the color
    vec2 direction = texture(vectorField2D_, gl_TexCoord[0].st).xy;
    if(direction != vec2(0.0))
        direction = (direction * (maxValue_ - minValue_)) + minValue_;
    FragData0 = getFlowColor(vec3(direction, 0.0));
#endif
}
