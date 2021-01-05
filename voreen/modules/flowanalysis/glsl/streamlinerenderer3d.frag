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

#include "modules/mod_transfunc.frag"
#include "modules/mod_shading.frag"

in vData {
    vec3 position;
    vec3 velocity;
    float radius;
    float time;
} frag;

uniform float timeWindowStart_;
uniform float timeWindowSize_;

#ifdef COLOR_VELOCITY
    uniform TF_SAMPLER_TYPE transFuncTex_;      //< defined in generate header
    uniform TransFuncParameters transFuncParam_;
#elif defined (COLOR_DIRECTION)
    uniform mat4 colorRotationMatrix_;

    //direction is not null
    vec4 applyDirection(vec3 direction) {
        return vec4((colorRotationMatrix_ * vec4(normalize(direction),1.0)).xyz/vec3(2.0) + vec3(0.5),1.0);
    }
#else
    #error No color mode has been set
#endif


void main() {

    if(frag.time < timeWindowStart_ || frag.time > timeWindowStart_ + timeWindowSize_)
        discard;

    vec4 color = vec4(1.0);
#ifdef COLOR_VELOCITY
    color =  applyTF(transFuncParam_, transFuncTex_, vec4(length(frag.velocity)));
#elif defined (COLOR_DIRECTION)
    if(length(frag.velocity) == 0.0)
        discard;
    color = applyDirection(frag.velocity);
#else
    #error No color mode has been set
#endif

    if(color.a == 0.0)
        discard;

    FragData0 = color;
}
