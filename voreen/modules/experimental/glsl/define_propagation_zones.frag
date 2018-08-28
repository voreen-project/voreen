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

uniform sampler2D entryPoints_;                      // ray entry points
uniform TextureParameters entryParameters_;
uniform sampler2D exitPoints_;                        // ray exit points
uniform TextureParameters exitParameters_;

uniform vec2 projectedLight_;

/***
 * Defines correct propagation zone, and represented zone color
 ***/
vec4 definePropagationZone(vec2 p) {
    vec2 dir = p - projectedLight_;
    vec4 result = vec4(0.0);

    if(abs(dir.x) > abs(dir.y)){
        if(dir.x < 0)
            result = vec4(1.0, 0.0, 0.0, 1.0);
        else
            result = vec4(0.0, 0.0, 1.0, 1.0);
    }
    else{
        if(dir.y < 0)
            result = vec4(0.0, 1.0, 0.0, 1.0);
        else
            result = vec4(1.0, 1.0, 0.0, 1.0);
    }

    return result;
}

/***
 * The main method.
 ***/
void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;

    vec3 frontPos = textureLookup2Dnormalized(entryPoints_, entryParameters_, p).rgb;
    vec3 backPos = textureLookup2Dnormalized(exitPoints_, exitParameters_, p).rgb;

    //determine whether the ray has to be casted
    if (frontPos == backPos) {
        //background need no raycasting
        discard;
    } else {
        //fragCoords are lying inside the boundingbox
        FragData0 = definePropagationZone(p);
    }
}
