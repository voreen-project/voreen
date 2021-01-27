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
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

uniform sampler2D backTex_;
uniform sampler2D frontTex_;

uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;

const float Samplings = 250.0;


/***
 * returns the coord of the first hit.
 ***/
vec4 hitpoints(in vec3 first, in vec3 last) {

    vec3 direction = last - first;
    int steps = int(floor(Samplings * length(direction)));
    vec3 diff1 = direction / float(steps);
    vec4 value;
    vec4 result = vec4(0.0);

    for (int i=0; i<steps; i++) {
        value.a = textureLookup3D(volume_, volumeStruct_ ,first).a;
        first += diff1;
        if (value.a != 0.0) {
            result=vec4(first,1.0);
            i = steps;
        }
    }
    return result;
}

/***
 * The main method.
 ***/
void main() {
    vec4 fragCoord = vec4(gl_FragCoord.xy, 0.0, 0.0);
    vec3 frontPos = textureLookup2D(frontTex_, fragCoord.xy).rgb;
    vec3 backPos = textureLookup2D(backTex_, fragCoord.xy).rgb;

    FragData0 = hitpoints(frontPos, backPos);
}
