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

#include "modules/mod_shading.frag"

struct FragStruct {
  vec4 posEye;
  vec4 color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
  vec3 normal;
#endif
};

in FragStruct fragData;

uniform bool highlight_;

uniform mat4 viewMatrix_;

void main() {
    vec4 result = fragData.color;
#if defined(TRIANGLE_VEC4_VEC3) || defined(TRIANGLE_VEC3)
    result.xyz = APPLY_SHADING(fragData.normal, fragData.posEye.xyz, lightSource_.position_, vec3(0.0), fragData.color.rgb, fragData.color.rgb, fragData.color.rgb);
#endif
    if(highlight_){
        result.xyz = result.xyz * 1.2;
    }
    //result = vec4(FragmentIn.normal, 1.0);
    FragData0 = result;
}
