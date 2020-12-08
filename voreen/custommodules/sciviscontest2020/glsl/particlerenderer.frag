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

uniform TF_SAMPLER_TYPE transFuncTex_;
uniform TransFuncParameters transFuncParam_;

flat in vec4 values;

uniform vec2 domain;
uniform ivec2 timesteps;

void main() {	
    if( values.y < float( timesteps.x ) - 0.5 || values.y > float( timesteps.y ) + 0.5 ) discard;
	
    float value = clamp( ( values.x - domain.x ) / ( domain.y - domain.x ), 0.0, 1.0 );
    vec4 color = texture( transFuncTex_, value );
    if( color.a == 0.0 ) discard;
	
    FragData0 = color;
}
