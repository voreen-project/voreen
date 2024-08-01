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

#version 330



layout(location = 0) in vec3 in_vertex;
layout(location = 1) in vec3 in_vel;
layout(location = 2) in vec3 in_angular_momenta;
layout(location = 3) in float in_spin_parameter;
layout(location = 4) in float in_radius;

uniform mat4 viewMatrix;

out vec4  geom_position;
out vec3  geom_vel;
out vec3  geom_angular_momenta;
out float geom_spin_parameter;
out float geom_radius;

void main(){
	vec4 pos = viewMatrix* vec4(in_vertex, 1);
	geom_position = pos;
	geom_vel = in_vel;
	geom_angular_momenta =(viewMatrix* vec4(in_angular_momenta,0)).xyz;
	geom_spin_parameter= in_spin_parameter;
	geom_radius= in_radius;	

	gl_Position =  pos ;
}