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

layout(points) in;
layout(line_strip, max_vertices=200) out;

in vec4  geom_position[1];
in vec3  geom_vel[1];
in vec3  geom_angular_momenta[1];
in float geom_spin_parameter[1];
in float geom_radius[1];

out float frag_spin;

uniform mat4 projectionMatrix;


void main() {

    frag_spin =geom_spin_parameter[0];
    float r =geom_radius[0];;

    vec3 angularMomenta = geom_angular_momenta[0];	
    vec4 pos = geom_position[0];  //introduce a single vertex at the origin
    vec3 nonCoLinearAngularMomenta = vec3(-angularMomenta.z, angularMomenta.x, angularMomenta.y);
    vec3 sinBase = normalize(cross(angularMomenta, nonCoLinearAngularMomenta))* r*1.2;
    vec3 cosBase = normalize(cross(angularMomenta, sinBase))*r*1.2;

    for(float i = 0; i < 6.38 ; i+=0.1)  //generate vertices at positions on the circumference from 0 to 2*pi
    {
	vec3 newPos= vec3(pos.xyz+cosBase*cos(i)+sinBase*sin(i));
	gl_Position = projectionMatrix* vec4(newPos,1);
	EmitVertex();
    }

    EndPrimitive();


}
