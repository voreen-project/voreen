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
layout(triangle_strip, max_vertices=4) out;

in vec4  geom_position[1];
in float geom_bodyValue[1];
in vec3  geom_vel[1];
in float geom_radius[1];

out vec2  frag_coord;
out float frag_bodyValue;
out vec3  frag_vel;
out vec3  position;

uniform mat4 projectionMatrix;

void main() {
    vec3 position= geom_position[0].xyz;
    float r = geom_radius[0];
    vec4 up = vec4(0, 1, 0, 0);
    vec4 right = vec4(1, 0, 0, 0);

    frag_coord     = vec2(1.0, 1.0);
    frag_bodyValue = geom_bodyValue[0];
    frag_vel 	   = geom_vel[0];
    gl_Position    = projectionMatrix*(r*right+r*up+geom_position[0]);
    EmitVertex();

    frag_coord     = vec2(-1.0, 1.0);
    frag_bodyValue = geom_bodyValue[0];
    frag_vel 	   = geom_vel[0];
    gl_Position    = projectionMatrix*(-r*right+r*up+geom_position[0]);
    EmitVertex();

    frag_coord     = vec2(1.0, -1.0);
    frag_bodyValue = geom_bodyValue[0];
    frag_vel 	   = geom_vel[0];
    gl_Position    = projectionMatrix*(r*right-r*up+geom_position[0]);
    EmitVertex();

    frag_coord     = vec2(-1.0, -1.0);
    frag_bodyValue = geom_bodyValue[0];
    frag_vel 	   = geom_vel[0];
    gl_Position    = projectionMatrix*(-r*right-r*up+geom_position[0]);
    EmitVertex();

    EndPrimitive();

}
