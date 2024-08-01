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
layout(lines) in;
layout(line_strip, max_vertices=6) out;

uniform mat4 projectionMatrix;

void main() {

    vec3 first = gl_in[0].gl_Position.xyz;
    vec3 second = gl_in[1].gl_Position.xyz;
    vec3 away = vec3(0,0,1);
    vec3 along = normalize(first - second);
    float l = 0.05*distance(first, second);
    vec3 orth = normalize(cross(away, along))*l;


    gl_Position = projectionMatrix*vec4(first,1);
    EmitVertex();
    gl_Position = projectionMatrix*vec4(second,1);
    EmitVertex();

    EndPrimitive();

    gl_Position = projectionMatrix*vec4(first+orth,1);
    EmitVertex();
    gl_Position = projectionMatrix*vec4(first-orth,1);
    EmitVertex();

    EndPrimitive();

    gl_Position = projectionMatrix*vec4(second+orth,1);
    EmitVertex();
    gl_Position = projectionMatrix*vec4(second-orth,1);
    EmitVertex();

    EndPrimitive();
}
