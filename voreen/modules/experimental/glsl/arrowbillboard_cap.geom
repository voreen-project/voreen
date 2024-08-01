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
//layout(points, max_vertices=6) out;

uniform mat4 P_;
uniform mat4 M_;

in VertexData {
    vec3 dir;
    vec3 pos;
    vec3 color;
} VertexIn[];

out VertexData{
    vec2 coord;
    vec3 normal;
    vec3 color;
} VertexOut;

void main(){
    float radius = 1;
    vec3 pos = VertexIn[0].pos+ VertexIn[0].dir;
    vec3 up = normalize(cross(VertexIn[0].dir, vec3(0, 0, 1)));
    vec3 left = normalize(cross(up, VertexIn[0].dir));



    VertexOut.normal = -VertexIn[0].dir;
    VertexOut.color = VertexIn[0].color;

    VertexOut.coord = vec2(1, 1);
    gl_Position = P_*M_*vec4(pos+radius*up+radius*left, 1);
    EmitVertex();
    VertexOut.coord = vec2(1, -1);
    gl_Position = P_*M_*vec4(pos+radius*up-radius*left, 1);
    EmitVertex();
    VertexOut.coord = vec2(-1, 1);
    gl_Position = P_*M_*vec4(pos-radius*up+radius*left, 1);
    EmitVertex();
    VertexOut.coord = vec2(-1, -1);
    gl_Position = P_*M_*vec4(pos-radius*up-radius*left, 1);
    EmitVertex();

    EndPrimitive();

}
