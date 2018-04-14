/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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
//layout(line_strip, max_vertices=6) out;

uniform mat4 P_;
uniform mat4 M_;

in VertexData {
    vec3 dir;
    vec3 pos;
    vec3 color;
} VertexIn[];

out VertexData{
    vec2 coord;
    float z;
    vec3 zdir;
    vec3 extdir;
    float radius;
    float depth;
    vec3 dir;
    vec3 color;
} VertexOut;

void main(){
    float radius = 1;
    float z = normalize((M_*vec4(VertexIn[0].dir, 0)).xyz).z;
    float perspective_factor = radius*inversesqrt(1-z*z);

    vec4 pos = M_*vec4(VertexIn[0].pos-perspective_factor*VertexIn[0].dir, 1.0);
    vec4 to = M_*vec4(VertexIn[0].pos+(1+perspective_factor)*VertexIn[0].dir, 1.0);
    //vec4 to = M_*vec4(VertexIn[0].pos+vec3(1, 0, 0), 1.0);
    vec3 dir = to.xyz-pos.xyz;
    vec3 up = normalize(cross(dir, vec3(0, 0, 1)));





    //gl_Position = P_*(pos);
    //EmitVertex();

    //gl_Position = P_*(to);
    //EmitVertex();
    VertexOut.color = VertexIn[0].color;
    VertexOut.z = normalize(dir).z;
    VertexOut.extdir = up;
    VertexOut.zdir = cross(up, dir);
    VertexOut.radius = radius;
    VertexOut.dir = dir;

    dir = vec3(0);
    //dir = -dir/2;

    VertexOut.coord = vec2(-perspective_factor, 1);
    VertexOut.depth = pos.z+radius*up.z-dir.z;
    gl_Position = P_*(pos+vec4(radius*up-dir, 0.0));
    EmitVertex();
    VertexOut.coord = vec2(-perspective_factor, -1);
    VertexOut.depth = pos.z-radius*up.z-dir.z;
    gl_Position = P_*(pos-vec4(radius*up-dir, 0.0));
    EmitVertex();
    VertexOut.coord = vec2(1+perspective_factor, 1);
    VertexOut.depth = to.z+radius*up.z+dir.z;
    gl_Position = P_*(to+vec4(radius*up+dir, 0.0));
    EmitVertex();
    VertexOut.coord = vec2(1+perspective_factor, -1);
    VertexOut.depth = to.z-radius*up.z+dir.z;
    gl_Position = P_*(to-vec4(radius*up+dir, 0.0));
    EmitVertex();

    EndPrimitive();

}
