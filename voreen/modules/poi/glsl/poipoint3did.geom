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

#version 330
layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

uniform mat4 P_;
uniform float radius_;

in vec3 vert_color[];
in vec2 vert_size[];
//in uvec2 vert_texture[];

out vec3 geom_color;
out vec2 geom_coord;
out vec2 geom_size;
//flat out uvec2 geom_texture;

void main() {
    float radius = radius_;
    geom_color = vert_color[0];
    geom_size = vert_size[0];
//    geom_texture = vert_texture[0];
    gl_Position = P_*(gl_in[0].gl_Position + radius*vec4(1.0, 1.0, 0.0, 0.0)+radius*0.0*vec4(0.0, 0.0, 1.0, 0.0));
    geom_coord = vec2(1, 1);
    EmitVertex();
    
    geom_color = vert_color[0];
    geom_size = vert_size[0];
//    geom_texture = vert_texture[0];
    gl_Position = P_*(gl_in[0].gl_Position + radius*vec4(1.0, -1.0, 0.0, 0.0)+radius*0.0*vec4(0.0, 0.0, 1.0, 0.0));
    geom_coord = vec2(1, 0);
    EmitVertex();
    
    geom_color = vert_color[0];
    geom_size = vert_size[0];
//    geom_texture = vert_texture[0];
    gl_Position = P_*(gl_in[0].gl_Position + radius*vec4(-1.0, 1.0, 0.0, 0.0)+radius*0.0*vec4(0.0, 0.0, 1.0, 0.0));
    geom_coord = vec2(0, 1);
    EmitVertex();

    geom_color = vert_color[0];
    geom_size = vert_size[0];
//    geom_texture = vert_texture[0];
    gl_Position = P_*(gl_in[0].gl_Position + radius*vec4(-1.0, -1.0, 0.0, 0.0)+radius*0.0*vec4(0.0, 0.0, 1.0, 0.0));
    geom_coord = vec2(0, 0);
    EmitVertex();


    EndPrimitive();
}