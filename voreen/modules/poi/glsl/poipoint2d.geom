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
layout(triangle_strip, max_vertices = 4) out;


uniform mat4 screenToNDC_;
uniform vec3 basex_;
uniform vec3 basey_;
uniform float radius_;

in vec3  vert_color[];
in float vert_active[];
in float vert_id[];

out vec3  geom_color;
out vec2  geom_coord;
out float geom_depth;
out float geom_active;
out float geom_id;

void main() {
    float radius = radius_;
    vec4 basex = vec4(1, 0, 0, 0);
    vec4 basey = vec4(0, 1, 0, 0);
    
    geom_color = vert_color[0];
    geom_active = vert_active[0];
    gl_Position = screenToNDC_*(gl_in[0].gl_Position + radius*basex + radius*basey);
    geom_depth = gl_in[0].gl_Position.z;
    geom_coord = vec2(1, 1);
    geom_id = vert_id[0];
    EmitVertex();
    
    geom_color = vert_color[0];
    geom_active = vert_active[0];
    gl_Position = screenToNDC_*(gl_in[0].gl_Position + radius*basex - radius*basey);
    geom_depth = gl_in[0].gl_Position.z;
    geom_coord = vec2(1, -1);
    geom_id = vert_id[0];
    EmitVertex();
    
    geom_color = vert_color[0];
    geom_active = vert_active[0];
    gl_Position = screenToNDC_*(gl_in[0].gl_Position - radius*basex + radius*basey);
    geom_depth = gl_in[0].gl_Position.z;
    geom_coord = vec2(-1, 1);
    geom_id = vert_id[0];
    EmitVertex();

    geom_color = vert_color[0];
    geom_active = vert_active[0];
    gl_Position = screenToNDC_*(gl_in[0].gl_Position - radius*basex - radius*basey);
    geom_depth = gl_in[0].gl_Position.z;
    geom_coord = vec2(-1, -1);
    geom_id = vert_id[0];
    EmitVertex();


    EndPrimitive();
}