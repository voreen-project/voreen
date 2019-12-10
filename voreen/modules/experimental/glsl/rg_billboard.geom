/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

//#extension GL_EXT_geometry_shader4: enable

layout(points) in;
layout(triangle_strip, max_vertices=4) out;

uniform float radius;

varying out vec3 vertex_light_position;
varying out vec4 eye_position;

void main() {
    //vertex_light_position = vec3(0.f,0.f,1.f);//normalize(gl_LightSource[0].position.xyz);
    eye_position = gl_PositionIn[0];

    float r = radius;

    gl_FrontColor = gl_FrontColorIn[0];

    gl_TexCoord[0].st = vec2(1.0, 1.0);
    gl_Position = gl_PositionIn[0];
    gl_Position += vec4(r, r, 0, 0);
    gl_Position = gl_ProjectionMatrix * gl_Position;
    EmitVertex();

    gl_TexCoord[0].st = vec2(-1.0, 1.0);
    gl_Position = gl_PositionIn[0];
    gl_Position += vec4(-r, r, 0, 0);
    gl_Position = gl_ProjectionMatrix * gl_Position;
    EmitVertex();

    gl_TexCoord[0].st = vec2(1.0, -1.0);
    gl_Position = gl_PositionIn[0];
    gl_Position += vec4(r, -r, 0, 0);
    gl_Position = gl_ProjectionMatrix * gl_Position;
    EmitVertex();

    gl_TexCoord[0].st = vec2(-1.0, -1.0);
    gl_Position = gl_PositionIn[0];
    gl_Position += vec4(-r, -r, 0, 0);
    gl_Position = gl_ProjectionMatrix * gl_Position;
    EmitVertex();
    EndPrimitive();
}
