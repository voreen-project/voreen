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
layout(location = 4) in float in_radius;
layout(location = 5) in float in_mass;

uniform mat4 viewMatrix;

out vec4  geom_position;
out vec3  geom_vel;
out float geom_radius;
out float geom_mass;
flat out int   geom_vertexID;


void main(){
    vec4 pos = viewMatrix* vec4(in_vertex, 1);
    geom_position = pos;
    geom_vel = in_vel;
    geom_radius= in_radius;
    geom_vertexID = gl_VertexID;
    geom_mass = in_mass;
    vec3 normVec = in_vertex / abs(in_vertex);
    gl_Position =  pos;
}
