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
layout(triangle_strip, max_vertices=4) out;

in float geom_radius[2];

uniform mat4 projectionMatrix;
uniform float baseWidth;
uniform float spacing;

void emitEdgeVertex(vec3 pos) {
    gl_Position = projectionMatrix * vec4(pos,1);
    EmitVertex();
}

void main() {
    vec3 p0 = gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_in[1].gl_Position.xyz;

    vec3 alongLine = normalize(p1 - p0);
    vec3 orthLine = normalize(cross(alongLine, vec3(0,0,1)))*baseWidth;

    emitEdgeVertex(p0+orthLine+alongLine*(geom_radius[0]*spacing));
    emitEdgeVertex(p0-orthLine+alongLine*(geom_radius[0]*spacing));
    emitEdgeVertex(p1         -alongLine*(geom_radius[1]*spacing));

    EndPrimitive();
}
