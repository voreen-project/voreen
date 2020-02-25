/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

layout(lines) in;
layout(line_strip, max_vertices = 2) out;

in vData
{
    vec3 position;
    vec3 velocity;
    float radius;
    float time;
}vertex[];

out fData
{
    vec3 position;
    vec3 velocity;
    float radius;
    float time;
}frag;

void main()
{
    frag.position = vertex[0].position;
    frag.velocity = vertex[0].velocity;
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    frag.position = vertex[1].position;
    frag.velocity = vertex[1].velocity;
    gl_Position = gl_in[1].gl_Position;
    EmitVertex();

    EndPrimitive();
}
