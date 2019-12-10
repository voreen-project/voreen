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

// This shader simply renders a screen filling quad without depending on any input.
// Can be invoked using glDrawArrays(GL_POINTS, 0, 1); using an empty vbo
// Obtained from http://stackoverflow.com/questions/2588875/whats-the-best-way-to-draw-a-fullscreen-quad-in-opengl-3-2
layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

out vec2 frag_texcoord;

void main() {
    gl_Position = vec4( 1.0, 1.0, 0.0, 1.0 );
    frag_texcoord = vec2(1.0, 1.0);
    EmitVertex();

    gl_Position = vec4(-1.0, 1.0, 0.0, 1.0 );
    frag_texcoord = vec2(0.0, 1.0);
    EmitVertex();

    gl_Position = vec4( 1.0,-1.0, 0.0, 1.0 );
    frag_texcoord = vec2(1.0, 0.0);
    EmitVertex();

    gl_Position = vec4(-1.0,-1.0, 0.0, 1.0 );
    frag_texcoord = vec2(0.0, 0.0);
    EmitVertex();

    EndPrimitive();
}
