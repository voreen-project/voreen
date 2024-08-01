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

layout(location = 0) in vec3 vert_pos;
layout(location = 1) in vec3 vert_normal;
 
out VertexData {
    vec3 normal;    
	vec4 color;
} VertexOut;

uniform mat4 MVP_;
uniform mat4 MV_;
uniform mat4 P_;
uniform mat3 MVtinv_;

void main() {    
	VertexOut.color = vec4(0.5*normalize(vert_normal)+vec3(0.5), 1.0);
    VertexOut.normal = MVtinv_*vert_normal;
    gl_Position = MV_*vec4(vert_pos, 1.0);
}
