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



layout(location = 0) in vec3 in_pos;
layout(location = 1) in float in_time;
layout(location = 2) in int in_matrix;
layout(location = 3) in int in_type;


uniform mat4 viewProjectionMatrix_;
uniform mat4 normMats_[101];
uniform vec4 typeColors_[6];

out float frag_time;
out vec4 frag_color;

void main(){
	frag_time = in_time;
	//frag_color = vec4(1.0f,1.0f,1.0f,1.0f);
	frag_color = typeColors_[in_type];
	gl_Position = viewProjectionMatrix_*vec4(in_pos, 1);
}