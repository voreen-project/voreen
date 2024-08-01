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
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4 out_color;

in float frag_time;
in vec4 frag_color;

uniform TransFuncParameters transferFunc_;
uniform sampler1D transferFuncTex_;
uniform vec2 interval_;
uniform float alphaFactor_;


void main(){
	if (frag_time < interval_.x  || frag_time > interval_.y) discard;

	float intensity = applyTF(transferFunc_, transferFuncTex_, frag_time).a;
	out_color = frag_color;
	//out_color = applyTF(transferFunc_, transferFuncTex_, frag_time);
	out_color.a *= alphaFactor_ * intensity;
}