/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

layout(location=0)in vec4 vert_position;
layout(location=1)in vec4 vert_color;
layout(location=2)in vec3 vert_normal;
layout(location=3)in vec4 vert_texcoord;

#ifdef CLIPPING_ENABLED
uniform vec4 clip_planes[NUM_CLIP_PLANES];
#endif

out vec3 frag_EyePosition;
out vec4 frag_color;
out vec3 frag_EyeNormal;
out vec4 frag_texcoord;
#ifdef CLIPPING_ENABLED
out float gl_ClipDistance[NUM_CLIP_PLANES];
#endif

/**
 * Simply pass through the provided vertex data.
 */
void main() {
    gl_Position = modelViewProjectionMatrixStack_ * vert_position;
    vec4 eyePosition = modelViewMatrixStack_ * vert_position;
    frag_EyePosition = eyePosition.xyz;
    frag_color = vert_color;
    frag_EyeNormal = mat3(transpose(modelViewMatrixInverseStack_))* vert_normal;
    frag_texcoord = vert_texcoord;
#ifdef CLIPPING_ENABLED
    for(int i = 0; i < NUM_CLIP_PLANES; ++i) {
        gl_ClipDistance[i] = dot(eyePosition, clip_planes[i]);
    }
#endif
}
