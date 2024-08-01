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

layout(location = 0) in vec4 vertexPosition;
#if defined(USE_COLOR)
layout(location = 1) in vec4 vertexColor;
#endif
#if defined(USE_NORMAL)
layout(location = 2) in vec3 vertexNormal;
#endif
#if defined(USE_TEXCOORD)
layout(location = 3) in vec2 vertexTexCoord;
layout(location = 4) in ivec2 vertexTexIndex;
#endif

uniform mat4 modelViewMatrix_;
uniform mat4 projectionMatrix_;
#if defined(USE_NORMAL)
uniform mat3 normalMatrix_;
#endif
#if !defined(USE_COLOR)
uniform vec4 solidColor_;
#endif

out VertexData {
    // Color and eye-pos attributes always exist.
    vec4 posEye;
    vec4 color;
#if defined(USE_NORMAL)
    vec3 normal;
#endif
#if defined(USE_TEXCOORD)
    vec2  texCoord;
    flat ivec2 texIndex;
#endif
} vertexData;

void main() {
#if defined(USE_COLOR)
    vertexData.color = vertexColor;
#else
    vertexData.color = solidColor_;
#endif
#if defined(USE_NORMAL)
    // Transform normals to eye space.
    vertexData.normal = normalize(normalMatrix_ * vertexNormal);
#endif
#if defined(USE_TEXCOORD)
    vertexData.texCoord = vertexTexCoord;
    vertexData.texIndex = vertexTexIndex;
#endif
    // Transform vertex positions to eye space...
    vertexData.posEye = modelViewMatrix_ * vertexPosition;
    // and project it for gl_Position.
    gl_Position = projectionMatrix_ * vertexData.posEye;
}
