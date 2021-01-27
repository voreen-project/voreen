/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#version 430 core

// vertex shader input
layout(location = 0) in vec4 position;
#if defined(USE_COLOR)
layout(location = 1) in vec4 color;
#endif
#if defined(USE_NORMAL)
layout(location = 2) in vec3 normal;
#endif
#if defined(USE_TEXCOORD)
layout(location = 3) in vec2 texCoord;
layout(location = 4) in ivec2 texIndex;
#endif

// model transformation
layout(location = 0) uniform mat4 modelMatrix;

// vertex shader output
out VertexData {
    vec4 position;
#if defined(USE_COLOR)
    vec4 color;
#endif
#if defined(USE_NORMAL)
    vec3 normal;
#endif
#if defined(USE_TEXCOORD)
    vec2  texCoord;
    ivec2 texIndex;
#endif
} vout;

void main(){
    //just pass through the vertex data in model space
    vout.position = position;
#if defined(USE_COLOR)    
    vout.color = color;
#endif
#if defined(USE_NORMAL)
    vout.normal = normal;
#endif
#if defined(USE_TEXCOORD)
    vout.texCoord = texCoord;
    vout.texIndex = texIndex;
#endif

    //compute vertex position in world space for geometry shader
    gl_Position = modelMatrix * position;
}
