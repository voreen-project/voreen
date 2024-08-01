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

layout(location=0)in vec3 vert_position;
layout(location=3)in vec4 vert_texcoord;

uniform mat4 textureMatrix_;
#if defined(SLICE_TEXTURE_MODE_2D)
uniform mat4 toSliceCoordMatrix;
#endif

out vec4 frag_volcoord;
#if defined(SLICE_TEXTURE_MODE_2D)
out vec2 frag_slicecoord;
#endif

void main() {
    frag_volcoord = textureMatrix_ * vert_texcoord;
    #if defined(SLICE_TEXTURE_MODE_2D)
    frag_slicecoord = (toSliceCoordMatrix * frag_volcoord).xy;
    #endif

    // set out vertex
    gl_Position = modelViewProjectionMatrixStack_ * vec4(vert_position, 1.0);
}
