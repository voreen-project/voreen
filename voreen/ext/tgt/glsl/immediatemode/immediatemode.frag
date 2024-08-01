/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#include "immediatemode/mod_immediatelighting.frag"

in vec3 frag_EyePosition;
in vec4 frag_color;
in vec3 frag_EyeNormal;
in vec4 frag_texcoord;

// 0 = no texturing, 1 = 1D, 2 = 2D, 3 = 3D
#if TEXDIM == 1
uniform sampler1D sampler_;
#elif TEXDIM == 2
uniform sampler2D sampler_;
#elif TEXDIM == 3
uniform sampler3D sampler_;
#endif

// lighting uniforms
uniform bool lightingEnabled_;
uniform bool materialColorAmbient_;     // take ambient color from attribute or texture instead from material
uniform bool materialColorDiffuse_;    // take diffuse color from attribute or texture instead from material
uniform bool materialColorSpecular_;   // take specular color from attribute or texture instead from material

out vec4 color;

void main() {
    if (!lightingEnabled_ || materialColorAmbient_ || materialColorDiffuse_ || materialColorSpecular_) {
#if TEXDIM == 0
        color = frag_color;
#elif TEXDIM == 1
        color = texture(sampler_, frag_texcoord.x);
#elif TEXDIM == 2
        color = texture(sampler_, frag_texcoord.xy);
#elif TEXDIM == 3
        color = texture(sampler_, frag_texcoord.xyz);
#endif
    }

    if (lightingEnabled_) {
        vec3 ma = materialColorAmbient_ ? color.xyz : material_.ambientColor_;
        vec3 md = materialColorDiffuse_ ? color.xyz : material_.diffuseColor_;
        vec3 ms = materialColorSpecular_ ? color.xyz : material_.specularColor_;

        color = vec4(phongShading(frag_EyeNormal, frag_EyePosition, ma, md, ms), 1);
    }
}
