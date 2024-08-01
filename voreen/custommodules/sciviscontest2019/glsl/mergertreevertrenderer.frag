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

out vec4 out_color;

in vec4 frag_color;
in vec2 frag_coord;
in vec3 frag_pos;
in float frag_radius;

uniform sampler2D colorTex_;

void main() {
    if (dot(frag_coord, frag_coord) > 1)
        discard;
#ifdef USE_TEXTURES
    vec4 textureColor = texture(colorTex_, 0.5f*(vec2(1.0f, 1.0f)+frag_coord));
    out_color = textureColor.x*frag_color;
#else
    float dz= sqrt(1 - frag_coord.x*frag_coord.x - frag_coord.y*frag_coord.y);
    float dcont = max(0.0,dz);
    out_color = vec4(pow(dcont,1)*frag_color.xyz, 1);
    gl_FragDepth = 1-(dz*frag_radius*0.1f);
#endif
}
