/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "modules/mod_sampler2d.frag"

uniform sampler2D backTex_;
uniform sampler2D frontTex_;
uniform TextureParameters texParameters_;
uniform float threshold_;

in vec2 frag_texcoord_supplied;

/**
 * Adds inverse of back texture to front texture modulated by red component. Used for rendering
 * text that stays visible independent of the background.
 */
void main() {
    vec2 coord = frag_texcoord_supplied;
    vec4 back = textureLookup2Dnormalized(backTex_, texParameters_, coord);
    vec4 front = textureLookup2Dnormalized(frontTex_, texParameters_, coord);

    vec4 backInv = vec4(1.0 - back.rgb, 1.0);
    backInv.r += (backInv.r >= 0.5 ? 1.0 : -1.0) * threshold_;
    backInv.g += (backInv.g >= 0.5 ? 1.0 : -1.0) * threshold_;
    backInv.b += (backInv.b >= 0.5 ? 1.0 : -1.0) * threshold_;

    FragData0 = mix(back, backInv, front.r); // do not use alpha component, as FTGL does not
                                                // generate correct transparency
}
