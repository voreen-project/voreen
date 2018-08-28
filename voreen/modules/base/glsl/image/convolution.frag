/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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
#include "modules/mod_filtering.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters textureParameters_;

uniform sampler2D filter_;
uniform TextureParameters filterParameters_;

uniform int filterSize_;

void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;

    vec4 result = vec4(0.0);
    float ffs = filterSize_;
    float d = 0.0;

    for(int x=-filterSize_; x<=filterSize_; ++x) {
        for(int y=-filterSize_; y<=filterSize_; ++y) {
            vec2 filterPos = ((vec2(x,y)/ffs)+1.0)*0.5;
            vec4 filt = textureLookup2Dnormalized(filter_, filterParameters_, filterPos);
            filt *= 2.0;
            filt -= 1.0;
            d += filt.r;
            vec4 orig = textureLookup2Dnormalized(colorTex_, textureParameters_, (gl_FragCoord.xy+vec2(x,y)) * screenDimRCP_);
            result += filt*orig;
        }
    }
    /*result /= d;*/

    FragData0 = result;
    gl_FragDepth = textureLookup2Dnormalized(depthTex_, textureParameters_, p).x;
}
