/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

uniform sampler2D colorTex0_;
uniform TextureParameters textureParameters0_;
uniform sampler2D colorTex1_;
uniform sampler2D depthTex1_;
uniform TextureParameters textureParameters1_;
uniform int option_;
uniform float threshold_ = 0.2;

void main() {
    vec2 fragCoord = gl_FragCoord.xy;

    // fetch input textures
    vec4 front = textureLookup2Dscreen(colorTex0_, textureParameters0_, fragCoord);
    vec4 back = textureLookup2Dscreen(colorTex1_, textureParameters1_, fragCoord);
    float fragDepth = textureLookup2Dscreen(depthTex1_, textureParameters1_, fragCoord).x;

    vec4 result;
    vec4 backInv;

    if(option_ == 0){
        backInv = vec4(1.0-back.r, 1.0-back.g, 1.0-back.b, 1.0);
        backInv.r = backInv.r >= 0.5
            ? max(backInv.r, backInv.r + threshold_)
            : min(backInv.r, backInv.r - threshold_);
        backInv.g = backInv.g >= 0.5
            ? max(backInv.g, backInv.g + threshold_)
            : min(backInv.g, backInv.g - threshold_);
        backInv.b = backInv.b >= 0.5
            ? max(backInv.b, backInv.b + threshold_)
            : min(backInv.b, backInv.b - threshold_);
        result = back*(1.0-front.a)+backInv*front.a;
    } else
        result = back*(1.0-front.a)+front;

    FragData0 = result;
    gl_FragDepth = fragDepth;
}
