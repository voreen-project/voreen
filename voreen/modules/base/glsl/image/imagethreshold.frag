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
#include "modules/mod_filtering.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters texParams_;

uniform float lowerThreshold_;
uniform float upperThreshold_;
uniform vec4 lowerMaskColor_;
uniform vec4 upperMaskColor_;

void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec4 fragColor = textureLookup2Dnormalized(colorTex_, texParams_, p);
    float gray = rgbToGrayScale(fragColor).x;

    if (gray <= lowerThreshold_)
        fragColor = lowerMaskColor_;
    else if (gray > upperThreshold_)
        fragColor = upperMaskColor_;

    FragData0 = fragColor;
    gl_FragDepth = textureLookup2Dnormalized(depthTex_, texParams_, p).x;
}
