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
uniform int kernelRadius_;
uniform int mode_; //0: dilation, 1: erosion
uniform int shape_; //0: square, 1: circle

void main() {
    vec2 fragCoord = gl_FragCoord.xy;

    if (mode_ == 0 && shape_ == 0) FragData0 = dilation(colorTex_, textureParameters_, fragCoord, kernelRadius_);
    else if (mode_ == 0 && shape_ == 1) FragData0 = dilation_circle(colorTex_, textureParameters_, fragCoord, kernelRadius_);
    else if (mode_ == 1 && shape_ == 0) FragData0 = erosion(colorTex_, textureParameters_, fragCoord, kernelRadius_);
    else if (mode_ == 1 && shape_ == 1) FragData0 = erosion_circle(colorTex_, textureParameters_, fragCoord, kernelRadius_);

    gl_FragDepth = textureLookup2Dscreen(depthTex_, textureParameters_, fragCoord).x;
}

