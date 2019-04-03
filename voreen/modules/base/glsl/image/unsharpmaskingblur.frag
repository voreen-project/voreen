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

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters textureParameters_;

uniform vec2 dir_;
uniform int kernelRadius_;
uniform float[26] gaussKernel_;
uniform float norm_;
uniform vec4 channelWeights_;

/***
 * Performs a Gaussian filter and returns the filtered value.
 *
 * @fragCoord - screen coordinates of the current fragment
 * @norm - the normalization factor of the kernel
 * @dir - since separability is exploited, the blur direction needs to be specified
 * @channelWeights - vec4 with 1 for a color channel the filter should be applied to, otherwise 0 (values in between blend between filtered and original image)
 ***/
vec4 applyColorGaussFilter(in vec2 fragCoord, in float norm, in vec2 dir, in vec4 channelWeights) {
    vec4 result = vec4(0.0);
    for (int i=-kernelRadius_; i<=kernelRadius_; i++) {
        int absi = i;
        if (absi < 0) absi = -absi;
        vec4 curColor = textureLookup2Dscreen(colorTex_, textureParameters_,
                               fragCoord+(dir*i));
        result +=  curColor * gaussKernel_[absi];
    }
    result /= norm;

    //change channels according to their weights
    vec4 final = result * channelWeights + textureLookup2Dscreen(colorTex_, textureParameters_,fragCoord)*(vec4(1.0,1.0,1.0,1.0) - channelWeights);
    return final;
}

/***
 * The main method.
 ***/
void main() {
    vec2 fragCoord = gl_FragCoord.xy;

    // apply Gaussian
    FragData0 = applyColorGaussFilter(fragCoord, norm_, dir_, channelWeights_);

    // retain depth values from the original
    float depth;
    depth = textureLookup2Dscreen(depthTex_, textureParameters_, fragCoord).x;
    gl_FragDepth = depth;
 }
