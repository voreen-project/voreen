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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_normdepth.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex0_;
uniform sampler2D depthTex1_;
uniform TextureParameters texParams_;

uniform vec2 dir_;
uniform float sigma_;
int halfKernelDim_;
float[25] gaussKernel_;

uniform float lambda_;

// computes the Gauss kernel and returns its norm factor
float computeGaussKernel() {
    // compute kernel
    for (int i=0; i<=halfKernelDim_; i++)
        gaussKernel_[i] = exp(-float(i*i)/(2.0*sigma_*sigma_));
    // compute norm
    float norm = 0.0;
    for (int i=1; i<=halfKernelDim_; i++)
        norm += gaussKernel_[i];
    // so far we have just computed norm for one half
    norm = 2.0 * norm + gaussKernel_[0];
    return norm;
}

/***
 * Performs a Gaussian filter on the depth buffer. And returns the
 * filtered value.
 *
 * @fragCoord - screen coordinates of the current fragment
 * @norm - the normalization factor of the kernel
 * @dir - since separability is exploited, the blur direction needs to be specified
 ***/
float applyDepthGaussFilter(in vec2 fragCoord, in float norm, in vec2 dir) {
    float result = 0.0;
    for (int i=-halfKernelDim_; i<=halfKernelDim_; i++) {
        int absi = i;
        if (absi < 0) absi = -absi;
        float curDepth = textureLookup2Dscreen(depthTex0_, texParams_,
                               fragCoord+(dir*i)).x;
        result += curDepth * gaussKernel_[absi];
    }
    result /= norm;
    return result;
}

/***
 * The main method.
 ***/
void main() {

    vec4 fragCoord = gl_FragCoord;
    // initialize and apply Gaussian filter
    halfKernelDim_ = int(2.5 * sigma_);
    float norm = computeGaussKernel();
    float filteredDepth = applyDepthGaussFilter(fragCoord.xy, norm, dir_);

    if (dir_.y != 1.0) {
        // only write out filtered depth value
        FragData0 = vec4(1.0);
        gl_FragDepth = filteredDepth;
    } else {
        // we are in the second vertical pass and have to modulate the color
        float curDepth = textureLookup2Dscreen(depthTex1_, texParams_, fragCoord.xy).x;
        float deltaD = normDepth(filteredDepth) - normDepth(curDepth);
        // apply depth darkening
        float deltaDNeg = 0.0;
        if (deltaD < 0.0) deltaDNeg = deltaD;
        vec4 curColor = textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy);
        curColor.rgb += deltaDNeg * lambda_;
        // write out modulated color and original depth value
        FragData0 = curColor;
        gl_FragDepth = curDepth;
    }
}
