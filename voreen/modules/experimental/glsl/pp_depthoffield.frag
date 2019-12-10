/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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
uniform sampler2D depthTex_;
uniform TextureParameters textureParameters_;

uniform float depthFocus_;
uniform float maxSigma_;

float sigmaX_;
float sigmaY_;
int halfKernelDimX_;
int halfKernelDimY_;
float[25] gaussKernelX_;
float[25] gaussKernelY_;

// computes the Gauss kernel and returns its norm factor
float computeGaussKernel() {
    // compute kernel
    for (int x=0; x<=halfKernelDimX_; x++)
        gaussKernelX_[x] = exp(-float(x*x)/(2.0*sigmaX_*sigmaX_));
    for (int y=0; y<=halfKernelDimY_; y++)
        gaussKernelY_[y] = exp(-float(y*y)/(2.0*sigmaY_*sigmaY_));
    // compute norm
    float norm = 0.0;
    for (int x=0; x<=halfKernelDimX_; x++) {
        for (int y=1; y<=halfKernelDimY_; y++)
            norm += gaussKernelX_[x]*gaussKernelY_[y];
    }
    // so far we have just computed norm for one quarter
    norm = 4.0 * norm + gaussKernelX_[0]*gaussKernelY_[0];
    return norm;
}

/***
 * Performs a Gaussian filter on the color buffer. And returns the
 * filtered value.
 *
 * @fragCoord - screen coordinates of the current fragment
 * @norm - the normalization factor of the kernel
 ***/
vec4 applyColorGaussFilter(in vec2 fragCoord, in float norm) {
    vec4 result = vec4(0.0);
    for (int x=-halfKernelDimX_; x<=halfKernelDimX_; x++) {
        int absx = x;
        if (absx < 0) absx = -absx;
        for (int y=-halfKernelDimY_; y<=halfKernelDimY_; y++) {
            int absy = y;
            if (absy < 0) absy = -absy;
            vec4 curColor = textureLookup2Dscreen(colorTex_, textureParameters_,
                                                  vec2(fragCoord.x+float(x), fragCoord.y+float(y)));
            result += curColor * gaussKernelX_[absx]*gaussKernelY_[absy];
        }
    }
    result /= norm;
    return result;
}

/***
 * The main method.
 ***/
void main() {
    vec2 fragCoord = gl_FragCoord.xy;

    float curDepth = textureLookup2Dscreen(depthTex_, textureParameters_, fragCoord).z;
    float curDepthNorm = normDepth(curDepth);

    float realDepth = clamp(abs(depthFocus_ - curDepthNorm), 0.0, 0.5) * 2.0;
    sigmaX_ = realDepth * maxSigma_;
    sigmaY_ = realDepth * maxSigma_;
    halfKernelDimX_ = int(2.5 * sigmaX_);
    halfKernelDimY_ = int(2.5 * sigmaY_);
    float norm = computeGaussKernel();
    FragData0 = applyColorGaussFilter(fragCoord, norm);
    gl_FragDepth = curDepth;
}
