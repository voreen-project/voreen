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

#include "gradientfilter.h"

#include "slicereader.h"

namespace voreen {

GradientFilter::GradientFilter(GradientType gradientType, const tgt::vec3& spacing, const SamplingStrategy<ParallelFilterValue1D>& samplingStrategy, const std::string& sliceBaseType)
    : ParallelVolumeFilter<ParallelFilterValue1D, ParallelFilterValue3D>(1, samplingStrategy, sliceBaseType)
    , gradientType_(gradientType)
    , spacing_(spacing)
{
    switch (gradientType_) {
    case VOG_CENTRAL_DIFFERENCE:
        gradientFunction_ = std::bind(&GradientFilter::calcGradientCentralDifferences, this, std::placeholders::_1, std::placeholders::_2);
        break;
    case VOG_LINEAR_REGRESSION:
        gradientFunction_ = std::bind(&GradientFilter::calcGradientLinearRegression, this, std::placeholders::_1, std::placeholders::_2);
        break;
    case VOG_SOBEL:
        gradientFunction_ = std::bind(&GradientFilter::calcGradientSobel, this, std::placeholders::_1, std::placeholders::_2);
        break;
    default:
        tgtAssert(false, "Unimplemented gradient type");
    }
}

ParallelFilterValue3D GradientFilter::getValue(const Sample& sample, const tgt::ivec3& pos) const {
    return gradientFunction_(sample, pos);
}

GradientType GradientFilter::getGradientType() const {
    return gradientType_;
}

ParallelFilterValue3D GradientFilter::calcGradientCentralDifferences(const Sample& sample, const tgt::ivec3& pos) {

    ParallelFilterValue1D v0 = sample(pos + tgt::ivec3(1, 0, 0));
    ParallelFilterValue1D v1 = sample(pos + tgt::ivec3(0, 1, 0));
    ParallelFilterValue1D v2 = sample(pos + tgt::ivec3(0, 0, 1));

    ParallelFilterValue1D v3 = sample(pos - tgt::ivec3(1, 0, 0));
    ParallelFilterValue1D v4 = sample(pos - tgt::ivec3(0, 1, 0));
    ParallelFilterValue1D v5 = sample(pos - tgt::ivec3(0, 0, 1));

    tgt::vec3 gradient(v3 - v0, v4 - v1, v5 - v2);
    gradient /= (2.0f * spacing_);
    return ParallelFilterValue3D(gradient);
}



ParallelFilterValue3D GradientFilter::calcGradientLinearRegression(const Sample& sample, const tgt::ivec3& pos) {

    // Euclidean weights for voxels with Manhattan distances of 1/2/3
    const float w_1 = 1.f;
    const float w_2 = 0.5f;
    const float w_3 = 1.f/3.f;

    const float w_A = 1.f / (8.f + 2.f/3.f);
    const float w_B = w_A;
    const float w_C = w_A;

    //left plane
    ParallelFilterValue1D v000 = sample(pos + tgt::ivec3(-1,-1,-1));
    ParallelFilterValue1D v001 = sample(pos + tgt::ivec3(-1, -1, 0));
    ParallelFilterValue1D v002 = sample(pos + tgt::ivec3(-1, -1, 1));
    ParallelFilterValue1D v010 = sample(pos + tgt::ivec3(-1, 0, -1));
    ParallelFilterValue1D v011 = sample(pos + tgt::ivec3(-1, 0, 0));
    ParallelFilterValue1D v012 = sample(pos + tgt::ivec3(-1, 0, 1));
    ParallelFilterValue1D v020 = sample(pos + tgt::ivec3(-1, 1, -1));
    ParallelFilterValue1D v021 = sample(pos + tgt::ivec3(-1, 1, 0));
    ParallelFilterValue1D v022 = sample(pos + tgt::ivec3(-1, 1, 1));

    //mid plane
    ParallelFilterValue1D v100 = sample(pos + tgt::ivec3(0, -1, -1));
    ParallelFilterValue1D v101 = sample(pos + tgt::ivec3(0, -1, 0));
    ParallelFilterValue1D v102 = sample(pos + tgt::ivec3(0, -1, 1));
    ParallelFilterValue1D v110 = sample(pos + tgt::ivec3(0, 0, -1));
    //ParallelFilterValue1D v111 = sample(pos + ivec3(0, 0, 0));
    ParallelFilterValue1D v112 = sample(pos + tgt::ivec3(0, 0, 1));
    ParallelFilterValue1D v120 = sample(pos + tgt::ivec3(0, 1, -1));
    ParallelFilterValue1D v121 = sample(pos + tgt::ivec3(0, 1, 0));
    ParallelFilterValue1D v122 = sample(pos + tgt::ivec3(0, 1, 1));

    //right plane
    ParallelFilterValue1D v200 = sample(pos + tgt::ivec3(1, -1, -1));
    ParallelFilterValue1D v201 = sample(pos + tgt::ivec3(1, -1, 0));
    ParallelFilterValue1D v202 = sample(pos + tgt::ivec3(1, -1, 1));
    ParallelFilterValue1D v210 = sample(pos + tgt::ivec3(1, 0, -1));
    ParallelFilterValue1D v211 = sample(pos + tgt::ivec3(1, 0, 0));
    ParallelFilterValue1D v212 = sample(pos + tgt::ivec3(1, 0, 1));
    ParallelFilterValue1D v220 = sample(pos + tgt::ivec3(1, 1, -1));
    ParallelFilterValue1D v221 = sample(pos + tgt::ivec3(1, 1, 0));
    ParallelFilterValue1D v222 = sample(pos + tgt::ivec3(1, 1, 1));

    tgt::vec3 gradient;

    gradient.x = w_1 * ( v211 - v011 )                +
                 w_2 * ( v201 + v210 + v212 + v221
                         -v001 - v010 - v012 - v021 ) +
                 w_3 * ( v200 + v202 + v220 + v222
                         -v000 - v002 - v020 - v022 );

    gradient.y = w_1 * ( v121 - v101 )                +
                 w_2 * ( v021 + v120 + v122 + v221
                         -v001 - v100 - v102 - v201 ) +
                 w_3 * ( v020 + v022 + v220 + v222
                         -v000 - v002 - v200 - v202 );

    gradient.z = w_1 * ( v112 - v110 )                +
                 w_2 * ( v012 + v102 + v122 + v212
                         -v010 - v100 - v120 - v210 ) +
                 w_3 * ( v002 + v022 + v202 + v222
                         -v000 - v020 - v200 - v220 );

    gradient.x *= w_A;
    gradient.y *= w_B;
    gradient.z *= w_C;

    gradient /= spacing_;
    gradient *= -1.f;

    return ParallelFilterValue3D(gradient);
}

ParallelFilterValue3D GradientFilter::calcGradientSobel(const Sample& sample, const tgt::ivec3& pos) {

    //left plane
    ParallelFilterValue1D v000 = sample(pos + tgt::ivec3(-1, -1, -1));
    ParallelFilterValue1D v001 = sample(pos + tgt::ivec3(-1, -1, 0));
    ParallelFilterValue1D v002 = sample(pos + tgt::ivec3(-1, -1, 1));
    ParallelFilterValue1D v010 = sample(pos + tgt::ivec3(-1, 0, -1));
    ParallelFilterValue1D v011 = sample(pos + tgt::ivec3(-1, 0, 0));
    ParallelFilterValue1D v012 = sample(pos + tgt::ivec3(-1, 0, 1));
    ParallelFilterValue1D v020 = sample(pos + tgt::ivec3(-1, 1, -1));
    ParallelFilterValue1D v021 = sample(pos + tgt::ivec3(-1, 1, 0));
    ParallelFilterValue1D v022 = sample(pos + tgt::ivec3(-1, 1, 1));
    //mid plane
    ParallelFilterValue1D v100 = sample(pos + tgt::ivec3(0, -1, -1));
    ParallelFilterValue1D v101 = sample(pos + tgt::ivec3(0, -1, 0));
    ParallelFilterValue1D v102 = sample(pos + tgt::ivec3(0, -1, 1));
    ParallelFilterValue1D v110 = sample(pos + tgt::ivec3(0, 0, -1));
    //ParallelFilterValue1D v111 = sample(pos + ivec3(0, 0, 0)); //not needed for calculation
    ParallelFilterValue1D v112 = sample(pos + tgt::ivec3(0, 0, 1));
    ParallelFilterValue1D v120 = sample(pos + tgt::ivec3(0, -1, -1));
    ParallelFilterValue1D v121 = sample(pos + tgt::ivec3(0, -1, 0));
    ParallelFilterValue1D v122 = sample(pos + tgt::ivec3(0, -1, 1));
    //right plane
    ParallelFilterValue1D v200 = sample(pos + tgt::ivec3(1, -1, -1));
    ParallelFilterValue1D v201 = sample(pos + tgt::ivec3(1, -1, 0));
    ParallelFilterValue1D v202 = sample(pos + tgt::ivec3(1, -1, 1));
    ParallelFilterValue1D v210 = sample(pos + tgt::ivec3(1, 0, -1));
    ParallelFilterValue1D v211 = sample(pos + tgt::ivec3(1, 0, 0));
    ParallelFilterValue1D v212 = sample(pos + tgt::ivec3(1, 0, 1));
    ParallelFilterValue1D v220 = sample(pos + tgt::ivec3(1, 1, -1));
    ParallelFilterValue1D v221 = sample(pos + tgt::ivec3(1, 1, 0));
    ParallelFilterValue1D v222 = sample(pos + tgt::ivec3(1, 1, 1));

    tgt::vec3 gradient = tgt::vec3::zero;

    //filter x-direction
    gradient.x += -1 * v000;
    gradient.x += -3 * v010;
    gradient.x += -1 * v020;
    gradient.x += 1 * v200;
    gradient.x += 3 * v210;
    gradient.x += 1 * v220;
    gradient.x += -3 * v001;
    gradient.x += -6 * v011;
    gradient.x += -3 * v021;
    gradient.x += +3 * v201;
    gradient.x += +6 * v211;
    gradient.x += +3 * v221;
    gradient.x += -1 * v002;
    gradient.x += -3 * v012;
    gradient.x += -1 * v022;
    gradient.x += +1 * v202;
    gradient.x += +3 * v212;
    gradient.x += +1 * v222;

    //filter y-direction
    gradient.y += -1 * v000;
    gradient.y += -3 * v100;
    gradient.y += -1 * v200;
    gradient.y += +1 * v020;
    gradient.y += +3 * v120;
    gradient.y += +1 * v220;
    gradient.y += -3 * v001;
    gradient.y += -6 * v101;
    gradient.y += -3 * v201;
    gradient.y += +3 * v021;
    gradient.y += +6 * v121;
    gradient.y += +3 * v221;
    gradient.y += -1 * v002;
    gradient.y += -3 * v102;
    gradient.y += -1 * v202;
    gradient.y += +1 * v022;
    gradient.y += +3 * v122;
    gradient.y += +1 * v222;

    //filter z-direction
    gradient.z += -1 * v000;
    gradient.z += -3 * v100;
    gradient.z += -1 * v200;
    gradient.z += +1 * v002;
    gradient.z += +3 * v102;
    gradient.z += +1 * v202;
    gradient.z += -3 * v010;
    gradient.z += -6 * v110;
    gradient.z += -3 * v210;
    gradient.z += +3 * v012;
    gradient.z += +6 * v112;
    gradient.z += +3 * v212;
    gradient.z += -1 * v020;
    gradient.z += -3 * v120;
    gradient.z += -1 * v220;
    gradient.z += +1 * v022;
    gradient.z += +3 * v122;
    gradient.z += +1 * v222;

    gradient /= 22.f;   // sum of all positive weights
    gradient /= 2.f;    // this mask has a step length of 2 voxels
    gradient /= spacing_;
    gradient *= -1.f;

    return ParallelFilterValue3D(gradient);
}

} // namespace voreen
