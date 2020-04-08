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

#include "vorticityfilter.h"

#include "slicereader.h"

namespace voreen {

VorticityFilter::VorticityFilter(GradientType gradientType, const tgt::vec3& spacing, const SamplingStrategy<ParallelFilterValue3D>& samplingStrategy, const std::string& sliceBaseType)
    : ParallelVolumeFilter<ParallelFilterValue3D, ParallelFilterValue3D>(1, samplingStrategy, sliceBaseType)
    , gradientX_(gradientType, spacing, SamplingStrategy<ParallelFilterValue1D>(samplingStrategy.type_, samplingStrategy.outsideVolumeValue_.val_[0]), sliceBaseType)
    , gradientY_(gradientType, spacing, SamplingStrategy<ParallelFilterValue1D>(samplingStrategy.type_, samplingStrategy.outsideVolumeValue_.val_[1]), sliceBaseType)
    , gradientZ_(gradientType, spacing, SamplingStrategy<ParallelFilterValue1D>(samplingStrategy.type_, samplingStrategy.outsideVolumeValue_.val_[2]), sliceBaseType)
    , gradientType_(gradientType)
{
}

ParallelFilterValue3D VorticityFilter::getValue(const Sample& sample, const tgt::ivec3& pos) const {

    // Calculate gradient for each, x, y and z.
    // TODO: this will sample the same position multiple times which could be improved
    //  by either caching or not relying on the Gradient Filter as it.
    GradientFilter::Sample sampleX = [&sample] (const tgt::ivec3& pos) {
        return ParallelFilterValue1D(sample(pos)[0]);
    };
    GradientFilter::Sample sampleY = [&sample] (const tgt::ivec3& pos) {
        return ParallelFilterValue1D(sample(pos)[1]);
    };
    GradientFilter::Sample sampleZ = [&sample] (const tgt::ivec3& pos) {
        return ParallelFilterValue1D(sample(pos)[2]);
    };

    tgt::vec3 gradientX = gradientX_.getValue(sampleX, pos).val_;
    tgt::vec3 gradientY = gradientY_.getValue(sampleY, pos).val_;
    tgt::vec3 gradientZ = gradientZ_.getValue(sampleZ, pos).val_;

    // Calculate vorticity from all gradients.
    tgt::vec3 vorticity;
    vorticity.x = gradientZ.y - gradientY.z;
    vorticity.y = gradientX.z - gradientZ.x;
    vorticity.z = gradientY.x - gradientX.y;

    return ParallelFilterValue3D(vorticity);
}

GradientType VorticityFilter::getGradientType() const {
    return gradientType_;
}

} // namespace voreen
