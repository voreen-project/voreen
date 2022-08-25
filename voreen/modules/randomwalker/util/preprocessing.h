/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_RANDOM_WALKER_PREPROCESSING_H
#define VRN_RANDOM_WALKER_PREPROCESSING_H

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

VolumeAtomic<float> meanFilter3x3x3(const VolumeAtomic<float>& img);
VolumeAtomic<float> medianFilter3x3x3(const VolumeAtomic<float>& img);
float estimateVariance(const VolumeAtomic<float>& img, const VolumeAtomic<float>& mean, int filter_extent);
VolumeAtomic<float> toVolumeAtomicFloat(const VolumeRAM& img);
VolumeAtomic<float> applyRWM(const VolumeAtomic<float>& vol, RealWorldMapping rwm);
VolumeAtomic<float> applyRWM(const VolumeRAM& vol, RealWorldMapping rwm);

template<int extent>
VolumeAtomic<float> meanFilterTemplate(const VolumeAtomic<float>& img);
VolumeAtomic<float> meanFilter(const VolumeAtomic<float>& img, int extent);
VolumeAtomic<float> variances(const VolumeAtomic<float>& img, const VolumeAtomic<float>& mean, int extent);

/// ----------------------------------------------------------------------------
/// Implementation -------------------------------------------------------------
/// ----------------------------------------------------------------------------

template<int extent>
VolumeAtomic<float> meanFilterTemplate(const VolumeAtomic<float>& img) {
    const tgt::ivec3 start(0);
    const tgt::ivec3 end(img.getDimensions());
    const size_t numVoxels = tgt::hmul(img.getDimensions());

    const int k = extent;
    const int N=2*k+1;
    const tgt::ivec3 neighborhoodSize(k);

    // mean
    auto conv = [&] (const VolumeAtomic<float>& input, VolumeAtomic<float>& output, int dim) {
        VRN_FOR_EACH_VOXEL(center, start, end) {
            tgt::ivec3 neigh(0);
            neigh[dim] = neighborhoodSize[dim];
            const tgt::ivec3 neighborhoodStart = tgt::max(start, center - neigh);
            const tgt::ivec3 neighborhoodEnd = tgt::min(end, center + neigh + tgt::ivec3(1));

            const int numNeighborhoodVoxels = tgt::hmul(neighborhoodEnd-neighborhoodStart);

            float sum=0.0f;
            VRN_FOR_EACH_VOXEL(pos, neighborhoodStart, neighborhoodEnd) {
                sum += input.voxel(pos);
            }
            float estimation = sum/numNeighborhoodVoxels;
            output.voxel(center) = estimation;
        }
    };
    VolumeAtomic<float> tmp(img.getDimensions());
    VolumeAtomic<float> tmp2(img.getDimensions());
    conv(img, tmp2, 0);
    conv(tmp2, tmp, 1);
    conv(tmp, tmp2, 2);

    return tmp2;
}

}

#endif
