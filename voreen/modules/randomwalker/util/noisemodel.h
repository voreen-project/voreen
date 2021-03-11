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

#ifndef VRN_RANDOM_WALKER_NOISEMODEL_H
#define VRN_RANDOM_WALKER_NOISEMODEL_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "preprocessing.h"

namespace voreen {

enum RWNoiseModel {
    RW_NOISE_GAUSSIAN,
    RW_NOISE_POISSON,
    RW_NOISE_GAUSSIAN_BIAN,
};

// Parameter estimation according to
//
// A. Bian, X Jiang: Statistical Modeling Based Adaptive Parameter Setting for Random Walk Segmentation
// https://link.springer.com/chapter/10.1007%2F978-3-319-48680-2_61
struct RWNoiseModelGaussianBian {
    VolumeAtomic<float> mean;
    float diff_variance_inv;

    RWNoiseModelGaussianBian(RWNoiseModelGaussianBian&&) = default;

    static RWNoiseModelGaussianBian prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        // Careful: This is _only_ valid for the mean filter as an estimator!!
        const int k = 1;
        const int N=2*k+1;
        const int N4=N*N*N*N;
        float correction_factor = 2.0f/N4;

        float diff_variance = variance * correction_factor;
        diff_variance = std::max(diff_variance, std::numeric_limits<float>::min());

        return RWNoiseModelGaussianBian {
            std::move(mean),
            1.0f/diff_variance,
        };
    }
    static RWNoiseModelGaussianBian prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelGaussianBian::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        float voxelIntensity = mean.voxel(voxel);
        float neighborIntensity = mean.voxel(neighbor);
        float beta = 2.0f * betaBias * diff_variance_inv;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};

struct RWNoiseModelGaussian {
    VolumeAtomic<float> values;
    float variance_inv;

    RWNoiseModelGaussian(RWNoiseModelGaussian&&) = default;

    static RWNoiseModelGaussian prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        return RWNoiseModelGaussian {
            vol.copy(),
            1.0f/variance,
        };
    }
    static RWNoiseModelGaussian prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelGaussian::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        float voxelIntensity = values.voxel(voxel);
        float neighborIntensity = values.voxel(neighbor);

        float beta = 0.125f * betaBias * variance_inv;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};
struct RWNoiseModelPoisson {
    VolumeAtomic<float> values;

    RWNoiseModelPoisson(RWNoiseModelPoisson&&) = default;

    static RWNoiseModelPoisson prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        return RWNoiseModelPoisson {
            applyRWM(vol, rwm),
        };
    }
    static RWNoiseModelPoisson prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelPoisson {
            applyRWM(vol, rwm),
        };
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        float voxelIntensity = values.voxel(voxel);
        float neighborIntensity = values.voxel(neighbor);

        float beta = 0.5f * betaBias;
        float weight;

        float d = sqrt(voxelIntensity) - sqrt(neighborIntensity);
        weight = exp(-beta * d * d);
        return weight;
    }
};

}

#endif
