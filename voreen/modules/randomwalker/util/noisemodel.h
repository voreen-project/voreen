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

#ifndef VRN_RANDOM_WALKER_NOISEMODEL_H
#define VRN_RANDOM_WALKER_NOISEMODEL_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "preprocessing.h"

namespace voreen {

enum RWNoiseModel {
    RW_NOISE_GAUSSIAN,
    RW_NOISE_POISSON,
};

struct RWNoiseModelGaussian {
    static VolumeAtomic<float> preprocess(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        return preprocessForAdaptiveParameterSetting(vol);
    }
    static VolumeAtomic<float> preprocess(const VolumeRAM& vol, RealWorldMapping rwm) {
        return preprocessForAdaptiveParameterSetting(vol);
    }
    static float getEdgeWeight(float voxelIntensity, float neighborIntensity, float betaBias) {
        float beta = 0.125f * betaBias;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};
struct RWNoiseModelPoisson {
    static VolumeAtomic<float> preprocess(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        tgtAssert(vol.getNumChannels() == 1, "Only volumes with one channel expected");
        size_t voxels = vol.getNumVoxels();
        VolumeAtomic<float> converted(vol.getDimensions());
        for(size_t i=0; i<voxels; ++i) {
            converted.voxel(i) = rwm.normalizedToRealWorld(vol.voxel(i));
        }
        return converted;
    }
    static VolumeAtomic<float> preprocess(const VolumeRAM& vol, RealWorldMapping rwm) {
        tgtAssert(vol.getNumChannels() == 1, "Only volumes with one channel expected");
        size_t voxels = vol.getNumVoxels();
        VolumeAtomic<float> converted(vol.getDimensions());
        for(size_t i=0; i<voxels; ++i) {
            converted.voxel(i) = rwm.normalizedToRealWorld(vol.getVoxelNormalized(i));
        }
        return converted;
    }
    static float getEdgeWeight(float voxelIntensity, float neighborIntensity, float betaBias) {
        float beta = 0.5f * betaBias;
        float weight;

        float d = sqrt(voxelIntensity) - sqrt(neighborIntensity);
        weight = exp(-beta * d * d);
        return weight;
    }
};

}

#endif
