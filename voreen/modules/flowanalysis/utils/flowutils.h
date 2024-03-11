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

#ifndef VRN_FLOWUTILS_H
#define VRN_FLOWUTILS_H

#include "voreen/core/datastructures/volume/volumeram.h"

namespace voreen {

class VelocitySampler {
public:

    virtual ~VelocitySampler() {};
    virtual tgt::vec3 sample(tgt::vec3 pos) const = 0;
};

class SpatialSampler : public VelocitySampler {
public:
    SpatialSampler(const VolumeRAM* volume,
                   const RealWorldMapping& rwm,
                   VolumeRAM::Filter filter,
                   const tgt::mat4& toVoxelMatrix = tgt::mat4::identity,
                   const tgt::mat4& velocityTransformationMatrix = tgt::mat4::identity);

    /**
     * Samples the given volume at the given specified position.
     * @param pos Position in world space.
     */
    tgt::vec3 sample(tgt::vec3 pos) const;

private:
    const tgt::mat4 toVoxelMatrix_;
    const bool toVoxelMatrixSet_;
    const tgt::mat4 velocityTransformationMatrix_;
    const bool velocityTransformationMatrixSet_;
    std::function<tgt::vec3(tgt::vec3)> sampleFunction_;
};

class SpatioTemporalSampler : public VelocitySampler {
public:

    /**
     * Initializes a spatio-temporal sampler that samples linearly between two volumes.
     * @param volume0 first volume
     * @param volume1 second volume
     * @param alpha value in range [0,1]. alpha = 0 will return sample in volume0, alpha = 1 in volume1, respectively
     * @param rwm
     * @param filter spacial sampling strategy
     * @param toVoxelMatrix
     */
    SpatioTemporalSampler(const VolumeRAM* volume0, const VolumeRAM* volume1,
                          float alpha,
                          const RealWorldMapping& rwm,
                          VolumeRAM::Filter filter,
                          const tgt::mat4& toVoxelMatrix = tgt::mat4::identity,
                          const tgt::mat4& velocityTransformationMatrix = tgt::mat4::identity);

    /**
     * Samples the given volume at the given specified position.
     * @param pos Position in world space.
     */
    tgt::vec3 sample(tgt::vec3 pos) const;

private:
    const SpatialSampler filter0_;
    const SpatialSampler filter1_;
    const float alpha_;
};

tgt::mat4 createTransformationMatrix(const tgt::vec3& position, const tgt::vec3& velocity);

}


#endif