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

#ifndef VRN_RANDOMWALKEREDGEWEIGHTS_H
#define VRN_RANDOMWALKEREDGEWEIGHTS_H

#include "randomwalkersolver.h"
#include "randomwalkerseeds.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include <string>

namespace voreen {

//---------------------------------------------------------------------------------------

struct RandomWalkerVoxelAccessor {
    virtual ~RandomWalkerVoxelAccessor() {}
    virtual float voxel(const tgt::svec3& voxel) = 0;
};

struct RandomWalkerVoxelAccessorVolume : public RandomWalkerVoxelAccessor {
    RandomWalkerVoxelAccessorVolume(const VolumeBase& vol);
    virtual float voxel(const tgt::svec3& voxel);
private:
    const VolumeRAM* vol_;
    RealWorldMapping rwm_;
};
struct RandomWalkerVoxelAccessorVolumeAtomic : public RandomWalkerVoxelAccessor {
    RandomWalkerVoxelAccessorVolumeAtomic(VolumeAtomic<float>&& vol, RealWorldMapping rwm);
    virtual float voxel(const tgt::svec3& voxel);
private:
    VolumeAtomic<float> vol_;
    RealWorldMapping rwm_;
};

//---------------------------------------------------------------------------------------

struct RandomWalkerEdgeWeight {
    virtual ~RandomWalkerEdgeWeight() {}
    virtual float edgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor) = 0;
};

/**
 * Derives edge weights from the intensity difference of two neighored voxels:
 *     w_{ij} = e^{-\beta (int_i - int_j)^2}
 *
 * The weights are clamped to the range [minWeight, maxWeight].
 */
struct RandomWalkerEdgeWeightIntensity : public RandomWalkerEdgeWeight {
    RandomWalkerEdgeWeightIntensity(VolumeAtomic<float> vol, tgt::vec2 intensityRange, float beta = 4000.f, float minWeight = 1e-6f, float maxWeight = 1.f);
    virtual float edgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor);
private:
    VolumeAtomic<float> vol;
    float beta;
    float minWeight;
    float maxWeight;
    float intensityScale;
};

template<typename NoiseModel>
struct RandomWalkerEdgeWeightAdaptive: public RandomWalkerEdgeWeight {
    RandomWalkerEdgeWeightAdaptive(NoiseModel model, float minWeight)
        : model_(std::move(model))
        , minWeight_(minWeight)
    {
    }
    float edgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor) {
        float weight = model_.getEdgeWeight(voxel, neighbor, 1.0f);
        return std::max(minWeight_, weight);
    }
private:
    NoiseModel model_;
    float minWeight_;
};

//---------------------------------------------------------------------------------------

/**
 * Derives edge weights by blending the intensity difference between neighored voxels
 * with their alpha-value difference obtained by applying a transfer function:
 *     w_{ij} = e^{-\beta ( (1.0-blendFactor)*(int_i - int_j)^2 + blendFactor*(tf(int_i) - tf(int_j))^2 )}
 *
 * The weights are clamped to the range [minWeight, maxWeight].
 */
struct RandomWalkerEdgeWeightTransfunc : public RandomWalkerEdgeWeight {
    RandomWalkerEdgeWeightTransfunc(VolumeAtomic<float> vol, const TransFunc1D* transFunc, tgt::vec2 intensityRange, float beta = 4000.f, float blendFactor = 0.5f, float minWeight = 1e-6f, float maxWeight = 1.f);
    virtual float edgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor);
private:
    VolumeAtomic<float> vol;
    const TransFunc1D* transFunc;
    float beta;
    float blendFactor;
    float minWeight;
    float maxWeight;
    std::vector<float> opacityBuffer;
    float intensityScale;
};

/**
 * Composite class for Random Walker edge weight computation,
 * used by RandomWalkerSolver.
 */
class RandomWalkerWeights {

public:
    RandomWalkerWeights(std::unique_ptr<RandomWalkerEdgeWeight> weightFun, tgt::ivec3 volDim);
    virtual ~RandomWalkerWeights() {}

    virtual void processVoxel(const tgt::ivec3& voxel, const RandomWalkerSeeds* seeds, EllpackMatrix<float>& mat, float* vec, const size_t* volumeIndexToRowTable, boost::mutex& mat_mutex, boost::mutex& vec_mutex);

protected:
    std::unique_ptr<RandomWalkerEdgeWeight> weightFun_;
    tgt::ivec3 volDim_;

    static const std::string loggerCat_;
};

} //namespace

#endif
