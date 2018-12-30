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

#ifndef VRN_RANDOMWALKEREDGEWEIGHTS_H
#define VRN_RANDOMWALKEREDGEWEIGHTS_H

#include "randomwalkersolver.h"
#include "randomwalkerseeds.h"

#include <string>

namespace voreen {

class TransFunc1D;

/**
 * Abstract base class for Random Walker edge weight computation,
 * used by RandomWalkerSolver.
 */
class RandomWalkerWeights {

public:
    virtual ~RandomWalkerWeights() {}

    virtual void initialize(const VolumeBase* volume, const RandomWalkerSeeds* seeds,
        const RandomWalkerSolver* solver);

    virtual void processVoxel(const tgt::ivec3& voxel, const RandomWalkerSeeds* seeds,
        EllpackMatrix<float>& mat, float* &vec, const RandomWalkerSolver* solver,
        const VolumeRAM* vol, const RealWorldMapping& rwm);

    virtual float getEdgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor,
        float voxelIntensity, float neighborIntensity) const = 0;

protected:
    const VolumeBase* volume_;
    size_t numVoxels_;
    tgt::ivec3 volDim_;
    float minIntensity_;
    float maxIntensity_;
    float intensityScale_;

    static const std::string loggerCat_;
};

//---------------------------------------------------------------------------------------

/**
 * Derives edge weights from the intensity difference of two neighored voxels:
 *     w_{ij} = e^{-\beta (int_i - int_j)^2}
 *
 * The weights are clamped to the range [minWeight, maxWeight].
 */
class RandomWalkerWeightsIntensity : public RandomWalkerWeights {

public:
    RandomWalkerWeightsIntensity(float beta = 4000.f,
        float minWeight = 1e-6f, float maxWeight = 1.f);

    virtual float getEdgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor,
        float voxelIntensity, float neighborIntensity) const;

protected:
    const float beta_;
    const float minWeight_;
    const float maxWeight_;
};

//---------------------------------------------------------------------------------------

/**
 * Derives edge weights by blending the intensity difference between neighored voxels
 * with their alpha-value difference obtained by applying a transfer function:
 *     w_{ij} = e^{-\beta ( (1.0-blendFactor)*(int_i - int_j)^2 + blendFactor*(tf(int_i) - tf(int_j))^2 )}
 *
 * The weights are clamped to the range [minWeight, maxWeight].
 */
class RandomWalkerWeightsTransFunc : public RandomWalkerWeights {

public:
    RandomWalkerWeightsTransFunc(const TransFunc1D* transFunc, float beta = 4000.f,
        float blendFactor = 0.5f, float minWeight = 1e-6f, float maxWeight = 1.f);

    ~RandomWalkerWeightsTransFunc();

    virtual void initialize(const VolumeBase* volume, const RandomWalkerSeeds* seeds,
        const RandomWalkerSolver* solver);

    virtual float getEdgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor,
        float voxelIntensity, float neighborIntensity) const;

private:
    const TransFunc1D* transFunc_;
    const float beta_;
    const float blendFactor_;
    const float minWeight_;
    const float maxWeight_;

    float* opacityBuffer_;
    size_t opacityBufferSize_;
};

} //namespace

#endif
