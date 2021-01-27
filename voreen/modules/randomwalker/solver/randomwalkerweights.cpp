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

#include "randomwalkerweights.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "tgt/vector.h"
#include "tgt/memory.h"

namespace {

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}

}

namespace voreen {

const std::string RandomWalkerWeights::loggerCat_("voreen.RandomWalker.RandomWalkerWeights");

RandomWalkerWeights::RandomWalkerWeights(std::unique_ptr<RandomWalkerVoxelAccessor> voxelFun, std::unique_ptr<RandomWalkerEdgeWeight> weightFun, tgt::ivec3 volDim)
    : voxelFun_(std::move(voxelFun))
    , weightFun_(std::move(weightFun))
    , volDim_(volDim)
{
}

void RandomWalkerWeights::processVoxel(const tgt::ivec3& voxel, const RandomWalkerSeeds* seeds, EllpackMatrix<float>& mat, float* vec, size_t* volumeIndexToRowTable)
{
    tgtAssert(seeds, "no seed definer passed");
    tgtAssert(volumeIndexToRowTable, "no volumeIndexToRowTable passed");
    tgtAssert(mat.isInitialized(), "matrix not initialized");

    const int x = voxel.x;
    const int y = voxel.y;
    const int z = voxel.z;

    size_t index = volumeCoordsToIndex(voxel, volDim_);
    if (seeds->isSeedPoint(index))
        return;

    size_t curRow = volumeIndexToRowTable[index];

    float curIntensity = voxelFun_->voxel(voxel);

    float weightSum = 0;

    // x-neighbors
    if (x > 0) {
        tgt::ivec3 neighbor(x-1, y, z);

        size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim_);
        float neighborIntensity = voxelFun_->voxel(neighbor);
        float weight = weightFun_->edgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = volumeIndexToRowTable[neighborIndex];
            //tgtAssert(nRow >= 0 && nRow < numUnseeded_, "Invalid row");
            mat.setValue(curRow, nRow, -weight);
        }
        else {
            vec[curRow] += weight * seeds->getSeedValue(neighbor);
        }

        weightSum += weight;
    }
    if (x < volDim_.x-1) {
        tgt::ivec3 neighbor(x+1, y, z);

        size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim_);
        float neighborIntensity = voxelFun_->voxel(neighbor);
        float weight = weightFun_->edgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = volumeIndexToRowTable[neighborIndex];
            //tgtAssert(nRow >= 0 && nRow < numUnseeded_, "Invalid row");
            mat.setValue(curRow, nRow, -weight);
        }
        else {
            vec[curRow] += weight * seeds->getSeedValue(neighbor);
        }

        weightSum += weight;
    }

    // y-neighbors
    if (y > 0) {
        tgt::ivec3 neighbor(x, y-1, z);

        size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim_);
        float neighborIntensity = voxelFun_->voxel(neighbor);
        float weight = weightFun_->edgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = volumeIndexToRowTable[neighborIndex];
            //tgtAssert(nRow >= 0 && nRow < numUnseeded_, "Invalid row");
            mat.setValue(curRow, nRow, -weight);
        }
        else {
            vec[curRow] += weight * seeds->getSeedValue(neighbor);
        }

        weightSum += weight;
    }
    if (y < volDim_.y-1) {
        tgt::ivec3 neighbor(x, y+1, z);

        size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim_);
        float neighborIntensity = voxelFun_->voxel(neighbor);
        float weight = weightFun_->edgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = volumeIndexToRowTable[neighborIndex];
            //tgtAssert(nRow >= 0 && nRow < numUnseeded_, "Invalid row");
            mat.setValue(curRow, nRow, -weight);
        }
        else {
            vec[curRow] += weight * seeds->getSeedValue(neighbor);
        }

        weightSum += weight;
    }

    // z-neighbors
    if (z > 0) {
        tgt::ivec3 neighbor(x, y, z-1);

        size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim_);
        float neighborIntensity = voxelFun_->voxel(neighbor);
        float weight = weightFun_->edgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = volumeIndexToRowTable[neighborIndex];
            //tgtAssert(nRow >= 0 && nRow < numUnseeded_, "Invalid row");
            mat.setValue(curRow, nRow, -weight);
        }
        else {
            vec[curRow] += weight * seeds->getSeedValue(neighbor);
        }

        weightSum += weight;
    }
    if (z < volDim_.z-1) {
        tgt::ivec3 neighbor(x, y, z+1);

        size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim_);
        float neighborIntensity = voxelFun_->voxel(neighbor);
        float weight = weightFun_->edgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = volumeIndexToRowTable[neighborIndex];
            //tgtAssert(nRow >= 0 && nRow < numUnseeded_, "Invalid row");
            mat.setValue(curRow, nRow, -weight);
        }
        else {
            vec[curRow] += weight * seeds->getSeedValue(neighbor);
        }

        weightSum += weight;
    }

    mat.setValue(curRow, curRow, weightSum);
}

RandomWalkerVoxelAccessorVolume::RandomWalkerVoxelAccessorVolume(const VolumeBase& volume)
    : vol_(volume.getRepresentation<VolumeRAM>())
    , rwm_(volume.getRealWorldMapping())
{
}

float RandomWalkerVoxelAccessorVolume::voxel(const tgt::svec3& voxel) {
    return rwm_.normalizedToRealWorld(vol_->getVoxelNormalized(voxel));
}

RandomWalkerVoxelAccessorVolumeAtomic::RandomWalkerVoxelAccessorVolumeAtomic(VolumeAtomic<float>&& volume, RealWorldMapping rwm)
    : vol_(std::move(volume))
    , rwm_(rwm)
{
}

float RandomWalkerVoxelAccessorVolumeAtomic::voxel(const tgt::svec3& voxel) {
    return rwm_.normalizedToRealWorld(vol_.voxel(voxel));
}

//---------------------------------------------------------------------------------------
RandomWalkerEdgeWeightTransfunc::RandomWalkerEdgeWeightTransfunc(const TransFunc1D* transFunc, tgt::vec2 intensityRange, float beta, float blendFactor, float minWeight, float maxWeight)
    : transFunc(transFunc)
    , beta(beta)
    , blendFactor(blendFactor)
    , minWeight(minWeight)
    , maxWeight(maxWeight)
    , opacityBuffer(transFunc->getDimensions().x, 0.0f)
    , intensityScale(1.f / (intensityRange.y - intensityRange.x))
{
    tgtAssert(transFunc, "null pointer passed as trans func");
    tgtAssert(beta >= 0.f, "beta must not be negative");
    tgtAssert(maxWeight >= 0.f, "max weight must not be negative");
    tgtAssert(minWeight <= maxWeight, "min weight must be less or equal max weight");
    tgtAssert(blendFactor >= 0.f && blendFactor <= 1.f, "blend factor must be between 0.0 and 1.0");

    // construct opacity buffer
    if (!transFunc) {
        throw VoreenException("No compatible transfer function. Abort.");
    }


    if(transFunc->getDataType() == TransFuncBase::TF_FLOAT) {
        for (size_t i=0; i<opacityBuffer.size(); i++) {
            opacityBuffer[i] = (const_cast<TransFunc1D*>(transFunc)->getTexture()->texel<tgt::Vector4<GLfloat> >(i)).a;
        }
    } else {
        for (size_t i=0; i<opacityBuffer.size(); i++) {
            opacityBuffer[i] = static_cast<float>((const_cast<TransFunc1D*>(transFunc)->getTexture()->texel<tgt::Vector4<GLubyte> >(i)).a) / 255.f;
        }
    }
}

float RandomWalkerEdgeWeightTransfunc::edgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor, float voxelIntensity, float neighborIntensity) {
    // intensity difference
    float intDiff = (voxelIntensity - neighborIntensity) * intensityScale;
    float intDiffSqr = intDiff*intDiff;

    // Map realworld intensity values to TF space:
    voxelIntensity = transFunc->realWorldToNormalized(voxelIntensity);
    neighborIntensity = transFunc->realWorldToNormalized(neighborIntensity);

    // opacity difference
    int opacityIndex = tgt::ifloor(voxelIntensity * (opacityBuffer.size()-1));
    int opacityIndexNeighbor = tgt::ifloor(neighborIntensity * (opacityBuffer.size()-1));
    tgtAssert(opacityIndex >= 0 && opacityIndex < (int)opacityBuffer.size(), "invalid opacity buffer index");
    tgtAssert(opacityIndexNeighbor >= 0 && opacityIndexNeighbor < (int)opacityBuffer.size(),
            "invalid opacity buffer index");
    float opacity = opacityBuffer[opacityIndex];
    float nOpacity = opacityBuffer[opacityIndexNeighbor];
    float opacityDiff = nOpacity - opacity;
    float opacityDiffSqr = opacityDiff*opacityDiff;

    // blend
    float grad = (1.f - blendFactor)*intDiffSqr + blendFactor*opacityDiffSqr;

    // final weight
    float weight = exp(-beta * grad);
    weight = tgt::clamp(weight, minWeight, maxWeight);
    return weight;
}

//---------------------------------------------------------------------------------------

RandomWalkerEdgeWeightIntensity::RandomWalkerEdgeWeightIntensity(tgt::vec2 intensityRange, float beta, float minWeight, float maxWeight)
    : beta(beta)
    , minWeight(minWeight)
    , maxWeight(maxWeight)
    , intensityScale(1.f / (intensityRange.y - intensityRange.x))
{
    tgtAssert(beta >= 0.f, "beta must not be negative");
    tgtAssert(maxWeight >= 0.f, "max weight must not be negative");
    tgtAssert(minWeight <= maxWeight, "min weight must be less or equal max weight");
}


float RandomWalkerEdgeWeightIntensity::edgeWeight(const tgt::ivec3& voxel, const tgt::ivec3& neighbor, float voxelIntensity, float neighborIntensity) {
    float intDiff = (voxelIntensity - neighborIntensity) * intensityScale;
    float intDiffSqr = intDiff*intDiff;
    float weight = exp(-beta * intDiffSqr);
    weight = tgt::clamp(weight, minWeight, maxWeight);

    return weight;
}

}   // namespace
