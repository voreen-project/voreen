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

#include "randomwalkerweights.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "tgt/vector.h"

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}

namespace voreen {

const std::string RandomWalkerWeights::loggerCat_("voreen.RandomWalker.RandomWalkerWeights");

void RandomWalkerWeights::initialize(const VolumeBase* volume, const RandomWalkerSeeds* /*seeds*/,
                                     const RandomWalkerSolver* /*solver*/) {
    volume_ = volume;
    tgtAssert(volume_, "null pointer passed");

    float min = volume_->getDerivedData<VolumeMinMax>()->getMinNormalized();
    float max = volume_->getDerivedData<VolumeMinMax>()->getMaxNormalized();
    RealWorldMapping rwm = volume_->getRealWorldMapping();
    minIntensity_ = rwm.normalizedToRealWorld(min);
    maxIntensity_ = rwm.normalizedToRealWorld(max);

    intensityScale_ = 1.f / (maxIntensity_ - minIntensity_);
    numVoxels_ = volume_->getNumVoxels();
    volDim_ = volume_->getDimensions();

    LDEBUG("intensity range: " << tgt::vec2(minIntensity_, maxIntensity_) << ", intensity scale: " << intensityScale_);
}

void RandomWalkerWeights::processVoxel(const tgt::ivec3& voxel, const RandomWalkerSeeds* seeds,
    EllpackMatrix<float>& mat, float* &vec, const RandomWalkerSolver* solver, const VolumeRAM* vol, const RealWorldMapping& rwm)
{
    tgtAssert(volume_, "no volume");
    tgtAssert(seeds, "no seed definer passed");
    tgtAssert(solver, "no solver passed");
    tgtAssert(mat.isInitialized(), "matrix not initialized");

    const int x = voxel.x;
    const int y = voxel.y;
    const int z = voxel.z;

    size_t index = volumeCoordsToIndex(voxel, volDim_);
    if (seeds->isSeedPoint(index))
        return;

    size_t curRow = solver->getRowIndex(index);

    // For performance reasons calculated outside of the function (and the loop)
    //const VolumeRAM* vol = volume_->getRepresentation<VolumeRAM>();
    //RealWorldMapping rwm = volume_->getRealWorldMapping();

    float curIntensity = rwm.normalizedToRealWorld(vol->getVoxelNormalized(voxel));

    float weightSum = 0;

    // x-neighbors
    if (x > 0) {
        tgt::ivec3 neighbor(x-1, y, z);

        size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim_);
        float neighborIntensity = rwm.normalizedToRealWorld(vol->getVoxelNormalized(neighbor));
        float weight = getEdgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = solver->getRowIndex(neighborIndex);
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
        float neighborIntensity = rwm.normalizedToRealWorld(vol->getVoxelNormalized(neighbor));
        float weight = getEdgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = solver->getRowIndex(neighborIndex);
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
        float neighborIntensity = rwm.normalizedToRealWorld(vol->getVoxelNormalized(neighbor));
        float weight = getEdgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = solver->getRowIndex(neighborIndex);
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
        float neighborIntensity = rwm.normalizedToRealWorld(vol->getVoxelNormalized(neighbor));
        float weight = getEdgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = solver->getRowIndex(neighborIndex);
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
        float neighborIntensity = rwm.normalizedToRealWorld(vol->getVoxelNormalized(neighbor));
        float weight = getEdgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = solver->getRowIndex(neighborIndex);
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
        float neighborIntensity = rwm.normalizedToRealWorld(vol->getVoxelNormalized(neighbor));
        float weight = getEdgeWeight(voxel, neighbor, curIntensity, neighborIntensity);

        if (!seeds->isSeedPoint(neighbor)) {
            size_t nRow = solver->getRowIndex(neighborIndex);
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

//---------------------------------------------------------------------------------------

RandomWalkerWeightsTransFunc::RandomWalkerWeightsTransFunc(const TransFunc1D* transFunc,
        float beta, float blendFactor, float minWeight, float maxWeight) :
   transFunc_(transFunc),
   beta_(beta),
   blendFactor_(blendFactor),
   minWeight_(minWeight),
   maxWeight_(maxWeight),
   opacityBuffer_(0),
   opacityBufferSize_(0)
{
    tgtAssert(transFunc, "null pointer passed as trans func");
    tgtAssert(beta >= 0.f, "beta must not be negative");
    tgtAssert(maxWeight_ >= 0.f, "max weight must not be negative");
    tgtAssert(minWeight <= maxWeight, "min weight must be less or equal max weight");
    tgtAssert(blendFactor_ >= 0.f && blendFactor <= 1.f, "blend factor must be between 0.0 and 1.0");
}

RandomWalkerWeightsTransFunc::~RandomWalkerWeightsTransFunc() {
    delete[] opacityBuffer_;
}

void RandomWalkerWeightsTransFunc::initialize(const VolumeBase* volume, const RandomWalkerSeeds* seeds,
        const RandomWalkerSolver* solver) {
    RandomWalkerWeights::initialize(volume, seeds, solver);

    // construct opacity buffer
    if (!transFunc_)
        throw VoreenException("No compatible transfer function. Abort.");

    opacityBufferSize_ = transFunc_->getDimensions().x;
    tgtAssert(opacityBufferSize_ > 0, "invalid transfunc dimensions");
    opacityBuffer_ = new float[opacityBufferSize_];
    if(transFunc_->getDataType() == TransFuncBase::TF_FLOAT) {
        for (size_t i=0; i<opacityBufferSize_; i++) {
            opacityBuffer_[i] = (const_cast<TransFunc1D*>(transFunc_)->getTexture()->texel<tgt::Vector4<GLfloat> >(i)).a;
        }
    } else {
        for (size_t i=0; i<opacityBufferSize_; i++) {
            opacityBuffer_[i] = static_cast<float>((const_cast<TransFunc1D*>(transFunc_)->getTexture()->texel<tgt::Vector4<GLubyte> >(i)).a) / 255.f;
        }
    }
}

float RandomWalkerWeightsTransFunc::getEdgeWeight(const tgt::ivec3& /*voxel*/, const tgt::ivec3& /*neighbor*/,
        float voxelIntensity, float neighborIntensity) const
{
    // intensity difference
    float intDiff = (voxelIntensity - neighborIntensity) * intensityScale_;
    float intDiffSqr = intDiff*intDiff;

    // Map realworld intensity values to TF space:
    voxelIntensity = transFunc_->realWorldToNormalized(voxelIntensity);
    neighborIntensity = transFunc_->realWorldToNormalized(neighborIntensity);

    // opacity difference
    int opacityIndex = tgt::ifloor(voxelIntensity * (opacityBufferSize_-1));
    int opacityIndexNeighbor = tgt::ifloor(neighborIntensity * (opacityBufferSize_-1));
    tgtAssert(opacityIndex >= 0 && opacityIndex < (int)opacityBufferSize_, "invalid opacity buffer index");
    tgtAssert(opacityIndexNeighbor >= 0 && opacityIndexNeighbor < (int)opacityBufferSize_,
            "invalid opacity buffer index");
    float opacity = opacityBuffer_[opacityIndex];
    float nOpacity = opacityBuffer_[opacityIndexNeighbor];
    float opacityDiff = nOpacity - opacity;
    float opacityDiffSqr = opacityDiff*opacityDiff;

    // blend
    float grad = (1.f - blendFactor_)*intDiffSqr + blendFactor_*opacityDiffSqr;

    // final weight
    float weight = exp(-beta_ * grad);
    weight = tgt::clamp(weight, minWeight_, maxWeight_);
    return weight;
}

//---------------------------------------------------------------------------------------

RandomWalkerWeightsIntensity::RandomWalkerWeightsIntensity(
        float beta, float minWeight, float maxWeight) :
    beta_(beta),
    minWeight_(minWeight),
    maxWeight_(maxWeight)
{
    tgtAssert(beta >= 0.f, "beta must not be negative");
    tgtAssert(maxWeight_ >= 0.f, "max weight must not be negative");
    tgtAssert(minWeight <= maxWeight, "min weight must be less or equal max weight");
}

float RandomWalkerWeightsIntensity::getEdgeWeight(const tgt::ivec3& /*voxel*/, const tgt::ivec3& /*neighbor*/,
        float voxelIntensity, float neighborIntensity) const {

    float intDiff = (voxelIntensity - neighborIntensity) * intensityScale_;
    float intDiffSqr = intDiff*intDiff;
    float weight = exp(-beta_ * intDiffSqr);
    weight = tgt::clamp(weight, minWeight_, maxWeight_);

    return weight;
}

}   // namespace
