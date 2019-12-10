/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "randomwalkersolver.h"

#include "randomwalkerseeds.h"
#include "randomwalkerweights.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "tgt/vector.h"

namespace voreen {

const std::string RandomWalkerSolver::loggerCat_("voreen.RandomWalker.RandomWalkerSolver");

RandomWalkerSolver::RandomWalkerSolver(const VolumeBase* volume,
        RandomWalkerSeeds* seeds, RandomWalkerWeights* edgeWeights) :
    volume_(volume),
    seeds_(seeds),
    edgeWeights_(edgeWeights),
    mat_(),
    vec_(0),
    volIndexToRow_(0),
    solution_(0),
    numSeeds_(0),
    state_(Initial)
{
    tgtAssert(volume_, "no volume");
    tgtAssert(seeds_, "no seed point definer");
    tgtAssert(edgeWeights_, "no edge weight calculator");

    volDim_ = volume->getDimensions();
    numVoxels_ = tgt::hmul(volDim_);
    volSpacing_ = volume->getSpacing();
    volTransformation_ = volume->getPhysicalToWorldMatrix();
}

RandomWalkerSolver::~RandomWalkerSolver() {
    delete[] vec_;
    delete[] volIndexToRow_;
    delete[] solution_;
    vec_ = 0;
    volIndexToRow_ = 0;
    solution_ = 0;

    delete seeds_;
    delete edgeWeights_;
    seeds_ = 0;
    edgeWeights_ = 0;
}

void RandomWalkerSolver::setupEquationSystem() {

    if (state_ != Initial)
        throw VoreenException("System has already been setup");

    tgtAssert(!mat_.isInitialized(), "matrix already created");
    tgtAssert(!vec_, "vector already created");
    tgtAssert(!volIndexToRow_, "volIndexToRow buffer vector already created");

    state_ = Failure;

    // initialize seeds
    try {
        seeds_->initialize(volume_->getRepresentation<VolumeRAM>());
    }
    catch (VoreenException& e) {
        throw VoreenException("Failed to initialize seeds: " + std::string(e.what()));
    }
    numSeeds_ = seeds_->getNumSeeds();
    if (numSeeds_ == 0)
        throw VoreenException("No seed points");

    size_t systemSize = numVoxels_ - numSeeds_;

    // initialize matrix
    mat_.setDimensions(systemSize, systemSize, 7);
    try {
        mat_.initializeBuffers();
    }
    catch (VoreenException& e) {
        throw VoreenException("Failed to set up matrix: " + std::string(e.what()));
    }

    // initialize vector
    try {
        vec_ = new float[systemSize];
    }
    catch (std::bad_alloc&) {
        throw VoreenException("Bad allocation during creation of vec buffer");
    }
    for (size_t i=0; i<systemSize; i++)
        vec_[i] = 0.f;

    // compute mapping from voxel index to row index
    // (leaving out seed voxels)
    computeVolIndexToRowMapping(seeds_);
    tgtAssert(volIndexToRow_, "volIndexToRowBuffer empty");

    // initialize edge weight computer
    try {
        edgeWeights_->initialize(volume_, seeds_, this);
    }
    catch (VoreenException& e) {
        throw VoreenException("Failed to initialize edge weights: " + std::string(e.what()));
    }

    // For performance reasons, get VolumeRAM and rwm now and pass it to processVoxel.
    const VolumeRAM* volram = volume_->getRepresentation<VolumeRAM>();
    const RealWorldMapping& rwm = volume_->getRealWorldMapping();
    // iterate over volume and compute edge weights for each voxel
    #ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
    #endif
    for (int z=0; z<volDim_.z; z++) {
        for (int y=0; y<volDim_.y; y++) {
            for (int x=0; x<volDim_.x; x++) {
                edgeWeights_->processVoxel(tgt::ivec3(x, y, z), seeds_, mat_, vec_, this, volram, rwm);
            }
        }
    }

    state_ = Setup;
}

int RandomWalkerSolver::solve(const VoreenBlas* voreenBlas, float* oldSystemSolution, VoreenBlas::ConjGradPreconditioner preConditioner, float errorThreshold, int maxIterations, ProgressReporter& progress) {

    tgtAssert(voreenBlas, "null pointer passed");

    if (state_ != Setup)
        throw VoreenException("System is not setup or has already been solved");
    tgtAssert(mat_.isInitialized(), "matrix not initialized");
    tgtAssert(vec_, "vector not created");
    tgtAssert(volIndexToRow_, "volIndexToRow buffer vector not created");
    tgtAssert(!solution_, "solution buffer already created");

    // create solution buffer
    try {
        solution_ = new float[getSystemSize()];
    }
    catch (std::bad_alloc&) {
        throw VoreenException("Bad allocation during creation of solution buffer");
    }

    int iterations = voreenBlas->sSpConjGradEll(mat_, vec_, solution_, oldSystemSolution,
        preConditioner, errorThreshold, maxIterations, progress);
    state_ = Solved;
    return iterations;
}

EllpackMatrix<float>& RandomWalkerSolver::getMatrix() {
    return mat_;
}

const EllpackMatrix<float>& RandomWalkerSolver::getMatrix() const {
    return mat_;
}

const float* RandomWalkerSolver::getVector() const {
    return vec_;
}

float* RandomWalkerSolver::getVector() {
    return vec_;
}

const float* RandomWalkerSolver::getSolution() const {
    return solution_;
}

float* RandomWalkerSolver::getSolution() {
    return solution_;
}

float RandomWalkerSolver::getProbabilityValue(size_t voxel) const {
    tgtAssert(state_ == Solved, "system has not been solved");

    if (seeds_->isSeedPoint(voxel))
        return seeds_->getSeedValue(voxel);
    else
        return solution_[volIndexToRow_[voxel]];
}

tgt::vec2 RandomWalkerSolver::getProbabilityRange() const {
    tgtAssert(state_ == Solved, "system has not been solved");

    float minProb = getProbabilityValue(0);
    float maxProb = getProbabilityValue(0);
    for (size_t voxel=0; voxel<numVoxels_; voxel++) {
        const float val = getProbabilityValue(voxel);
        minProb = std::min(minProb, val);
        maxProb = std::max(maxProb, val);
    }

    return tgt::vec2(minProb, maxProb);
}

bool RandomWalkerSolver::isSeedPoint(size_t voxel) const {
    return seeds_->isSeedPoint(voxel);
}

float RandomWalkerSolver::getSeedValue(size_t voxel) const {
    return seeds_->getSeedValue(voxel);
}

size_t RandomWalkerSolver::getNumSeeds() const {
    return numSeeds_;
}

size_t RandomWalkerSolver::getSystemSize() const {
    return numVoxels_ - numSeeds_;
}

tgt::ivec3 RandomWalkerSolver::getVolumeDimensions() const {
    return volDim_;
}

size_t RandomWalkerSolver::getNumVoxels() const {
    return numVoxels_;
}

size_t RandomWalkerSolver::getRowIndex(size_t voxel) const {
    return volIndexToRow_[voxel];
}

void RandomWalkerSolver::computeVolIndexToRowMapping(const RandomWalkerSeeds* seeds) {

    try {
        volIndexToRow_ = new size_t[numVoxels_];
    }
    catch (std::bad_alloc&) {
        throw VoreenException("Bad allocation during creation of volIndexToRow buffer");
    }

    // compute volIndexToRow values
    size_t curRow = 0;
    for (size_t i=0; i<numVoxels_; i++) {
        if (!seeds->isSeedPoint(i)) {
            volIndexToRow_[i] = curRow;
            curRow++;
        }
        else {
            volIndexToRow_[i] = -1;
        }
    }
}

RandomWalkerSolver::SystemState RandomWalkerSolver::getSystemState() const {
    return state_;
}

tgt::vec2 RandomWalkerSolver::getSeedRange() const {
    tgtAssert(state_ >= Setup, "system has not been setup");
    return seeds_->getSeedRange();
}


}   // namespace
