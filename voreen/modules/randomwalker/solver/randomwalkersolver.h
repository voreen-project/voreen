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

#ifndef VRN_RANDOMWALKERSOLVER_H
#define VRN_RANDOMWALKERSOLVER_H

#include "voreen/core/utils/voreenblas/voreenblas.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/utils/exception.h"
#include "tgt/vector.h"

#include <string>
#include <ostream>

namespace voreen {

class Volume;
class RandomWalkerSeeds;
class RandomWalkerWeights;

/**
 * Framework for computing a 3D random walker solution on a voxel grid
 * using a conjugate gradient solver.
 *
 * For setting up the random walker equation system, the solver makes use of helper objects
 * that define the seeds and the edge weights. You may either implement the respective
 * interfaces yourself, or use the basic implementations already provided.
 *
 * For solving the equation system the solver employs a VoreenBlas helper object
 * that provides level 1 and 2 BLAS operations as well as a conjugate gradient solver.
 * Besides a basic CPU implementation, OpenMP and OpenCL versions are available.
 * The latter two are part of the 'openmp' and 'opencl' modules, respectively.
 * Note: Due to the heavy computational load involved in solving the random walker problem
 * on a 3D grid, we strongly recommend to use the OpenCL implementation!
 *
 * The solver undergoes the life cycle initialization -> problem setup -> solved. Some operations
 * are only permitted in certain stages. If an error occurs in a transition between two stages,
 * the solver enters the failure state.
 *
 * Usage:
 * 1. construct the solver by passing the volume along with the seed and edge weight definers
 * 2. call setupEquationSystem() on the solver
 * 3. if setup successful, call solve() and pass the VoreenBlas instance to use
 * 4. if computation successful, either retrieve the random walker solution as binary buffer
 *    or call generateBinarySegmentation() / generateProbabilityVolume() to retrieve the segmented volume
 *
 */
class RandomWalkerSolver {

public:
    enum SystemState {
        Failure,
        Initial,
        Setup,
        Solved
    };

    /**
     * @param volume the volume on whose voxel grid the random walker solution is to be computed
     * @param seeds helper object that specifies which voxels of the volume are seeds
     * @param edgeWeights helper object that defines weight of an edge between two voxels
     *
     * @see RandomWalkerTwoLabelSeeds, RandomWalkerWeightsIntensity, RandomWalkerWeightsTransFunc
     */
    RandomWalkerSolver(const VolumeBase* volume, RandomWalkerSeeds* seeds, RandomWalkerWeights* edgeWeights);

    /// Frees the allocated resources.
    ~RandomWalkerSolver();

    /**
     * Constructs the random walker equation system. On success, the solver enters the state 'Setup'.
     *
     * @throw VoreenException if setup has failed
     */
    void setupEquationSystem();

    /**
     * Computes the random walker solution using a conjugate gradient solver provided by the passed
     * VoreenBlas helper object. This operation is only allowed, if the equation system has been setup (state 'Setup').
     *
     * The computation finishes when either the equation error falls below the specified threshold or the
     * maximum number of conjugate gradient iterations is reached. The solver then enters the state 'Solved'.
     *
     * @see VoreenBlasCPU, VoreenBlasMP, VoreenBlasCL
     *
     * @param voreenBlas VoreenBlas implementation to use for solving the equation system
     * @param preConditioner Preconditioner to be used by the conjugate gradient solver
     * @param errorThreshold error value at which the equation system is considered to be solved
     * @param maxIterations maximum number of iterations, pass 0 for no limit
     *
     * @return the number of iterations required to compute the solution
     * @throw VoreenException if the computation has failed
     */
    int solve(const VoreenBlas* voreenBlas, float* oldSystemSolution, VoreenBlas::ConjGradPreconditioner preConditioner, float errorThreshold, int maxIterations, ProgressReporter& progress);

    /// Returns the current state of the solver.
    SystemState getSystemState() const;

    /**
     * Returns the quadratic matrix representing the random walker equation system.
     * The matrix contains one row/col per unseeded voxel.
     * Only allowed in state >= 'Setup'.
     */
    const EllpackMatrix<float>& getMatrix() const;
    EllpackMatrix<float>& getMatrix();

    /**
     * Returns the vector of the random walker equation system.
     * The number of vector elements equals the number of unseeded voxels.
     * Only allowed in state >= 'Setup'.
     */
    const float* getVector() const;
    float* getVector();

    /**
     * Returns the random walker solution consisting of a vector that contains
     * one probability value for each unseeded voxel.
     * Only allowed in state 'Solved'.
     */
    const float* getSolution() const;
    float* getSolution();

    /**
     * Returns the computed probability value of the voxel with the passed index.
     * The probability value of seed voxels equals the seed value.
     * Operation only allowed in state 'Solved'.
     */
    float getProbabilityValue(size_t voxel) const;
    tgt::vec2 getProbabilityRange() const;

    /// Returns whether the voxel at the passed index is defined to be a seed.
    bool isSeedPoint(size_t voxel) const;

    /// Returns the seed value of the specified voxel, or -1.0 if the voxel is not a seed.
    float getSeedValue(size_t voxel) const;

    /// Returns the total number of seed points.
    size_t getNumSeeds() const;

    /// Returns the minimum and maximum seed value.
    tgt::vec2 getSeedRange() const;

    /**
     * Returns the number of equations in the random walker equation system,
     * which equals the number of unseeded voxels.
     */
    size_t getSystemSize() const;

    /// Returns the dimensions of the volume the random walker solution is computed for.
    tgt::ivec3 getVolumeDimensions() const;

    /// Returns the number of voxel of the volume the random walker solution is computed for.
    size_t getNumVoxels() const;

    /**
     * Returns the index of the equation that represents the voxel with the specified id.
     * For seed voxels -1 is returned.
     */
    size_t getRowIndex(size_t voxel) const;

    /**
     * Generates a volume containing the computed probability values.
     * Operation only allowed in state 'Solved'.
     */
    template<class T>
    T* generateProbabilityVolume() const;

    /**
     * Generates a binary segmentation volume by thresholding the computed probability values.
     * Operation only allowed in state 'Solved'.
     */
    template<class T>
    T* generateBinarySegmentation(float threshold) const;

private:
    void computeVolIndexToRowMapping(const RandomWalkerSeeds* seeds);

    const VolumeBase* volume_;
    RandomWalkerSeeds* seeds_;
    RandomWalkerWeights* edgeWeights_;

    //const float* opacityBuffer_;
    //size_t opacityBufferSize_;

    EllpackMatrix<float> mat_;
    float* vec_;
    size_t* volIndexToRow_;
    float* solution_;

    tgt::ivec3 volDim_;
    size_t numVoxels_;
    tgt::vec3 volSpacing_;
    tgt::mat4 volTransformation_;

    size_t numSeeds_;

    SystemState state_;

    static const std::string loggerCat_;
};


//-----------------------------------------------------------------------------
// template definitions

template<class T>
T* RandomWalkerSolver::generateProbabilityVolume() const {

    if (state_ != Solved)
        throw VoreenException("System has not been solved");

    // check seed range for validity
    const tgt::vec2 seedRange = getSeedRange();
    if (seedRange.y - seedRange.x <= 0.f)  {
        std::ostringstream msg;
        msg << "Invalid seed range: " << seedRange;
        throw VoreenException(msg.str());
    }

    // create output volume
    T* result = 0;
    try {
        result = new T(volDim_/*, volSpacing_, volTransformation_*/); //FIXME
    }
    catch (std::bad_alloc&) {
        throw VoreenException("Bad allocation during creation of output volume");
    }

    // write normalized probability values to output volume
    const float seedRangeSizeRCP = 1.f / (seedRange.y - seedRange.x);
    for (size_t i=0; i<numVoxels_; i++) {
        float probValue = getProbabilityValue(i);
        result->setVoxelNormalized(tgt::clamp((probValue - seedRange.x) * seedRangeSizeRCP, 0.f, 1.f), i);
    }

    return result;
}

template<class T>
T* RandomWalkerSolver::generateBinarySegmentation(float threshold) const {

    if (state_ != Solved)
        throw VoreenException("System has not been solved");

    // create output volume
    T* result = 0;
    try {
        result = new T(volDim_/*, volSpacing_, volTransformation_*/);//FIXME
    }
    catch (std::bad_alloc&) {
        throw VoreenException("Bad allocation during creation of output volume");
    }

    // write normalized probability values to output volume
    for (size_t i=0; i<numVoxels_; i++) {
        if (getProbabilityValue(i) <= threshold)
            result->setVoxelNormalized(0.f, i);
        else
            result->setVoxelNormalized(1.f, i);
    }

    return result;
}


} //namespace

#endif
