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

#include "parallelvectors.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "tgt/matrix.h"

#include "Eigen/Eigenvalues"

#include <vector>

namespace voreen {

ParallelVectors::ParallelVectors()
    : AsyncComputeProcessor<ParallelVectorsInput, ParallelVectorsOutput>()
    , _inV(Port::INPORT, "in_v", "Vector field V")
    , _inW(Port::INPORT, "in_w", "Vector field W")
    , _inJacobian(Port::INPORT, "in_jacobi", "Jacobi Matrix (Optional)")
    , _inMask(Port::INPORT, "inportMask", "Mask (Optional)")
    , _out(Port::OUTPORT, "outport", "Parallel Vector Solution Points")
    , _sujudiHaimes("sujudiHaimes", "Use Sujudi-Haimes method for filtering", false)
{
    addPort(_inV);
    _inV.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(_inW);
    _inW.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(_inJacobian);
    _inJacobian.addCondition(new PortConditionVolumeType("Matrix3(float)", "Volume_Mat3Float"));
    ON_CHANGE(_inJacobian, ParallelVectors, onChangedJacobianData);
    addPort(_inMask);
    addPort(_out);

    addProperty(_sujudiHaimes);
}

void ParallelVectors::onChangedJacobianData() {
    _sujudiHaimes.setReadOnlyFlag(!_inJacobian.hasData());
}

bool ParallelVectors::isReady() const
{
    if(!_inV.isReady() || !_inW.isReady()) {
        setNotReadyErrorMessage("V and W must both be defined");
        return false;
    }

    // Note: Jacobian and mask are both optional.

    return true;
}

void ParallelVectors::Process(const VolumeRAM& v, const VolumeRAM& w, const VolumeRAM_Mat3Float* jacobian, const VolumeRAM* mask, ParallelVectorSolutions& outSolution, const RealWorldMapping& rwmV, const RealWorldMapping& rwmW, ProgressReporter* progress) {
    const auto dim = v.getDimensions();
    auto triangleSolutions = std::vector<tgt::vec3>();
    auto triangleSolutionIndices = std::vector<int32_t>((dim.x - 1) * (dim.y - 1) * (dim.z - 1) * TetrahedraPerCube * TrianglesPerTetrahedron, -1);

    // We are very pessimistic to check for complex numbers.
    const double epsilon = 0;

    const int32_t trianglesPerXInc = 24;
    const int32_t trianglesPerYInc = (dim.x - 1) * TetrahedraPerCube * TrianglesPerTetrahedron;
    const int32_t trianglesPerZInc = (dim.x - 1) * (dim.y - 1) * TetrahedraPerCube * TrianglesPerTetrahedron;
    const std::array<int32_t, 24> partnerTriangleOffsets
    {
        -trianglesPerZInc + 19, 20, 4, trianglesPerYInc + 5,
        -trianglesPerZInc + 11, 4, -4, trianglesPerXInc + 13,
        -trianglesPerYInc - 5, -4, 4, trianglesPerXInc + 5,
        -trianglesPerYInc + 11, 4, -4, trianglesPerZInc - 11,
        -trianglesPerXInc - 5, -4, 4, trianglesPerZInc - 19,
        -trianglesPerXInc - 13, -20, -4, trianglesPerYInc - 11
    };

    auto voxels = std::vector<tgt::svec3>();
    voxels.reserve(v.getNumVoxels());
    for (size_t x = 0; x < dim.x - 1; ++x) {
        for (size_t y = 0; y < dim.y - 1; ++y) {
            for (size_t z = 0; z < dim.z - 1; ++z) {
                if (!mask || mask->getVoxelNormalized(x, y, z) != 0.0f)
                    voxels.emplace_back(tgt::svec3(x, y, z));
            }
        }
    }
    voxels.shrink_to_fit();

    std::unique_ptr<ThreadedTaskProgressReporter> progressReporter;
    if(progress) {
        progressReporter.reset(new ThreadedTaskProgressReporter(*progress, voxels.size()));
    }
    bool aborted = false;

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (long voxelIndex = 0; voxelIndex < static_cast<long>(voxels.size()); ++voxelIndex) {
        if (aborted) {
            continue;
        }

        const auto x = voxels[voxelIndex].x;
        const auto y = voxels[voxelIndex].y;
        const auto z = voxels[voxelIndex].z;

        const std::array<const Tet, TetrahedraPerCube> cubeTets{
            Tet{                                                                              // front top left tet 0
                Triangle{tgt::svec3{x, y, z}, {x, y + 1, z}, {x + 1, y + 1, z}},              // front 0
                Triangle{tgt::svec3{x, y, z}, {x, y + 1, z}, {x + 1, y + 1, z + 1}},          // back 1
                Triangle{tgt::svec3{x, y, z}, {x + 1, y + 1, z}, {x + 1, y + 1, z + 1}},      // right 2
                Triangle{tgt::svec3{x, y + 1, z}, {x + 1, y + 1, z}, {x + 1, y + 1, z + 1}}}, // top 3
                Tet{                                                                              // front bottom right tet 1
                Triangle{tgt::svec3{x, y, z}, {x + 1, y, z}, {x + 1, y + 1, z}},              // front 4
                Triangle{tgt::svec3{x, y, z}, {x + 1, y, z}, {x + 1, y + 1, z + 1}},          // back 5
                Triangle{tgt::svec3{x, y, z}, {x + 1, y + 1, z}, {x + 1, y + 1, z + 1}},      // left 6
                Triangle{tgt::svec3{x + 1, y, z}, {x + 1, y + 1, z}, {x + 1, y + 1, z + 1}}}, // right 7
                Tet{                                                                              // middle right bottom tet 2
                Triangle{tgt::svec3{x, y, z}, {x + 1, y, z}, {x + 1, y, z + 1}},              // bottom 8
                Triangle{tgt::svec3{x, y, z}, {x + 1, y, z}, {x + 1, y + 1, z + 1}},          // front 9
                Triangle{tgt::svec3{x, y, z}, {x + 1, y, z + 1}, {x + 1, y + 1, z + 1}},      // left 10
                Triangle{tgt::svec3{x + 1, y, z}, {x + 1, y, z + 1}, {x + 1, y + 1, z + 1}}}, // right 11
                Tet{                                                                              // back right bottom tet 3
                Triangle{tgt::svec3{x, y, z}, {x, y, z + 1}, {x + 1, y, z + 1}},              // bottom 12
                Triangle{tgt::svec3{x, y, z}, {x, y, z + 1}, {x + 1, y + 1, z + 1}},          // left 13
                Triangle{tgt::svec3{x, y, z}, {x + 1, y, z + 1}, {x + 1, y + 1, z + 1}},      // right 14
                Triangle{tgt::svec3{x, y, z + 1}, {x + 1, y, z + 1}, {x + 1, y + 1, z + 1}}}, // back 15
                Tet{                                                                              // back left bottom tet 4
                Triangle{tgt::svec3{x, y, z}, {x, y, z + 1}, {x, y + 1, z + 1}},              // left 16
                Triangle{tgt::svec3{x, y, z}, {x, y, z + 1}, {x + 1, y + 1, z + 1}},          // right 17
                Triangle{tgt::svec3{x, y, z}, {x, y + 1, z + 1}, {x + 1, y + 1, z + 1}},      // front 18
                Triangle{tgt::svec3{x, y, z + 1}, {x, y + 1, z + 1}, {x + 1, y + 1, z + 1}}}, // back 19
                Tet{                                                                              // middle left top tet 5
                Triangle{tgt::svec3{x, y, z}, {x, y + 1, z}, {x, y + 1, z + 1}},              // left 20
                Triangle{tgt::svec3{x, y, z}, {x, y + 1, z}, {x + 1, y + 1, z + 1}},          // front 21
                Triangle{tgt::svec3{x, y, z}, {x, y + 1, z + 1}, {x + 1, y + 1, z + 1}},      // back 22
                Triangle{tgt::svec3{x, y + 1, z}, {x, y + 1, z + 1}, {x + 1, y + 1, z + 1} /* top 23 */}} };


        auto applyRwm = [&cubeTets] (size_t tetIndexInCube, size_t triIndexInTet, const VolumeRAM& v, const RealWorldMapping& rwm) {
            const auto triangle = cubeTets[tetIndexInCube][triIndexInTet];
            const auto vol1Voxel0 = Eigen::Vector3d(rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[0], 0)), rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[0], 1)), rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[0], 2)));
            const auto vol1Voxel1 = Eigen::Vector3d(rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[1], 0)), rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[1], 1)), rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[1], 2)));
            const auto vol1Voxel2 = Eigen::Vector3d(rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[2], 0)), rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[2], 1)), rwm.normalizedToRealWorld(v.getVoxelNormalized(triangle[2], 2)));
            return (Eigen::Matrix3d() << vol1Voxel0, vol1Voxel1, vol1Voxel2).finished();
        };

        for (size_t tetIndexInCube = 0; tetIndexInCube < cubeTets.size(); ++tetIndexInCube) {
            for (size_t triIndexInTet = 0; triIndexInTet < TrianglesPerTetrahedron; ++triIndexInTet) {
                const auto triangleSolutionIndex = TrianglesPerTetrahedron * (TetrahedraPerCube * ((dim.x - 1) * ((dim.y - 1) * z + y) + x) + tetIndexInCube) + triIndexInTet;
                const auto partnerTriangleSolutionIndex = triangleSolutionIndex + partnerTriangleOffsets[tetIndexInCube * TrianglesPerTetrahedron + triIndexInTet];

                // We already have a solution.
                if (triangleSolutionIndices[triangleSolutionIndex] != -1) {
                    continue;
                }
                if (partnerTriangleSolutionIndex >= 0 && partnerTriangleSolutionIndex < triangleSolutionIndices.size()) {
                    triangleSolutionIndices[partnerTriangleSolutionIndex] = -2;
                }

                //auto doBreak = false;
                //for (auto i = 0; i < 3; ++i) { // Check if vectors are too small at triangle vertices 
                //    if (tgt::lengthSq(V.voxel(triangle[i])) < 0.00000000000001 || tgt::lengthSq(W.voxel(triangle[i])) < 0.00000000000001) {
                //        doBreak = true;
                //        break;
                //    }
                //}
                //if (doBreak) {
                //    break;
                //}

                const auto V = applyRwm(tetIndexInCube, triIndexInTet, v, rwmV);
                const auto W = applyRwm(tetIndexInCube, triIndexInTet, w, rwmW);

                Eigen::Matrix3d M, inverse;
                double ignore_det;
                bool invertible;

                V.computeInverseAndDetWithCheck(inverse, ignore_det, invertible, 0.00000001);
                if (invertible)
                    M = inverse * W;
                else
                {
                    W.computeInverseAndDetWithCheck(inverse, ignore_det, invertible, 0.00000001);
                    if (invertible)
                        M = inverse * V;
                    else
                    {
                        continue;
                    }
                }

                // V or W is invertible => find a solution now

                const auto eigensolver = Eigen::EigenSolver<Eigen::Matrix3d>(M);
                const auto& eigenvalues = eigensolver.eigenvalues();
                auto eigenvectors = eigensolver.eigenvectors();

                for (int eigenVectorIndex = 0; eigenVectorIndex < 3; ++eigenVectorIndex) {
                    const auto isReal = std::abs(eigenvalues(eigenVectorIndex).imag()) <= epsilon;
                    const auto sameSign = std::signbit(eigenvectors.col(eigenVectorIndex).x().real()) == std::signbit(eigenvectors.col(eigenVectorIndex).y().real()) && std::signbit(eigenvectors.col(eigenVectorIndex).y().real()) == std::signbit(eigenvectors.col(eigenVectorIndex).z().real());

                    if (!isReal || !sameSign)
                        continue;

                    eigenvectors.col(eigenVectorIndex) /= eigenvectors.col(eigenVectorIndex).sum();

                    auto addSolution = !jacobian;

                    if (jacobian) {
                        // Interpolate jacobian at solution (barycentric coordinates)
                        const auto jacobianAtSolution =
                                static_cast<float>(eigenvectors.col(eigenVectorIndex).x().real()) * jacobian->voxel(cubeTets[tetIndexInCube][triIndexInTet][0]) +
                                static_cast<float>(eigenvectors.col(eigenVectorIndex).y().real()) * jacobian->voxel(cubeTets[tetIndexInCube][triIndexInTet][1]) +
                                static_cast<float>(eigenvectors.col(eigenVectorIndex).z().real()) * jacobian->voxel(cubeTets[tetIndexInCube][triIndexInTet][2]);

                        using EigenMat3fRowMajor = Eigen::Matrix<float, 3, 3, Eigen::StorageOptions::RowMajor | Eigen::StorageOptions::AutoAlign>;
                        const auto jacobianEigenvalues = Eigen::EigenSolver<EigenMat3fRowMajor>(reinterpret_cast<const EigenMat3fRowMajor&>(jacobianAtSolution), false).eigenvalues();

                        int numberOfComplexEigenvalues = 0;
                        for (int j = 0; j < 3; ++j) { // count complex eigenvalues of Jacobian
                            if (std::abs(jacobianEigenvalues(j).imag()) > epsilon) {
                                ++numberOfComplexEigenvalues;
                            }
                        }
                        addSolution = numberOfComplexEigenvalues == 2;
                    }

                    if (addSolution) {
                        const auto pos =
                                static_cast<float>(eigenvectors.col(eigenVectorIndex).x().real()) * tgt::vec3(cubeTets[tetIndexInCube][triIndexInTet][0]) +
                                static_cast<float>(eigenvectors.col(eigenVectorIndex).y().real()) * tgt::vec3(cubeTets[tetIndexInCube][triIndexInTet][1]) +
                                static_cast<float>(eigenvectors.col(eigenVectorIndex).z().real()) * tgt::vec3(cubeTets[tetIndexInCube][triIndexInTet][2]);

#ifdef VRN_MODULE_OPENMP
                        #pragma omp critical
#endif
                        {
                            triangleSolutionIndices[triangleSolutionIndex] = triangleSolutions.size();
                            if (partnerTriangleSolutionIndex >= 0 && partnerTriangleSolutionIndex < triangleSolutionIndices.size())
                                triangleSolutionIndices[partnerTriangleSolutionIndex] = triangleSolutions.size();
                            triangleSolutions.push_back(pos);
                        }
                    }

                    if(progressReporter && progressReporter->reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
                        #pragma omp critical
                        aborted = true;
#else
                        aborted = true;
                        break;
#endif
                    }
                    break;
                }
            }
        }
    }

    if (aborted) {
        throw boost::thread_interrupted();
    }

    outSolution.dimensions = dim;
    outSolution.solutions.swap(triangleSolutions);
    outSolution.triangleSolutionIndices.swap(triangleSolutionIndices);
}

ParallelVectorsInput ParallelVectors::prepareComputeInput() {

    if(!_inV.isReady() || !_inW.isReady()) {
        throw InvalidInputException("V and W must both be defined", InvalidInputException::S_ERROR);
    }

    const auto* volumeV = _inV.getData();
    const auto* volumeW = _inW.getData();
    const auto* volumeJacobi = _inJacobian.getData();
    const auto* mask = _inMask.getData();

    if (volumeV->getDimensions() != volumeW->getDimensions()
        || (volumeJacobi && volumeV->getDimensions() != volumeJacobi->getDimensions())
        || (mask && volumeV->getDimensions() != mask->getDimensions()))
    {
        throw InvalidInputException("Input dimensions do not match", InvalidInputException::S_ERROR);
    }

    const auto dim = volumeV->getDimensions();
    if (std::min({dim.x, dim.y, dim.z}) < 2) {
        throw InvalidInputException("Input dimensions must be greater than 1 in each dimension", InvalidInputException::S_ERROR);
    }

    return {
            _inV.getThreadSafeData(),
            _inW.getThreadSafeData(),
            _inJacobian.getThreadSafeData(),
            _inMask.getThreadSafeData(),
            _sujudiHaimes.get()
    };
}
ParallelVectorsOutput ParallelVectors::compute(ParallelVectorsInput input, ProgressReporter& progressReporter) const {

    // Check for optional jacobian matrix.
    std::unique_ptr<VolumeRAMRepresentationLock> volumeJacobian;
    if(input.sujudiHaimes && input.jacobian) {
        volumeJacobian.reset(new VolumeRAMRepresentationLock(input.jacobian));
    }
    const auto* jacobian = volumeJacobian ? dynamic_cast<const VolumeRAM_Mat3Float*>(**volumeJacobian) : nullptr;

    // Check for optional mask.
    std::unique_ptr<VolumeRAMRepresentationLock> volumeMask;
    if(input.mask) {
        volumeMask.reset(new VolumeRAMRepresentationLock(input.mask));
    }
    const auto* mask = volumeMask ? **volumeMask : nullptr;

    // Calculate solutions.
    auto solutions = std::unique_ptr<ParallelVectorSolutions>( new ParallelVectorSolutions() );
    ParallelVectors::Process(**VolumeRAMRepresentationLock(input.v),
                             **VolumeRAMRepresentationLock(input.w),
                             jacobian, mask, *solutions,
                             input.v->getRealWorldMapping(),
                             input.w->getRealWorldMapping(),
                             &progressReporter
    );
    solutions->voxelToWorldMatrix = _inV.getData()->getVoxelToWorldMatrix();

    return { std::move(solutions) };
}
void ParallelVectors::processComputeOutput(ParallelVectorsOutput output) {
    _out.setData(output.solutions.release());
}

} // namespace voreen
