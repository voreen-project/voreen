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

#include "vortexcriterion.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "Eigen/Eigenvalues"

namespace {

// --- Compute Delta Criterion --- //
float criterionDelta (const tgt::mat3& jacobian) {
    auto J = Eigen::Matrix3f::Map(jacobian.elem);
    auto eigenvalues = J.eigenvalues();
    return std::max({eigenvalues.x().imag(), eigenvalues.y().imag(), eigenvalues.z().imag()});
}

// --- Compute Q Criterion --- //
float criterionQ (const tgt::mat3& jacobian) {

    auto J = Eigen::Matrix3f::Map(jacobian.elem);

    // Strain Rate and Vorticity //
    const auto J_t = J.transpose();
    const auto S = 0.5f * (J + J_t);
    const auto O = 0.5f * (J - J_t);

    auto eigenvalues = (S.transpose() * S).eigenvalues();
    const auto spectralNormS = std::max({eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real()});
    eigenvalues = (O.transpose() * O).eigenvalues();
    const auto spectralNormO = std::max({eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real()});

    float q = 0.5f * (spectralNormO - spectralNormS);

    // Discard negative values.
    return std::max(0.0f, q);
}

// --- Compute Lambda2 Criterion --- //
float criterionLambda2 (const tgt::mat3& jacobian) {

    auto J = Eigen::Matrix3f::Map(jacobian.elem);

    // Strain Rate and Vorticity //
    const auto J_t = J.transpose();
    const auto S = 0.5f * (J + J_t);
    const auto O = 0.5f * (J - J_t);

    float lambda2 = 0.0f;

    auto eigenvalues = (S * S + O * O).eigenvalues();
    if (eigenvalues.x().real() < eigenvalues.y().real()) {
        if (eigenvalues.y().real() < eigenvalues.z().real()) {
            lambda2 = eigenvalues.y().real();
        } else if (eigenvalues.x().real() < eigenvalues.z().real()) {
            lambda2 = eigenvalues.z().real();
        } else {
            lambda2 = eigenvalues.x().real();
        }
    }
    else {
        if (eigenvalues.x().real() < eigenvalues.z().real()) {
            lambda2 = eigenvalues.x().real();
        } else if (eigenvalues.y().real() < eigenvalues.z().real()) {
            lambda2 = eigenvalues.z().real();
        } else {
            lambda2 = eigenvalues.y().real();
        }
    }

    // Discard positive values and invert values for MIP.
    if(lambda2 > 0.0f) {
        return 0.0f;
    }
    else {
        return -lambda2;
    }
}

// --- Compute LambdaCi Criterion --- //
float criterionLambdaCi (const tgt::mat3& jacobian) {

    const float epsilon = 0; // Be quite pessimistic..

    auto eigenvalues = Eigen::Matrix3f::Map(jacobian.elem).eigenvalues();

    std::vector<int> real;
    std::vector<int> imag;
    for(size_t i=0; i<3; i++) {
        if(std::abs(eigenvalues(i).imag()) > epsilon) {
            imag.push_back(i);
        }
        else {
            real.push_back(i);
        }
    }

    if(real.size() != 1) {
        return 0.0f;
    }

    if(std::abs(eigenvalues(imag[0]).real() - eigenvalues(imag[1]).real()) > epsilon) {
        return 0.0f;
    }

    if(std::abs(eigenvalues(imag[0]).imag() + eigenvalues(imag[1]).imag()) > epsilon) {
        return 0.0f;
    }

    return 1.0f / std::abs(eigenvalues(imag[0]).imag());
}

}

namespace voreen {

const std::string VortexCriterion::loggerCat_("voreen.flowanalysis.VortexCriterion");

VortexCriterion::VortexCriterion()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inputJacobian_(Port::INPORT, "VortexCriterion.inputJacobian", "Jacobi Matrix")
    , outputCriterion_(Port::OUTPORT, "VortexCriterion.outputCriterion", "Vortex Criterion")
    , criterion_("criterion", "Criterion")
{
    addPort(inputJacobian_);
    inputJacobian_.addCondition(new PortConditionVolumeType("Matrix3(float)", "Volume_Mat3Float"));
    addPort(outputCriterion_);

    addProperty(criterion_);
    criterion_.addOption("Q", "Q");
    criterion_.addOption("Delta", "Delta");
    criterion_.addOption("Lambda2", "Lambda2");
    //criterion_.addOption("LambdaCi", "LambdaCi"); // Needs some polishing.
}

Processor* VortexCriterion::create() const {
    return new VortexCriterion();
}


VortexCriterionIO VortexCriterion::prepareComputeInput() {

    auto inputJacobian = inputJacobian_.getThreadSafeData();
    if(!inputJacobian) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    const tgt::svec3 dim = inputJacobian->getDimensions();

    std::unique_ptr<VolumeRAM_Float> outputCriterion(new VolumeRAM_Float( dim ));

    std::function<float(const tgt::mat3&)> criterion;
    if(criterion_.get() == "Q") {
        criterion = criterionQ;
    }
    else if(criterion_.get() == "Delta") {
        criterion = criterionDelta;
    }
    else if(criterion_.get() == "Lambda2") {
        criterion = criterionLambda2;
    }
    else if(criterion_.get() == "LambdaCi") {
        criterion = criterionLambdaCi;
    }
    else {
        tgtAssert(false, "unhandled criterion");
        throw InvalidInputException("unhandled criterion", InvalidInputException::S_ERROR);
    }

    return VortexCriterionIO{
            std::move(inputJacobian),
            std::move(outputCriterion),
            criterion
    };
}

VortexCriterionIO VortexCriterion::compute(VortexCriterionIO input, ProgressReporter& progressReporter) const {

    VolumeRAMRepresentationLock inputVolumeData(input.inputJacobian);
    const auto* jacobian = dynamic_cast<const VolumeRAM_Mat3Float*>(*inputVolumeData);
    tgt::svec3 dimensions = inputVolumeData->getDimensions();

    ThreadedTaskProgressReporter progress(progressReporter, dimensions.z);
    bool aborted = false;

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (long z = 0; z < static_cast<long>(dimensions.z); z++) {
        if (aborted) {
            continue;
        }
        for (size_t y = 0; y < dimensions.y; y++) {
            for (size_t x = 0; x < dimensions.x; x++) {
                tgt::svec3 pos(x, y, z);
                input.outputCriterion->voxel(pos) = input.criterion(jacobian->voxel(pos));
            }
        }

        if (progress.reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
#pragma omp critical
            aborted = true;
#else
            aborted = true;
            break;
#endif
        }
    }

    if (aborted) {
        throw boost::thread_interrupted();
    }

    progressReporter.setProgress(1.0f);

    return input;
}

void VortexCriterion::processComputeOutput(VortexCriterionIO output) {
    tgt::vec3 spacing = output.inputJacobian->getSpacing();
    tgt::vec3 offset = output.inputJacobian->getOffset();

    auto* volume = new Volume(output.outputCriterion.release(), spacing, offset);
    volume->setModality(Modality(criterion_.get()));
    outputCriterion_.setData(volume);
}

} // namespace voreen
