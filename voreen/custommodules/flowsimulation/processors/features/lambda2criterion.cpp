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

#include "lambda2criterion.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"


#include "Eigen/Eigenvalues"

namespace voreen {

const std::string Lambda2Criterion::loggerCat_("voreen.flowsimulation.lambda2criterion");

Lambda2Criterion::Lambda2Criterion()
        : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
        , inputVolume_(Port::INPORT, "lambda2criterion.inputVolume", "Volume Input")
        , outputVolume_(Port::OUTPORT, "lambda2criterion.outputVolume", "Volume Output")
{
    addPort(inputVolume_);
    inputVolume_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(outputVolume_);
}

Processor* Lambda2Criterion::create() const {
    return new Lambda2Criterion();
}

Lambda2CriterionInput Lambda2Criterion::prepareComputeInput() {

    auto inputVolume = inputVolume_.getThreadSafeData();
    if(!inputVolume) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<VolumeRAM_Float> outputVolume(new VolumeRAM_Float(inputVolume->getDimensions()));

    return Lambda2CriterionInput{
            std::move(inputVolume),
            std::move(outputVolume)
    };
}

Lambda2CriterionOutput Lambda2Criterion::compute(Lambda2CriterionInput input, ProgressReporter& progressReporter) const {

    auto inputVolume = std::move(input.inputVolume);
    VolumeRAMRepresentationLock inputVolumeData(inputVolume);

    auto outputVolume = std::move(input.outputVolume);
    outputVolume->clear();

    tgt::svec3 dimensions = inputVolumeData->getDimensions();
    tgt::vec3 spacing = inputVolume->getSpacing();

    // Very ugly, but allows to easily switch to different gradient function.
    //VolumeOperatorGradient operatorGradient;
    auto gradientFunc = std::bind(
            static_cast<tgt::vec3(&)(const VolumeRAM*, const tgt::vec3&, const tgt::svec3&, size_t)>(
                        VolumeOperatorGradient::calcGradientCentralDifferences
                    ),
            std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4
    );

    ThreadedTaskProgressReporter progress(progressReporter, dimensions.z);

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (size_t z = 0; z < dimensions.z; z++) {
        for (size_t y = 0; y < dimensions.y; y++) {
            for (size_t x = 0; x < dimensions.x; x++) {
                auto jacobi = Eigen::Matrix3f();
                tgt::svec3 pos(x, y, z);
                for(size_t i = 0; i < 3; i++) {
                    tgt::vec3 gradient = gradientFunc(*inputVolumeData, spacing, pos, i);
                    for(size_t j = 0; j < 3; j++) {
                        jacobi(i, j) = gradient[j];
                    }
                }

                const auto jacobiT = jacobi.transpose();
                const auto S = 0.5f * ( jacobi + jacobiT );
                const auto O = 0.5f * ( jacobi - jacobiT );

                auto eigenvalues = ( S * S + O * O ).eigenvalues();
                if(eigenvalues.x().real() < eigenvalues.y().real()) {
                    if (eigenvalues.y().real() < eigenvalues.z().real()) {
                        outputVolume->voxel(pos) = eigenvalues.y().real();
                    } else if (eigenvalues.x().real() < eigenvalues.z().real()) {
                        outputVolume->voxel(pos) = eigenvalues.z().real();
                    } else {
                        outputVolume->voxel(pos) = eigenvalues.x().real();
                    }
                }
                else {
                    if (eigenvalues.x().real() < eigenvalues.z().real()) {
                        outputVolume->voxel(pos) = eigenvalues.x().real();
                    }
                    else if (eigenvalues.y().real() < eigenvalues.z().real()) {
                        outputVolume->voxel(pos) = eigenvalues.z().real();
                    }
                    else {
                        outputVolume->voxel(pos) = eigenvalues.y().real();
                    }
                }

                // Negative values indicate vortex region, so remove positive regions.
                if(outputVolume->voxel(pos) > 0.0f) {
                    outputVolume->voxel(pos) = 0.0f;
                }
            }
        }

        progress.reportStepDone();
    }

    std::unique_ptr<Volume> volume(new Volume(outputVolume.release(), inputVolume));
    volume->setMetaDataValue<StringMetaData>("name", "lambda2");

    progressReporter.setProgress(1.0f);

    return Lambda2CriterionOutput{
            std::move(volume)
    };
}

void Lambda2Criterion::processComputeOutput(Lambda2CriterionOutput output) {
    outputVolume_.setData(output.volume.release());
}

} // namespace voreen
