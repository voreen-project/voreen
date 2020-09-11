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

#include "lambdacicriterion.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"


#include "Eigen/Eigenvalues"

namespace voreen {

const std::string LambdaCiCriterion::loggerCat_("voreen.flowsimulation.lambdacicriterion");

LambdaCiCriterion::LambdaCiCriterion()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inputVolume_(Port::INPORT, "lambdacicriterion.inputVolume", "Volume Input")
    , outputVolume_(Port::OUTPORT, "lambdacicriterion.outputVolume", "Volume Output")
{
    addPort(inputVolume_);
    inputVolume_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(outputVolume_);
}

Processor* LambdaCiCriterion::create() const {
    return new LambdaCiCriterion();
}

LambdaCiCriterionInput LambdaCiCriterion::prepareComputeInput() {

    auto inputVolume = inputVolume_.getThreadSafeData();
    if(!inputVolume) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<VolumeRAM_Float> outputVolume(new VolumeRAM_Float(inputVolume->getDimensions()));

    return LambdaCiCriterionInput{
        std::move(inputVolume),
        std::move(outputVolume)
    };
}

LambdaCiCriterionOutput LambdaCiCriterion::compute(LambdaCiCriterionInput input, ProgressReporter& progressReporter) const {

    auto inputVolume = std::move(input.inputVolume);
    VolumeRAMRepresentationLock inputVolumeData(inputVolume);
    auto outputVolume = std::move(input.outputVolume);
    outputVolume->clear();

    tgt::svec3 dimensions = inputVolume->getDimensions();
    tgt::vec3 spacing = inputVolume->getSpacing();

    auto jacobi = Eigen::Matrix3f();

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
    for (long z = 0; z < static_cast<long>(dimensions.z); z++) {
        for (size_t y = 0; y < dimensions.y; y++) {
            for (size_t x = 0; x < dimensions.x; x++) {
                tgt::svec3 pos(x, y, z);
                for(size_t i = 0; i < 3; i++) {
                    tgt::vec3 gradient = gradientFunc(*inputVolumeData, spacing, pos, i);
                    for(size_t j = 0; j < 3; j++) {
                        jacobi(i, j) = gradient[j];
                    }
                }

                const float epsilon = 0.00001f;

                auto eigenvalues = jacobi.eigenvalues();

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
                    //std::cout << "Only one real allowed!" << std::endl;
                    continue;
                }

                if(std::abs(eigenvalues(imag[0]).real() - eigenvalues(imag[1]).real()) > epsilon) {
                    //std::cout << "Real parts of imaginary eigenvaules differ!" << std::endl;
                    continue;
                }

                if(std::abs(eigenvalues(imag[0]).imag() + eigenvalues(imag[1]).imag()) > epsilon) {
                    //std::cout << "Image parts of imaginary eigenvaules differ!" << std::endl;
                    continue;
                }

                outputVolume->voxel(pos) = 1.0f / std::abs(eigenvalues(imag[0]).imag());
                //std::cout << "Found!" << std::endl;
            }
        }

        progress.reportStepDone();
    }

    std::unique_ptr<Volume> volume(new Volume(outputVolume.release(), inputVolume));
    volume->setMetaDataValue<StringMetaData>("name", "lambdaCi");

    progressReporter.setProgress(1.0f);

    return LambdaCiCriterionOutput{
        std::move(volume)
    };
}

void LambdaCiCriterion::processComputeOutput(LambdaCiCriterionOutput output) {
    outputVolume_.setData(output.volume.release());
}

} // namespace voreen
