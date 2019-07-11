/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "similaritymatrixcreator.h"

#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"

#include <random>

#ifdef VRN_MODULE_VESSELTOPOLOGY
#include "custommodules/vesseltopology/datastructures/diskarraystorage.h"
#endif

namespace voreen {

const std::string SimilarityMatrixCreator::loggerCat_("voreen.ensembleanalysis.SimilarityMatrixCreator");

SimilarityMatrixCreator::SimilarityMatrixCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "inport", "Ensemble Datastructure Input", false)
    , seedMask_(Port::INPORT, "seedmask", "Seed Mask Input (optional)")
    , outport_(Port::OUTPORT, "outport", "Similarity Matrix Output", false)
    , singleChannelSimilarityMeasure_("singleChannelSimilarityMeasure", "Single Field Similarity Measure")
    , isoValue_("isoValue", "Iso-Value", 0.5f, 0.0f, 1.0f)
    , multiChannelSimilarityMeasure_("multiChannelSimilarityMeasure", "Multi Field Similarity Measure")
    , angleThreshold_("angleThreshold", "Angle Threshold", 20.0f, 0.0f, 180.0f)
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 8192, 1, 131072)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
{
    // Ports
    addPort(inport_);
    addPort(seedMask_);
    addPort(outport_);

    // Calculation
    addProperty(singleChannelSimilarityMeasure_);
    singleChannelSimilarityMeasure_.addOption("isovalue", "Iso-Surface", MEASURE_ISOSURFACE);
    singleChannelSimilarityMeasure_.addOption("multifield", "Multi-Field", MEASURE_MULTIFIELD);
    singleChannelSimilarityMeasure_.set("multifield");
    ON_CHANGE_LAMBDA(singleChannelSimilarityMeasure_, [this] {
       isoValue_.setVisibleFlag(singleChannelSimilarityMeasure_.getValue() == MEASURE_ISOSURFACE);
    });

    addProperty(isoValue_);
    isoValue_.setVisibleFlag(false);

    addProperty(multiChannelSimilarityMeasure_);
    multiChannelSimilarityMeasure_.addOption("magnitude", "Magnitude", MEASURE_MAGNITUDE);
    multiChannelSimilarityMeasure_.addOption("angleDifference", "Angle Difference", MEASURE_ANGLEDIFFERENCE);
    multiChannelSimilarityMeasure_.addOption("crossproduct", "Crossproduct Magnitude", MEASURE_CROSSPRODUCT);
    ON_CHANGE_LAMBDA(multiChannelSimilarityMeasure_, [this] {
        angleThreshold_.setVisibleFlag(multiChannelSimilarityMeasure_.getValue() == MEASURE_ANGLEDIFFERENCE);
    });

    addProperty(angleThreshold_);
    angleThreshold_.setVisibleFlag(false);

    addProperty(numSeedPoints_);
    addProperty(seedTime_);
}

Processor* SimilarityMatrixCreator::create() const {
    return new SimilarityMatrixCreator();
}


bool SimilarityMatrixCreator::isReady() const {
    if (!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if (!inport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }

    // Note: Seed Mask is optional!

    return true;
}

std::vector<std::reference_wrapper<Port>> SimilarityMatrixCreator::getCriticalPorts() {
    auto criticalPorts = AsyncComputeProcessor<ComputeInput, ComputeOutput>::getCriticalPorts();
    criticalPorts.erase(std::remove_if(criticalPorts.begin(), criticalPorts.end(), [this] (const std::reference_wrapper<Port>& port){
       return port.get().getID() == seedMask_.getID();
    }));
    return criticalPorts;
}

void SimilarityMatrixCreator::adjustPropertiesToInput() {
    //TODO: implement heuristic for auto-selecting number of seed points depending on dataset
    //numSeedPoints_.set(32768);
}

SimilarityMatrixCreatorInput SimilarityMatrixCreator::prepareComputeInput() {
    const EnsembleDataset* inputPtr = inport_.getThreadSafeData();
    if (!inputPtr)
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);

    const EnsembleDataset& input = *inputPtr;

    const tgt::Bounds& roi = input.getRoi(); // ROI is defined in physical coordinates.
    if (!roi.isDefined()) {
        throw InvalidInputException("ROI is not defined", InvalidInputException::S_ERROR);
    }

    for(const std::string& fieldName : input.getCommonFieldNames()) {
        size_t numChannels = input.getNumChannels(fieldName);
        if(numChannels != 1 && numChannels != 3) {
            throw InvalidInputException("Only volumes with 1 or 3 channels supported", InvalidInputException::S_ERROR);
        }
    }

    const VolumeBase* seedMask = seedMask_.getThreadSafeData();
    tgt::Bounds seedMaskBounds;
    tgt::mat4 seedMaskPhysicalToVoxelMatrix;
    std::unique_ptr<VolumeRAMRepresentationLock> seedMaskLock;
    if (seedMask) {
        seedMaskLock.reset(new VolumeRAMRepresentationLock(seedMask));
        seedMaskBounds = seedMask->getBoundingBox(false).getBoundingBox();
        seedMaskPhysicalToVoxelMatrix = seedMask->getPhysicalToVoxelMatrix();
        LINFO("Restricting seed points to volume mask");
    }

    std::unique_ptr<SimilarityMatrixList> outputMatrices(new SimilarityMatrixList(input));

    std::function<float()> rnd(
            std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    size_t maxTries = 500; // TODO: choose a user defined approach
    std::vector<tgt::vec3> seedPoints;
    seedPoints.reserve(numSeedPoints_.get());
    for (int k = 0; k<numSeedPoints_.get(); k++) {
        tgt::vec3 seedPoint;

        size_t tries = 0;
        do {
            seedPoint = tgt::vec3(rnd(), rnd(), rnd());
            seedPoint = tgt::vec3(roi.getLLF()) + seedPoint * tgt::vec3(roi.diagonal());
            tries++;
        } while (tries < maxTries && seedMask && (!seedMaskBounds.containsPoint(seedPoint) ||
                (*seedMaskLock)->getVoxelNormalized(seedMaskPhysicalToVoxelMatrix*seedPoint) == 0.0f));

        if(tries < maxTries) {
            seedPoints.push_back(seedPoint);
        }
    }

    return SimilarityMatrixCreatorInput{
            input,
            std::move(outputMatrices),
            std::move(seedPoints),
            singleChannelSimilarityMeasure_.getValue(),
            isoValue_.get(),
            multiChannelSimilarityMeasure_.getValue(),
            angleThreshold_.get()
    };
}

SimilarityMatrixCreatorOutput SimilarityMatrixCreator::compute(SimilarityMatrixCreatorInput input, ProgressReporter& progress) const {

    std::unique_ptr<SimilarityMatrixList> similarityMatrices = std::move(input.outputMatrices);
    std::vector<tgt::vec3> seedPoints = std::move(input.seedPoints);

    progress.setProgress(0.0f);

    const std::vector<std::string>& fieldNames = input.dataset.getCommonFieldNames();
    for (size_t fi=0; fi<fieldNames.size(); fi++) {

        const std::string& fieldName = fieldNames[fi];
        const tgt::vec2& valueRange = input.dataset.getValueRange(fieldName);
        size_t numChannels = input.dataset.getNumChannels(fieldName);

        // Init empty flags.
#ifdef VRN_MODULE_VESSELTOPOLOGY
        DiskArrayStorage<float> Flags(VoreenApplication::app()->getUniqueTmpFilePath());
#else
        std::vector<float> Flags;
        Flags.reserve(input.dataset.getTotalNumTimeSteps() * seedPoints.size() * numChannels);
#endif

        SubtaskProgressReporter runProgressReporter(progress, tgt::vec2(fi, 0.9f*(fi+1))/tgt::vec2(fieldNames.size()));
        float progressPerTimeStep = 1.0f / (input.dataset.getTotalNumTimeSteps());
        size_t index = 0;
        for (const EnsembleDataset::Run& run : input.dataset.getRuns()) {
            for (const EnsembleDataset::TimeStep& timeStep : run.timeSteps_) {

                const VolumeBase* volume = timeStep.fieldNames_.at(fieldName);
                tgt::mat4 physicalToVoxelMatrix = volume->getPhysicalToVoxelMatrix();
                RealWorldMapping rwm = volume->getRealWorldMapping();

                VolumeRAMRepresentationLock lock(volume);

                for (const tgt::vec3& seedPoint : seedPoints) {
                    for(size_t ch = 0; ch < numChannels; ch++) {
                        float value = lock->getVoxelNormalizedLinear(
                                physicalToVoxelMatrix * seedPoint, ch);

                        value = rwm.normalizedToRealWorld(value);
                        value = mapRange(value, valueRange.x, valueRange.y, 0.0f, 1.0f);

#ifdef VRN_MODULE_VESSELTOPOLOGY
                        Flags.storeElement(value);
#else
                        Flags.push_back(value);
#endif
                    }
                }

                // Update progress.
                runProgressReporter.setProgress(index * progressPerTimeStep);
                index++;
            }
        }

        ////////////////////////////////////////////////////////////
        //Calculate distances for up-right corner and reflect them//
        ////////////////////////////////////////////////////////////

        auto calcIndex = [&seedPoints, numChannels] (size_t timeStepIndex, size_t seedIndex, size_t channel = 0) {
            return timeStepIndex * seedPoints.size() * numChannels + seedIndex * numChannels + channel;
        };

        SimilarityMatrix& DistanceMatrix = similarityMatrices->getSimilarityMatrix(fieldName);

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for shared(Flags), shared(DistanceMatrix)
#endif
        for (long i = 0; i < static_cast<long>(DistanceMatrix.getSize()); i++) {
            for (long j = 0; j <= i; j++) {

                float scaleSum = 0.0f; // Can be interpret as the number of (equal) samples
                float resValue = 0.0f; // Can be interpret as the number of differing samples

                if(numChannels == 1) {
                    for (size_t k = 0; k < seedPoints.size(); k++) {

                        float a = Flags[calcIndex(i, k)];
                        float b = Flags[calcIndex(j, k)];

                        if(input.singleChannelSimilarityMeasure == MEASURE_ISOSURFACE) {
                            a = a < input.isoValue ? 1.0f : 0.0f;
                            b = b < input.isoValue ? 1.0f : 0.0f;
                        }

                        scaleSum += (1.0f - (a < b ? a : b));
                        resValue += (1.0f - (a > b ? a : b));
                    }
                }
                else {

                    for (size_t k = 0; k < seedPoints.size(); k++) {

                        tgt::vec4 direction_i = tgt::vec4::zero;
                        tgt::vec4 direction_j = tgt::vec4::zero;

                        for (size_t ch = 0; ch < numChannels; ch++) {
                            direction_i[ch] = Flags[calcIndex(i, k, ch)];
                            direction_j[ch] = Flags[calcIndex(j, k, ch)];
                        }

                        if(input.multiChannelSimilarityMeasure == MEASURE_MAGNITUDE) {
                            float a = tgt::length(direction_i);
                            float b = tgt::length(direction_j);
                            float diff = std::abs(a - b);

                            scaleSum += 1.0f;
                            resValue += diff;

                            //scaleSum += (1.0f - (a < b ? a : b));
                            //resValue += (1.0f - (a > b ? a : b));
                        }
                        else if(input.multiChannelSimilarityMeasure == MEASURE_ANGLEDIFFERENCE) {
                            scaleSum += 1.0f; // add one more sample
                            if (direction_i == tgt::vec4::zero && direction_j == tgt::vec4::zero) {
                                resValue += 1.0f;
                            }
                            else if (direction_i != tgt::vec4::zero && direction_j != tgt::vec4::zero) {
                                tgt::vec4 normDirection_i = tgt::normalize(direction_i);
                                tgt::vec4 normDirection_j = tgt::normalize(direction_j);

                                float dot = tgt::dot(normDirection_i, normDirection_j);
                                float angle = tgt::rad2deg(std::acos(tgt::clamp(dot, -1.0f, 1.0f)));

                                if(angle < input.angleThreshold) {
                                    resValue += angle / 180.0f;
                                }
                            }
                            else {
                                //resValue += 0.0f;
                            }
                        }
                        else if(input.multiChannelSimilarityMeasure == MEASURE_CROSSPRODUCT) {
                            scaleSum += 1.0f; // add one more sample
                            if (direction_i == tgt::vec4::zero && direction_j == tgt::vec4::zero) {
                                resValue += 1.0f;
                            }
                            else if (direction_i != tgt::vec4::zero && direction_j != tgt::vec4::zero) {
                                float area = tgt::length(tgt::cross(direction_i.xyz(), direction_j.xyz()));
                                // In case area is 0, we have to account for colinear vectors.
                                if(area < std::numeric_limits<float>::epsilon()) {
                                    float length_i = tgt::length(direction_i);
                                    float length_j = tgt::length(direction_j);

                                    tgt::vec4 normDirection_i = direction_i / length_i;
                                    tgt::vec4 normDirection_j = direction_j / length_j;

                                    // Determine direction of collinearity.
                                    float dot = tgt::dot(normDirection_i, normDirection_j);
                                    float angle = std::acos(tgt::clamp(dot, -1.0f, 1.0f));
                                    if(angle > tgt::PIf*0.5f) {
                                        resValue += tgt::abs(-length_i - length_j);
                                    }
                                    else {
                                        resValue += tgt::abs(length_i - length_j);
                                    }
                                }
                                else {
                                    resValue += area;
                                }
                            }
                            else {
                                //resValue += 0.0f;
                            }
                        }
                    }
                }

                if (scaleSum > 0.0f)
                    DistanceMatrix(i, j) = (scaleSum - resValue) / scaleSum;
                else
                    DistanceMatrix(i, j) = 1.0f;
            }
        }
    }

    progress.setProgress(1.0f);

    return SimilarityMatrixCreatorOutput{
            std::move(similarityMatrices)
    };
}

void SimilarityMatrixCreator::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.outputMatrices.release(), true);
}


} // namespace
