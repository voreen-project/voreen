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

namespace voreen {

static const int MAX_NUM_DIMENSIONS = 3;

const std::string SimilarityMatrixCreator::loggerCat_("voreen.ensembleanalysis.SimilarityMatrixCreator");

SimilarityMatrixCreator::SimilarityMatrixCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "inport", "Ensemble Datastructure Input", false)
    , seedMask_(Port::INPORT, "seedmask", "Seed Mask Input (optional)")
    , outport_(Port::OUTPORT, "outport", "Similarity Matrix Output", false)
    , fieldSimilarityMeasure_("fieldSimilarityMeasure", "Field Similarity Measure")
    , isoValue_("isovalue", "Iso-Value", 0.5f, 0.0f, 1.0f)
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 0, 0, 0)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
{
    // Ports
    addPort(inport_);
    addPort(seedMask_);
    addPort(outport_);

    // Calculation
    addProperty(fieldSimilarityMeasure_);
    fieldSimilarityMeasure_.addOption("isovalue", "Iso-Surface", MEASURE_ISOSURFACE);
    fieldSimilarityMeasure_.addOption("multifield", "Multi-Field", MEASURE_MULTIFIELD);
    fieldSimilarityMeasure_.set("multifield");
    ON_CHANGE_LAMBDA(fieldSimilarityMeasure_, [this] {
        isoValue_.setVisibleFlag(fieldSimilarityMeasure_.getValue() == MEASURE_ISOSURFACE);
    });
    addProperty(isoValue_);
    isoValue_.setVisibleFlag(false);

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

void SimilarityMatrixCreator::adjustPropertiesToInput() {

    const EnsembleDataset* dataset = inport_.getData();
    if (!dataset)
        return;

    numSeedPoints_.setMinValue(1);
    numSeedPoints_.setMaxValue(131072);
    numSeedPoints_.set(32768);
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

    const VolumeBase* seedMask = seedMask_.getData();
    tgt::Bounds seedMaskBounds;
    tgt::mat4 seedMaskPhysicalToVoxelMatrix;
    if (seedMask) {
        seedMaskBounds = seedMask->getBoundingBox(false).getBoundingBox();
        seedMaskPhysicalToVoxelMatrix = seedMask->getPhysicalToVoxelMatrix();
        LINFO("Restricting seed points to volume mask");
    }

    std::unique_ptr<SimilarityMatrixList> outputMatrices(new SimilarityMatrixList(input));

    std::function<float()> rnd(
            std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    std::vector<tgt::vec3> seedPoints;
    for (int k = 0; k < numSeedPoints_.get(); k++) {
        tgt::vec3 seedPoint(rnd(), rnd(), rnd());
        seedPoint = tgt::vec3(roi.getLLF()) + seedPoint * tgt::vec3(roi.diagonal());

        // TODO: very rough and dirty restriction, implement something more intelligent.
        if (!seedMask || (seedMaskBounds.containsPoint(seedPoints[k]) &&
                          seedMask->getRepresentation<VolumeRAM>()->getVoxelNormalized(seedMaskPhysicalToVoxelMatrix*seedPoints[k]) != 0.0f)) {
            seedPoints.push_back(seedPoint);
        }
    }

    return SimilarityMatrixCreatorInput{
            input,
            std::move(outputMatrices),
            std::move(seedPoints),
            fieldSimilarityMeasure_.getValue(),
            isoValue_.get()
    };
}

SimilarityMatrixCreatorOutput SimilarityMatrixCreator::compute(SimilarityMatrixCreatorInput input, ProgressReporter& progress) const {

    std::unique_ptr<SimilarityMatrixList> similarityMatrices = std::move(input.outputMatrices);
    std::vector<tgt::vec3> seedPoints = std::move(input.seedPoints);

    // Init empty flags.
    std::vector<std::vector<float>> Flags(input.dataset.getTotalNumTimeSteps(), std::vector<float>(seedPoints.size(), 0.0f));

    progress.setProgress(0.0f);

    const std::vector<std::string>& channels = input.dataset.getCommonChannels();
    for (size_t i=0; i<channels.size(); i++) {
        const std::string& channel = channels[i];
        const tgt::vec2& valueRange = input.dataset.getValueRange(channel);

        SubtaskProgressReporter runProgressReporter(progress, tgt::vec2(i, 0.9f*(i+1))/tgt::vec2(channels.size()));
        size_t index = 0;
        for (const EnsembleDataset::Run& run : input.dataset.getRuns()) {
            float progressPerTimeStep = 1.0f / (run.timeSteps_.size() * input.dataset.getRuns().size());
            for (const EnsembleDataset::TimeStep& timeStep : run.timeSteps_) {

                const VolumeBase* volume = timeStep.channels_.at(channel);
                tgt::mat4 physicalToVoxelMatrix = volume->getPhysicalToVoxelMatrix();

                for (size_t k=0; k<seedPoints.size(); k++) {
                    float value = volume->getRepresentation<VolumeRAM>()->getVoxelNormalizedLinear(physicalToVoxelMatrix * seedPoints[k]);
                    value = mapRange(value, valueRange.x, valueRange.y, 0.0f, 1.0f);

                    switch(input.fieldSimilarityMeasure) {
                    case MEASURE_ISOSURFACE:
                    {
                        bool inside = value < input.isoValue;
                        Flags[index][k] = inside ? 1.0f : 0.0f;
                        break;
                    }
                    case MEASURE_MULTIFIELD:
                        Flags[index][k] = value;
                        break;
                    default:
                        break;
                    }
                }
                index++;

                // Update progress.
                runProgressReporter.setProgress(index * progressPerTimeStep);
            }
        }

        ////////////////////////////////////////////////////////////
        //Calculate distances for up-right corner and reflect them//
        ////////////////////////////////////////////////////////////

        SimilarityMatrix& DistanceMatrix = similarityMatrices->getSimilarityMatrix(channel);

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for shared(Flags)
#endif
        for (long i=0; i<static_cast<long>(DistanceMatrix.getSize()); i++) {
            for (long j=0; j<i+1; j++) {
                float ScaleSum = 0.0f;
                float resValue = 0.0f;

                for (size_t k=0; k<seedPoints.size(); k++) {

                    float a = Flags[i][k];
                    float b = Flags[j][k];

                    ScaleSum += (1.0f - (a < b ? a : b));
                    resValue += (1.0f - (a > b ? a : b));
                }

                if (ScaleSum > 0.0f)
                    DistanceMatrix(i, j) = (ScaleSum-resValue) / ScaleSum;
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
