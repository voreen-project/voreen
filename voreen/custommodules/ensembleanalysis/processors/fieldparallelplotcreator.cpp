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

#include "fieldparallelplotcreator.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/textureunit.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/utils/glsl.h"

#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

#include <random>

namespace voreen {

const std::string FieldParallelPlotCreator::loggerCat_("voreen.ensembleanalysis.FieldParallelPlotCreator");

const std::string FieldParallelPlotCreator::META_DATA_HASH("EnsembleHash");

FieldParallelPlotCreator::FieldParallelPlotCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , seedMask_(Port::INPORT, "seedmask", "Seed Mask Input (optional)")
    , outport_(Port::OUTPORT, "fpp.representation", "FieldPlotData Port")
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 8192, 1, 131072)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , verticalResolution_("verticalResolution", "Vertical Resolution", 128, 10, 16384)
    , horizontalResolutionPerTimeUnit_("horizontalResolutionPerTimeUnit", "Horizontal Resolution (Per Time Unit)", 10, 1, 100)
{    
    addPort(inport_);
    addPort(seedMask_);
    addPort(outport_);

    addProperty(numSeedPoints_);
    addProperty(seedTime_);
    addProperty(verticalResolution_);
    addProperty(horizontalResolutionPerTimeUnit_);
}

FieldParallelPlotCreator::~FieldParallelPlotCreator() {
}

Processor* FieldParallelPlotCreator::create() const {
    return new FieldParallelPlotCreator();
}

std::vector<std::reference_wrapper<Port>> FieldParallelPlotCreator::getCriticalPorts() {
    auto criticalPorts = AsyncComputeProcessor<ComputeInput, ComputeOutput>::getCriticalPorts();
    criticalPorts.erase(std::remove_if(criticalPorts.begin(), criticalPorts.end(),[this] (const std::reference_wrapper<Port>& port){
        return port.get().getID() == seedMask_.getID();
    }));
    return criticalPorts;
}

FieldParallelPlotCreatorInput FieldParallelPlotCreator::prepareComputeInput() {
    const EnsembleDataset* inputPtr = inport_.getThreadSafeData();
    if (!inputPtr)
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);

    if (inputPtr->getMaxNumTimeSteps() < 2)
        throw InvalidInputException("Num Time Steps is 1 or less, no need to aggregate over time.", InvalidInputException::S_WARNING);

    const EnsembleDataset& input = *inputPtr;

    const tgt::Bounds& roi = input.getRoi(); // ROI is defined in physical coordinates.
    if(!roi.isDefined()) {
        throw InvalidInputException("ROI is not defined", InvalidInputException::S_ERROR);
    }

    const VolumeBase* seedMask = seedMask_.getThreadSafeData();
    tgt::Bounds seedMaskBounds;
    tgt::mat4 seedMaskPhysicalToVoxelMatrix;
    std::unique_ptr<VolumeRAMRepresentationLock> seedMaskLock;
    if(seedMask) {
        tgt::Bounds roiBounds = roi;
        seedMaskBounds = seedMask->getBoundingBox(false).getBoundingBox(false);

        roiBounds.intersectVolume(seedMaskBounds);
        if(!roiBounds.isDefined()) {
            throw InvalidInputException("Seed Mask does not overlap with ensemble ROI", InvalidInputException::S_ERROR);
        }

        seedMaskLock.reset(new VolumeRAMRepresentationLock(seedMask));

        VolumeMinMax* vmm = seedMask->getDerivedData<VolumeMinMax>();
        if(vmm->getMinNormalized() == 0.0f && vmm->getMaxNormalized() == 0.0f) {
            throw InvalidInputException("Seed Mask is empty", InvalidInputException::S_ERROR);
        }

        seedMaskPhysicalToVoxelMatrix = seedMask->getPhysicalToVoxelMatrix();
        LINFO("Restricting seed points to volume mask");
    }

    size_t height = static_cast<size_t>(verticalResolution_.get());
    size_t width  = static_cast<size_t>(input.getMaxTotalDuration() * horizontalResolutionPerTimeUnit_.get()) + 1;
    size_t depth  = static_cast<size_t>(input.getCommonFieldNames().size() * input.getRuns().size());

    std::unique_ptr<FieldPlotData> plotData(new FieldPlotData(width, height, depth));

    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    const size_t maxTries = -1; // TODO: choose a user defined approach
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
                std::abs((*seedMaskLock)->getVoxelNormalized(seedMaskPhysicalToVoxelMatrix*seedPoint)) <
                    std::numeric_limits<float>::epsilon()));

        if(tries < maxTries) {
            seedPoints.push_back(seedPoint);
        }
    }

    return FieldParallelPlotCreatorInput{
            input,
            std::move(plotData),
            std::move(seedPoints),
    };
}
FieldParallelPlotCreatorOutput FieldParallelPlotCreator::compute(FieldParallelPlotCreatorInput input, ProgressReporter& progress) const {

    progress.setProgress(0.0f);

    const EnsembleDataset& data = input.dataset;
    std::unique_ptr<FieldPlotData> plotData = std::move(input.outputPlot);
    std::vector<tgt::vec3> seedPoints = std::move(input.seedPoints);
    const std::vector<std::string>& fieldNames = data.getCommonFieldNames();

    const int pixelPerTimeUnit = horizontalResolutionPerTimeUnit_.get();
    float timeOffset = -data.getStartTime();

    size_t sliceNumber = 0;
    for (size_t i=0; i<fieldNames.size(); i++) {

        const std::string& field = fieldNames[i];
        SubtaskProgressReporter fieldProgressReporter(progress, tgt::vec2(i, i+1) / tgt::vec2(fieldNames.size()));

        for (size_t j=0; j<data.getRuns().size(); j++) {

            const EnsembleDataset::Run& run = data.getRuns()[j];
            SubtaskProgressReporter runProgressReporter(fieldProgressReporter, tgt::vec2(j, j+1) / tgt::vec2(data.getRuns().size()));

            const tgt::vec2& valueRange = data.getValueRange(field);
            float pixelOffset = pixelPerTimeUnit * (timeOffset + run.timeSteps_[0].time_);
            float pixel = pixelPerTimeUnit * run.timeSteps_[0].duration_;

            const VolumeBase* volumePrev = run.timeSteps_[0].fieldNames_.at(field);
            tgt::mat4 physicalToVoxelMatrixPrev = volumePrev->getPhysicalToVoxelMatrix();
            RealWorldMapping rwmPrev = volumePrev->getRealWorldMapping();
            VolumeRAMRepresentationLock lockPrev(volumePrev);

            float progressPerTimeStep = 1.0f / (input.dataset.getTotalNumTimeSteps());
            for (size_t t = 1; t < run.timeSteps_.size(); t++) {

                const VolumeBase* volumeCurr = run.timeSteps_[t].fieldNames_.at(field);
                tgt::mat4 physicalToVoxelMatrixCurr = volumeCurr->getPhysicalToVoxelMatrix();
                RealWorldMapping rwmCurr = volumeCurr->getRealWorldMapping();
                VolumeRAMRepresentationLock lockCurr(volumeCurr);

                // Determine pixel positions.
                size_t x1 = static_cast<size_t>(pixelOffset);
                size_t x2 = static_cast<size_t>(pixelOffset + pixel);

                for (size_t k = 0; k<seedPoints.size(); k++) {

                    float voxelPrev = lockPrev->getVoxelNormalizedLinear(physicalToVoxelMatrixPrev * seedPoints[k]);
                    voxelPrev = rwmPrev.normalizedToRealWorld(voxelPrev);
                    voxelPrev = mapRange(voxelPrev, valueRange.x, valueRange.y, 0.0f, 1.0f);
                    float voxelCurr = lockCurr->getVoxelNormalizedLinear(physicalToVoxelMatrixCurr * seedPoints[k]);
                    voxelCurr = rwmCurr.normalizedToRealWorld(voxelCurr);
                    voxelCurr = mapRange(voxelCurr, valueRange.x, valueRange.y, 0.0f, 1.0f);

                    plotData->drawConnection(x1, x2, voxelPrev, voxelCurr, sliceNumber);
                }

                volumePrev = volumeCurr;
                physicalToVoxelMatrixPrev = physicalToVoxelMatrixCurr;
                rwmPrev = rwmCurr;
                lockPrev = volumePrev;

                pixelOffset = pixelOffset + pixel;
                pixel = pixelPerTimeUnit * run.timeSteps_[t].duration_;

                // Update progress.
                runProgressReporter.setProgress(t*progressPerTimeStep);
            }
            sliceNumber++;
        }
    }

    // Add ensemble hash.
    plotData->getVolume()->getMetaDataContainer().addMetaData(FieldParallelPlotCreator::META_DATA_HASH, new StringMetaData(EnsembleHash(data).getHash()));
    
    // We're done here.
    progress.setProgress(1.0f);
    return FieldParallelPlotCreatorOutput{std::move(plotData)};
}

void FieldParallelPlotCreator::processComputeOutput(FieldParallelPlotCreatorOutput output) {
    outport_.setData(output.plotData.release(), true);
}

bool FieldParallelPlotCreator::isReady() const {
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

void FieldParallelPlotCreator::adjustPropertiesToInput() {
    //TODO: implement heuristic for auto-selecting number of seed points depending on dataset
    //numSeedPoints_.set(32768);
}

} // namespace voreen
