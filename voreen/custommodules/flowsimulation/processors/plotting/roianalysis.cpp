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

#include "roianalysis.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

#include "../../utils/utils.h"

namespace voreen {

const std::string RoiAnalysis::loggerCat_("voreen.flowsimulation.RoiAnalysis");

RoiAnalysis::RoiAnalysis()
    : AsyncComputeProcessor()
    , volumeListPort_(Port::INPORT, "input.volumelist", "Volume List Port")
    , maskPort_(Port::INPORT, "input.parameter", "Parameter Port")
    , outport_(Port::OUTPORT,"output.plot", "Plot Port")
    , outputQuantity_("outputQuantity", "Output Quantity")
{
    addPort(volumeListPort_);
    volumeListPort_.addCondition(new PortConditionVolumeListEnsemble());
    addPort(maskPort_);
    addPort(outport_);

    addProperty(outputQuantity_);
    outputQuantity_.addOption("maxMagnitude", "Max Magnitude");
    outputQuantity_.addOption("meanMagnitude", "Mean Magnitude");
    outputQuantity_.addOption("meanComponents", "Mean Components");
    outputQuantity_.addOption("integration", "Integration");
    outputQuantity_.setGroupID("output");
}

RoiAnalysisInput RoiAnalysis::prepareComputeInput() {

    auto volumes = volumeListPort_.getThreadSafeData();
    if(!volumes || volumes->empty()) {
        throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
    }
    auto maskVolume = maskPort_.getData();
    if(!maskVolume) {
        throw InvalidInputException("No mask", InvalidInputException::S_IGNORE);
    }

    VolumeRAMRepresentationLock mask(maskVolume);
    tgt::svec3 dim = mask->getDimensions();
    tgt::mat4 voxelToWorldMatrix = maskVolume->getVoxelToWorldMatrix();

    std::vector<tgt::vec3> seedPoints;
    for(size_t z=0; z<dim.z; z++) {
        for(size_t y=0; y<dim.y; y++) {
            for(size_t x=0; x<dim.x; x++) {
                seedPoints.push_back(voxelToWorldMatrix * tgt::vec3(x, y, z));
            }
        }
    }

    const VolumeBase* reference = volumes->first();
    const size_t numChannels = reference->getNumChannels();
    if(numChannels != 1 && numChannels != 3) {
        throw InvalidInputException("Only 1 and 3 channel volumes supported", InvalidInputException::S_ERROR);
    }

    std::function<std::vector<float>(const std::vector<tgt::vec3>&)> outputFunc;
    std::unique_ptr<PlotData> output;

    if(outputQuantity_.get() == "meanComponents") {
        outputFunc = [numChannels] (const std::vector<tgt::vec3>& samples) {
            tgt::vec3 mean = tgt::vec3::zero;
            for(const auto& sample : samples) {
                mean += sample;
            }

            std::vector<float> output(numChannels);
            for(size_t i=0; i<numChannels; i++) {
                output[i] = mean[i] / static_cast<float>(samples.size());
            }
            return output;
        };

        output.reset(new PlotData(1, numChannels));
        if(numChannels == 1) {
            output->setColumnLabel(1, "mean");
        }
        else if(numChannels == 3) {
            output->setColumnLabel(1, "x");
            output->setColumnLabel(2, "y");
            output->setColumnLabel(3, "z");
        }
    }
    else {

        output.reset(new PlotData(1, 1));
        output->setColumnLabel(0, "Time [s]");
        output->setColumnLabel(1, "Mask");

        if (outputQuantity_.get() == "maxMagnitude") {
            outputFunc = [](const std::vector<tgt::vec3>& samples) {
                float maxMagnitudeSq = 0.0f;
                for (const auto& sample : samples) {
                    maxMagnitudeSq = std::max(maxMagnitudeSq, tgt::lengthSq(sample));
                }
                float magnitude = std::sqrt(maxMagnitudeSq);
                return std::vector<float>(1, magnitude);
            };
        }
        else if (outputQuantity_.get() == "meanMagnitude") {
            outputFunc = [](const std::vector<tgt::vec3>& samples) {
                float meanMagnitude = 0.0f;
                for (const auto& sample : samples) {
                    meanMagnitude += tgt::length(sample) / samples.size();
                }
                return std::vector<float>(1, meanMagnitude);
            };
        }
        else if(outputQuantity_.get() == "integration") {
            tgt::vec3 spacing = maskVolume->getSpacing();
            float volumeSize = seedPoints.size() * tgt::hmul(spacing);

            outputFunc = [volumeSize](const std::vector<tgt::vec3>& samples) {
                float accumulatedSum = 0.0f;
                for (const auto& sample : samples) {
                    accumulatedSum += tgt::length(sample);
                }
                float integration = accumulatedSum * volumeSize;
                return std::vector<float>(1, integration);
            };
        }
        else {
            tgtAssert(false, "unhandled output quantity");
        }
    }

    return RoiAnalysisInput {
            std::move(volumes),
            std::move(seedPoints),
            std::move(output),
            std::move(outputFunc),
    };
}

RoiAnalysisOutput RoiAnalysis::compute(RoiAnalysisInput input, ProgressReporter& progressReporter) const {

    auto volumes = std::move(input.volumes);
    auto seedPoints = std::move(input.seedPoints);
    std::unique_ptr<PlotData> data = std::move(input.output);
    auto outputFunc = std::move(input.outputFunc);

    progressReporter.setProgress(0.0f);

    bool timeSeries = isTimeSeries(volumes);
    if(timeSeries) {
        data->setColumnLabel(0, "Time [s]");
    }
    else {
        data->setColumnLabel(0, "Time Step");
    }

    for(size_t t = 0; t < volumes->size(); t++) {
        const VolumeBase* volume = volumes->at(t);
        RealWorldMapping rwm = volume->getRealWorldMapping();
        tgt::mat4 worldToVoxelMatrix = volume->getWorldToVoxelMatrix();

        std::vector<PlotCellValue> values;
        if(timeSeries) {
            values.push_back(PlotCellValue(volume->getTimestep()));
        }
        else {
            values.push_back(PlotCellValue(t+1));
        }

        std::vector<tgt::vec3> samples;
        samples.reserve(seedPoints.size());

        // Gather samples.
        VolumeRAMRepresentationLock lock(volume);
        for(tgt::vec3 seedPoint : seedPoints) {
            seedPoint = worldToVoxelMatrix * seedPoint;

            tgt::vec3 sample = tgt::vec3::zero;
            for(size_t channel=0; channel < lock->getNumChannels(); channel++) {
                sample[channel] = rwm.normalizedToRealWorld(lock->getVoxelNormalizedLinear(seedPoint, channel));
            }
            samples.push_back(sample);
        }

        // Generate output.
        std::vector<float> output = outputFunc(samples);

        // Add to plot data.
        for(size_t i=0; i<output.size(); i++) {
            values.push_back(PlotCellValue(output[i]));
        }

        data->insert(values);

        progressReporter.setProgress((t+1.0f) / volumes->size());
    }

    // If we only have a single row, copy and shift to render the line visible.
    if(data->getRowsCount() == 1) {
        std::vector<PlotCellValue> copy = data->getRow(0).getCells();
        copy[0] = PlotCellValue(copy[0].getValue() + 1.0); // Add one second.
        data->insert(copy);
    }

    progressReporter.setProgress(1.0f);

    return RoiAnalysisOutput {
            std::move(data)
    };
}

void RoiAnalysis::processComputeOutput(RoiAnalysisOutput output) {
    outport_.setData(output.plotData.release(), true);
}

bool RoiAnalysis::isTimeSeries(const VolumeList* list) {
    std::set<float> timestamps;
    for(size_t t = 0; t < list->size(); t++) {
        if(!list->at(t)->hasMetaData(VolumeBase::META_DATA_NAME_TIMESTEP)) {
            return false;
        }
        else {
            timestamps.insert(list->at(t)->getTimestep());
        }
    }

    return timestamps.size() == list->size();
}

}
