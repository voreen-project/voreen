/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "flowindicatoranalysis.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

#include "../../utils/utils.h"

namespace voreen {

const std::string FlowIndicatorAnalysis::loggerCat_("voreen.flowsimulation.FlowIndicatorAnalysis");

FlowIndicatorAnalysis::FlowIndicatorAnalysis()
    : AsyncComputeProcessor()
    , volumeListPort_(Port::INPORT, "input.volumelist", "Volume List Port")
    , parameterPort_(Port::INPORT, "input.parameter", "Parameter Port")
    , maskPort_(Port::INPORT, "input.mask", "Sample Mask Port (Optional)")
    , outport_(Port::OUTPORT,"output.plot", "Plot Port")
    , outputQuantity_("outputQuantity", "Output Quantity")
    , indicator_("indicator", "Indicator", Processor::INVALID_RESULT, true)
    , transformSamples_("transformSamples", "Transform Samples to Disk", false)
    , exportCurvePath_("exportCurvePath", "Export Curve Path", "Choose Path", "", ".csv", FileDialogProperty::SAVE_FILE, Processor::VALID)
    , saveButton_("saveButton", "Save", Processor::VALID)
{
    addPort(volumeListPort_);
    addPort(parameterPort_);
    ON_CHANGE(parameterPort_, FlowIndicatorAnalysis, onParametersChange);
    volumeListPort_.addCondition(new PortConditionVolumeListEnsemble());
    addPort(maskPort_);
    addPort(outport_);

    addProperty(outputQuantity_);
    outputQuantity_.addOption("maxMagnitude", "Max Magnitude");
    outputQuantity_.addOption("meanMagnitude", "Mean Magnitude");
    outputQuantity_.addOption("medianMagnitude", "Median Magnitude");
    outputQuantity_.addOption("meanComponents", "Mean Components");
    outputQuantity_.setGroupID("output");

    addProperty(indicator_);
    indicator_.setGroupID("output");

    addProperty(transformSamples_);
    transformSamples_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

    addProperty(exportCurvePath_);
    exportCurvePath_.setGroupID("export");
    ON_CHANGE(exportCurvePath_, FlowIndicatorAnalysis, exportVelocityCurve);
    addProperty(saveButton_);
    saveButton_.setGroupID("export");
    ON_CHANGE(saveButton_, FlowIndicatorAnalysis, exportVelocityCurve);
    setPropertyGroupGuiName("export", "Export");
    //setPropertyGroupVisible("export", false); //TODO: can essentially be replaced by PlotDataExport
}

bool FlowIndicatorAnalysis::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Processor not initialized");
        return false;
    }

    auto* volumeList = volumeListPort_.getData();
    if(!volumeList || volumeList->empty()) {
        setNotReadyErrorMessage("No volume list");
        return false;
    }

    auto* config = parameterPort_.getData();
    if(!config) {
        setNotReadyErrorMessage("No parameter config");
        return false;
    }

    // Note: Sample mask is optional.

    return true;
}

void FlowIndicatorAnalysis::onParametersChange() {
    if(!parameterPort_.hasData()) {
        return;
    }

    std::string selected = indicator_.get();
    indicator_.setOptions(std::deque<Option<std::string>>());

    // Add "all" as extra option.
    // TODO: make sure that no indicator has been named "all".
    indicator_.addOption("all", "all");

    for(const FlowIndicator& indicator : parameterPort_.getData()->getFlowIndicators()) {
        std::string id = indicator.name_;//std::to_string(indicator.id_); // Name seems more intuitive than id!
        indicator_.addOption(id, indicator.name_);

        // Select old selected entry, if it is still available.
        if(id == selected) {
            indicator_.select(selected);
        }
    }
}

FlowIndicatorAnalysisInput FlowIndicatorAnalysis::prepareComputeInput() {

    auto volumes = volumeListPort_.getThreadSafeData();
    if(!volumes || volumes->empty()) {
        throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
    }
    auto parameterSetEnsemble = parameterPort_.getThreadSafeData();
    if(!parameterSetEnsemble || parameterSetEnsemble->getFlowIndicators().empty()) {
        throw InvalidInputException("No indicators", InvalidInputException::S_IGNORE);
    }

    const VolumeBase* reference = volumes->first();
    const size_t numChannels = reference->getNumChannels();
    if(numChannels != 1 && numChannels != 3) {
        throw InvalidInputException("Only 1 and 3 channel volumes supported", InvalidInputException::S_ERROR);
    }

    std::vector<FlowIndicator> indicators;
    for(const FlowIndicator& indicator : parameterSetEnsemble->getFlowIndicators()) {
        if(indicator.name_ == indicator_.getValue() || indicator_.getValue() == "all") {
            indicators.push_back(indicator);
        }
    }

    std::function<std::vector<float>(const std::vector<tgt::vec3>&)> outputFunc;
    std::unique_ptr<PlotData> output;

    if(outputQuantity_.get() == "meanComponents") {
        outputFunc = [] (const std::vector<tgt::vec3>& samples) {
            tgt::vec3 mean = tgt::vec3::zero;
            for(const auto& sample : samples) {
                mean += sample;
            }

            mean /= static_cast<float>(samples.size());

            std::vector<float> output(2);
            output[0] = std::sqrt(mean.x*mean.x + mean.y*mean.y);
            output[1] = mean.z;

            return output;
        };

        auto numIndicators = indicators.size();
        tgtAssert(numIndicators > 0, "invalid state");

        int numComponents = numChannels;
        if (numChannels == 3 && transformSamples_.get()) {
            numComponents = 2;
        }

        output.reset(new PlotData(1, numComponents * numIndicators));
        for (size_t i=0; i<indicators.size(); i++) {
            auto name = indicators[i].name_;
            if(numChannels == 1) {
                output->setColumnLabel(i + 1, name);
            }
            else if(numChannels == 3) {
                auto suffix = " - " + name;
                if(transformSamples_.get()) {
                    output->setColumnLabel(i * 2 + 1, "InP" + suffix); // in plane
                    output->setColumnLabel(i * 2 + 2, "ThP" + suffix); // through plane
                }
                else {
                    output->setColumnLabel(i * 3 + 1, "x" + suffix);
                    output->setColumnLabel(i * 3 + 2, "y" + suffix);
                    output->setColumnLabel(i * 3 + 3, "z" + suffix);
                }
            }
        }
    }
    else {

        output.reset(new PlotData(1, indicators.size()));
        output->setColumnLabel(0, "Time [s]");
        for(size_t i=0; i<indicators.size(); i++) {
            output->setColumnLabel(i+1, indicators[i].name_);
        }

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
                    meanMagnitude += tgt::length(sample) / static_cast<float>(samples.size());
                }
                return std::vector<float>(1, meanMagnitude);
            };
        }
        else if(outputQuantity_.get() == "medianMagnitude") {
            outputFunc = [](const std::vector<tgt::vec3>& samples) {
                std::vector<float> magnitudes;
                magnitudes.reserve(samples.size());
                for (const auto& sample : samples) {
                    magnitudes.emplace_back(tgt::length(sample));
                }
                std::nth_element(magnitudes.begin(), magnitudes.begin() + magnitudes.size()/2, magnitudes.end());
                return std::vector<float>(1, magnitudes[magnitudes.size() / 2]);
            };
        }
        else {
            tgtAssert(false, "unhandled output quantity");
        }
    }

    auto sampleMask = maskPort_.getThreadSafeData();

    return FlowIndicatorAnalysisInput {
        std::move(volumes),
        std::move(sampleMask),
        std::move(indicators),
        std::move(output),
        std::move(outputFunc),
        transformSamples_.get()
    };
}

FlowIndicatorAnalysisOutput FlowIndicatorAnalysis::compute(FlowIndicatorAnalysisInput input, ProgressReporter& progressReporter) const {

    auto volumes = std::move(input.volumes);
    auto sampleMask = std::move(input.sampleMask);
    auto indicators = std::move(input.indicators);
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
        std::vector<PlotCellValue> values;
        if(timeSeries) {
            values.emplace_back(PlotCellValue(volume->getTimestep()));
        }
        else {
            values.emplace_back(PlotCellValue(t+1));
        }
        for(const FlowIndicator& indicator : indicators) {
            // Gather samples.
            std::vector<tgt::vec3> samples = utils::sampleCylinder(volume, indicator.center_, indicator.normal_, indicator.radius_, input.transformSamples, 0, sampleMask);

            // Generate output.
            std::vector<float> output = outputFunc(samples);

            // Add to plot data.
            for(size_t i=0; i<output.size(); i++) {
                values.emplace_back(PlotCellValue(output[i]));
            }
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

    return FlowIndicatorAnalysisOutput {
        std::move(data)
    };
}

void FlowIndicatorAnalysis::processComputeOutput(FlowIndicatorAnalysisOutput output) {
    outport_.setData(output.plotData.release(), true);
}

void FlowIndicatorAnalysis::exportVelocityCurve() {

    if(!outport_.hasData()) {
        LWARNING("No data available");
        return;
    }

    const std::string& file = exportCurvePath_.get();
    if(file.empty()) {
        LWARNING("Empty path");
        return;
    }

    std::ofstream lineStream(file.c_str());
    if (lineStream.fail()) {
        LWARNING("CSV file could not be opened");
        return;
    }

    const PlotData* data = dynamic_cast<const PlotData*>(outport_.getData());
    tgtAssert(data, "Must be of type PlotData");

    lineStream << data->getColumnLabel(0);
    for(int i=1; i<data->getColumnCount(); i++) {
        lineStream << "," << data->getColumnLabel(i);
    }
    lineStream << std::endl;

    for(int i=0; i<data->getRowsCount(); i++) {
        auto cells = data->getRow(i).getCells();
        auto back = cells.back();
        cells.pop_back(); // Handle back separately.
        for(const auto& cell : cells) {
            lineStream << cell.getValue() << ",";
        }
        lineStream << back.getValue() << std::endl;
    }
}

bool FlowIndicatorAnalysis::isTimeSeries(const VolumeList* list) {
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
