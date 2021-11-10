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

#include "flowcenterlineanalysis.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

#include "../../utils/utils.h"

namespace voreen {

const std::string FlowCenterlineAnalysis::loggerCat_("voreen.flowsimulation.FlowCenterlineAnalysis");

FlowCenterlineAnalysis::FlowCenterlineAnalysis()
    : AsyncComputeProcessor()
    , volumeListPort_(Port::INPORT, "input.volumelist", "Volume List Port")
    , vesselGraphPort_(Port::INPORT, "input.vesselgraph", "Vessel Graph Port")
    , outport_(Port::OUTPORT,"output.plot", "Plot Port")
    , outputQuantity_("outputQuantity", "Output Quantity")
    , numSamples_("numSamples", "Num. Samples", 10, 1, 100)
    , transformSamples_("transformSamples", "Transform Samples to Disk", false)
    , neighborSampleOverlap_("neighborSampleOverlap", "Neighbor Sample Overlap", 0.0f, 0.0f, 1.0f)
    , exportCurvePath_("exportCurvePath", "Export Curve Path", "Choose Path", "", ".csv", FileDialogProperty::SAVE_FILE, Processor::VALID)
    , saveButton_("saveButton", "Save", Processor::VALID)
{
    addPort(volumeListPort_);
    addPort(vesselGraphPort_);
    volumeListPort_.addCondition(new PortConditionVolumeListEnsemble());
    addPort(outport_);

    addProperty(outputQuantity_);
    outputQuantity_.addOption("maxMagnitude", "Max Magnitude");
    outputQuantity_.addOption("meanMagnitude", "Mean Magnitude");
    outputQuantity_.addOption("medianMagnitude", "Median Magnitude");
    outputQuantity_.addOption("meanComponents", "Mean Components");
    outputQuantity_.setGroupID("output");
    outputQuantity_.invalidate();

    addProperty(numSamples_);
    numSamples_.setGroupID("output");

    addProperty(transformSamples_);
    transformSamples_.setGroupID("output");

    addProperty(neighborSampleOverlap_);
    neighborSampleOverlap_.setNumDecimals(3);
    neighborSampleOverlap_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

    addProperty(exportCurvePath_);
    exportCurvePath_.setGroupID("export");
    ON_CHANGE(exportCurvePath_, FlowCenterlineAnalysis, exportVelocityCurve);
    addProperty(saveButton_);
    saveButton_.setGroupID("export");
    ON_CHANGE(saveButton_, FlowCenterlineAnalysis, exportVelocityCurve);
    setPropertyGroupGuiName("export", "Export");
    //setPropertyGroupVisible("export", false); //TODO: can essentially be replaced by PlotDataExport
}

bool FlowCenterlineAnalysis::isReady() const {
    if(!volumeListPort_.isReady()) {
        setNotReadyErrorMessage("No volume list connected");
        return false;
    }

    if(!vesselGraphPort_.isReady()) {
        setNotReadyErrorMessage("No vessel graph connected");
        return false;
    }

    // Note: parametrization is optional

    return true;
}

MatchingResult FlowCenterlineAnalysis::performMatching() {
    if(!vesselGraphPort_.hasData()) {
        throw InvalidInputException("No vessel graph", InvalidInputException::S_ERROR);
    }

    auto vesselGraph = vesselGraphPort_.getData();

    size_t numNodes = vesselGraph->getNodes().size();
    if(numNodes != 4) {
        throw InvalidInputException("Num nodes has to be 4", InvalidInputException::S_ERROR);
    }
    size_t numEdges = vesselGraph->getEdges().size();
    if(numEdges != 3) {
        throw InvalidInputException("Num edges has to be 3", InvalidInputException::S_ERROR);
    }

    std::vector<VGNodeID> onering;
    VGNodeID bifurcationNode;
    for(auto& node : vesselGraph->getNodes()) {
        if(node.getDegree() == 1) {
            onering.push_back(node.getID());
        }
        else if(node.getDegree() == 3) {
            if(bifurcationNode != -1) {
                throw InvalidInputException("Vessel Graph degenerated", InvalidInputException::S_ERROR);
            }
            bifurcationNode = node.getID();
        }
        else {
            throw InvalidInputException("Vessel Graph degenerated", InvalidInputException::S_ERROR);
        }
    }

    tgtAssert(onering.size() == 3, "Invalid one ring size");

    // Perform the mapping.
    // The mapping is performed as follows: (Note that we can only map the pulmonary artery).
    // 1. Create the plane all the nodes reside in.
    // 2. Find MPA (largest avg. radius)
    // 3. Test on which side of the plane the bifurcation node resides.
    //

    MatchingResult result;

    tgt::vec3 p0 = vesselGraph->getNode(onering[0]).pos_;
    tgt::vec3 p1 = vesselGraph->getNode(onering[1]).pos_;
    tgt::vec3 p2 = vesselGraph->getNode(onering[2]).pos_;

    auto e1 = vesselGraph->getNode(onering[0]).getEdges().front();
    auto e2 = vesselGraph->getNode(onering[1]).getEdges().front();
    auto e3 = vesselGraph->getNode(onering[2]).getEdges().front();

    // Find MPA edge.
    if(e1.get().getAvgRadiusAvg() > e2.get().getAvgRadiusAvg()) {
        if(e1.get().getAvgRadiusAvg() > e3.get().getAvgRadiusAvg()) {
            result.edgeMapping[MatchingResult::MPA] = e1.get().getID();
        }
        else {
            result.edgeMapping[MatchingResult::MPA] = e3.get().getID();
        }
    }
    else {
        if(e2.get().getAvgRadiusAvg() > e3.get().getAvgRadiusAvg()) {
            result.edgeMapping[MatchingResult::MPA] = e2.get().getID();
        }
        else {
            result.edgeMapping[MatchingResult::MPA] = e3.get().getID();
        }
    }

    // Find MPA Node.
    if(vesselGraph->getEdge(result.edgeMapping[MatchingResult::MPA]).getNode1().getID() == bifurcationNode) {
        result.nodeMapping[MatchingResult::MPA] = vesselGraph->getEdge(result.edgeMapping[MatchingResult::MPA]).getNodeID2();
    }
    else {
        result.nodeMapping[MatchingResult::MPA] = vesselGraph->getEdge(result.edgeMapping[MatchingResult::MPA]).getNodeID1();
    }


    // Both mpa nodes.
    const auto& n1 = vesselGraph->getNode(result.nodeMapping[MatchingResult::MPA]);
    const auto& n2 = vesselGraph->getNode(bifurcationNode);
    VGNodeID n3, n4;
    if(onering[0] == result.nodeMapping[MatchingResult::MPA]) {
        n3 = onering[1];
        n4 = onering[2];
    }
    else if(onering[1] == result.nodeMapping[MatchingResult::MPA]) {
        n3 = onering[2];
        n4 = onering[0];
    }
    else if(onering[2] == result.nodeMapping[MatchingResult::MPA]) {
        n3 = onering[0];
        n4 = onering[1];
    }

    tgt::vec3 pos = vesselGraph->getNode(n3).pos_ + 0.5f * (vesselGraph->getNode(n4).pos_ - vesselGraph->getNode(n3).pos_);

    tgt::plane mpaPlane(n1.pos_, n2.pos_, pos);
    tgt::plane nodePlane(p0, p1, p2);

    if(nodePlane.distance(vesselGraph->getNode(bifurcationNode).pos_) > 0) {
        if(mpaPlane.distance(vesselGraph->getNode(n3).pos_) > 0) {
            result.nodeMapping[MatchingResult::RPA] = n3;
            result.nodeMapping[MatchingResult::LPA] = n4;
        }
        else {
            result.nodeMapping[MatchingResult::RPA] = n4;
            result.nodeMapping[MatchingResult::LPA] = n3;
        }
    }
    else {
        if(mpaPlane.distance(vesselGraph->getNode(n4).pos_) > 0) {
            result.nodeMapping[MatchingResult::RPA] = n4;
            result.nodeMapping[MatchingResult::LPA] = n3;
        }
        else {
            result.nodeMapping[MatchingResult::RPA] = n3;
            result.nodeMapping[MatchingResult::LPA] = n4;
        }
    }

    result.edgeMapping[MatchingResult::RPA] = vesselGraph->getNode(result.nodeMapping[MatchingResult::RPA]).getEdges().front().get().getID();
    result.edgeMapping[MatchingResult::LPA] = vesselGraph->getNode(result.nodeMapping[MatchingResult::LPA]).getEdges().front().get().getID();

    // Matching done.
    return result;
}

FlowCenterlineAnalysisInput FlowCenterlineAnalysis::prepareComputeInput() {

    auto volumes = volumeListPort_.getThreadSafeData();
    if(!volumes || volumes->empty()) {
        throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
    }

    const VolumeBase* reference = volumes->first();
    const size_t numChannels = reference->getNumChannels();
    if(numChannels != 1 && numChannels != 3) {
        throw InvalidInputException("Only 1 and 3 channel volumes supported", InvalidInputException::S_ERROR);
    }

    // Check the vessel graph and perform the matching.
    MatchingResult matchingResult = performMatching();

    // Vessel graph is safe to use now.
    auto vesselGraph = vesselGraphPort_.getThreadSafeData();

    std::string labels [] = {"MPA", "RPA", "LPA"};
    const size_t numSamples = numSamples_.get();

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

        output.reset(new PlotData(1, MatchingResult::NUM*numSamples*numChannels));
        output->setColumnLabel(0, "Time [s]");
        if(numChannels == 1) {
            for(size_t i=0; i<MatchingResult::NUM; i++) {
                for(size_t j=0; j<numSamples; j++) {
                    output->setColumnLabel(1 + numSamples * i + j, labels[i] + "_" + std::to_string(j));
                }
            }
        }
        else if(numChannels == 3) {
            if(transformSamples_.get()) {
                for(size_t i=0; i<MatchingResult::NUM; i++) {
                    for(size_t j=0; j<numSamples; j++) {
                        output->setColumnLabel(1 + 3 * (numSamples * i + j) + 0, "InP 1"); // in plane
                        output->setColumnLabel(1 + 3 * (numSamples * i + j) + 1, "InP 2"); // in plane
                        output->setColumnLabel(1 + 3 * (numSamples * i + j) + 2, "ThP"); // through plane
                    }
                }
            }
            else {
                for(size_t i=0; i<MatchingResult::NUM; i++) {
                    for(size_t j=0; j<numSamples; j++) {
                        output->setColumnLabel(1 + 3 * (numSamples * i + j) + 0, "x");
                        output->setColumnLabel(1 + 3 * (numSamples * i + j) + 1, "y");
                        output->setColumnLabel(1 + 3 * (numSamples * i + j) + 2, "z");
                    }
                }
            }
        }
    }
    else {

        output.reset(new PlotData(1, MatchingResult::NUM*numSamples));
        output->setColumnLabel(0, "Time [s]");
        for(size_t i=0; i<MatchingResult::NUM; i++) {
            for(size_t j=0; j<numSamples; j++) {
                output->setColumnLabel(1 + numSamples * i + j, labels[i] + "_" + std::to_string(j));
            }
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

    return FlowCenterlineAnalysisInput {
            std::move(volumes),
            std::move(vesselGraph),
            matchingResult,
            std::move(output),
            std::move(outputFunc),
            numSamples,
            transformSamples_.get(),
            neighborSampleOverlap_.get()
    };
}

FlowCenterlineAnalysisOutput FlowCenterlineAnalysis::compute(FlowCenterlineAnalysisInput input, ProgressReporter& progressReporter) const {

    auto volumes = std::move(input.volumes);
    auto vesselGraph = std::move(input.vesselGraph);
    auto matching = input.matchingResult;
    std::unique_ptr<PlotData> data = std::move(input.output);
    auto outputFunc = std::move(input.outputFunc);
    size_t numSamples = input.numSamples;
    float overlap = input.neighborSampleOverlap;

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

        for(size_t i=0; i<MatchingResult::NUM; i++) {

            const auto& node = vesselGraph->getNode(matching.nodeMapping[i]);
            const auto& edge = vesselGraph->getEdge(matching.edgeMapping[i]);

            const auto& voxels = edge.getVoxels();
            size_t numVoxels = voxels.size();

            std::function<size_t (size_t)> index;
            if(edge.getNode1().getID() == node.getID()) {
                index = [](size_t i) { return i; };
            }
            else {
                index = [numVoxels](size_t i) { return numVoxels - 1 - i;  };
            }

            //// Parametrize edge.
            float length = edge.getLength();

            // Calculate the length of a single sample.
            float spacing = length / (numSamples+1);
            float lengthPerSample = (length / numSamples) * (1.0f + 2.0f * overlap);

            for(size_t j=0; j<numSamples; j++) {

                size_t centerlinePosition = j * spacing;

                size_t mid = std::min<size_t>(centerlinePosition, numVoxels - 1);
                size_t num = 2; // Number of reference nodes in both directions.

                size_t frontIdx = mid > num ? (mid - num) : 0;
                size_t backIdx = std::min(mid + num, numVoxels - 1);

                const VesselSkeletonVoxel* ref   = &edge.getVoxels().at(index(mid));
                const VesselSkeletonVoxel* front = &edge.getVoxels().at(index(frontIdx));
                const VesselSkeletonVoxel* back  = &edge.getVoxels().at(index(backIdx));

                // Calculate average radius.
                float radius = 0.0f;
                for (size_t k = frontIdx; k <= backIdx; k++) {
                    radius += edge.getVoxels().at(index(k)).avgDistToSurface_;
                }
                radius /= (backIdx - frontIdx + 1);

                tgt::vec3 center = ref->pos_;
                tgt::vec3 normal = tgt::normalize(back->pos_ - front->pos_);

                std::vector<tgt::vec3> samples = utils::sampleCylinder(volume, center, normal, radius, lengthPerSample,
                                                                       input.transformSamples);

                // Generate output.
                std::vector<float> output = outputFunc(samples);

                // Add to plot data.
                for (size_t k = 0; k < output.size(); k++) {
                    values.emplace_back(PlotCellValue(output[k]));
                }
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

    return FlowCenterlineAnalysisOutput {
            std::move(data)
    };
}

void FlowCenterlineAnalysis::processComputeOutput(FlowCenterlineAnalysisOutput output) {
    outport_.setData(output.plotData.release(), true);
}

void FlowCenterlineAnalysis::exportVelocityCurve() {

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

bool FlowCenterlineAnalysis::isTimeSeries(const VolumeList* list) {
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
