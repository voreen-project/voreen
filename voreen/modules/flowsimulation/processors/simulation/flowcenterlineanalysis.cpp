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

#include "flowcenterlineanalysis.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

#include "../../utils/utils.h"

namespace voreen {

const std::string FlowCenterlineAnalysis::loggerCat_("voreen.flowsimulation.FlowCenterlineAnalysis");

FlowCenterlineAnalysis::FlowCenterlineAnalysis()
    : AsyncComputeProcessor()
    , ensembleDatasetPort_(Port::INPORT, "input.ensemble", "Ensemble Dataset Port")
    , outport_(Port::OUTPORT,"output.plot", "Plot Port")
    , vesselGraphFolder_("vesselGraphFolder", "Vessel Graph Folder", "Choose Path", "", "", FileDialogProperty::DIRECTORY)
    , outputQuantity_("outputQuantity", "Output Quantity")
    , numSamples_("numSamples", "Num. Samples", 10, 1, 100)
    , transformSamples_("transformSamples", "Transform Samples to Disk", false)
    , neighborSampleOverlap_("neighborSampleOverlap", "Neighbor Sample Overlap", 0.0f, 0.0f, 1.0f)
    , exportCurvePath_("exportCurvePath", "Export Curve Path", "Choose Path", "", "", FileDialogProperty::DIRECTORY, Processor::VALID)
    , saveButton_("saveButton", "Save", Processor::VALID)
{
    addPort(ensembleDatasetPort_);
    addPort(outport_);

    addProperty(vesselGraphFolder_);
    vesselGraphFolder_.setGroupID("input");
    setPropertyGroupGuiName("input", "Input");

    addProperty(outputQuantity_);
    outputQuantity_.addOption("maxMagnitude", "Max Magnitude");
    outputQuantity_.addOption("meanMagnitude", "Mean Magnitude");
    outputQuantity_.addOption("medianMagnitude", "Median Magnitude");
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
    //ON_CHANGE(exportCurvePath_, FlowCenterlineAnalysis, exportVelocityCurve);
    //addProperty(saveButton_);
    saveButton_.setGroupID("export");
    //ON_CHANGE(saveButton_, FlowCenterlineAnalysis, exportVelocityCurve);
    setPropertyGroupGuiName("export", "Export");
    //setPropertyGroupVisible("export", false); //TODO: can essentially be replaced by PlotDataExport
}

MatchingResult FlowCenterlineAnalysis::performMatching(std::unique_ptr<VesselGraph> vesselGraph) {

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

    tgt::plane nodePlane(vesselGraph->getNode(result.nodeMapping[MatchingResult::MPA]).pos_,
                         vesselGraph->getNode(n3).pos_,
                         vesselGraph->getNode(n4).pos_);

    if(tgt::dot(vesselGraph->getNode(bifurcationNode).pos_-vesselGraph->getNode(result.nodeMapping[MatchingResult::MPA]).pos_, nodePlane.n) > 0) {
        result.nodeMapping[MatchingResult::LPA] = n3;
        result.nodeMapping[MatchingResult::RPA] = n4;
    }
    else {
        result.nodeMapping[MatchingResult::RPA] = n3;
        result.nodeMapping[MatchingResult::LPA] = n4;
    }

    result.edgeMapping[MatchingResult::RPA] = vesselGraph->getNode(result.nodeMapping[MatchingResult::RPA]).getEdges().front().get().getID();
    result.edgeMapping[MatchingResult::LPA] = vesselGraph->getNode(result.nodeMapping[MatchingResult::LPA]).getEdges().front().get().getID();
    result.vesselGraph = std::move(vesselGraph);

    // Matching done.
    return result;
}

FlowCenterlineAnalysisInput FlowCenterlineAnalysis::prepareComputeInput() {

    if(vesselGraphFolder_.get().empty()) {
        throw InvalidInputException("Empty vessel graph path", InvalidInputException::S_ERROR);
    }

    if(exportCurvePath_.get().empty()) {
        throw InvalidInputException("Empty export path", InvalidInputException::S_ERROR);
    }

    auto ensemble = ensembleDatasetPort_.getThreadSafeData();
    if(!ensemble || ensemble->getMembers().empty()) {
        throw InvalidInputException("No ensemble", InvalidInputException::S_IGNORE);
    }

    std::map<std::string, MatchingResult> vesselGraphs;

    std::string vesselGraphPath = vesselGraphFolder_.get();
    std::vector<std::string> vesselGraphFiles = tgt::FileSystem::listFiles(vesselGraphPath, true);
    for(const std::string& file : vesselGraphFiles) {

        std::string path = vesselGraphPath + "/" + file;

        VesselGraphBuilder builder;
        auto output = std::move(builder).finalize();

        JsonDeserializer deserializer;
        try {
            std::fstream f(path, std::ios::in);
            bool compressed = tgt::FileSystem::fileExtension(path) == "gz";
            deserializer.read(f, compressed);
            deserializer.deserialize("graph", *output);

        } catch(SerializationException& e) {
            LERROR("Could not deserialize graph: " << e.what());

        } catch(...) {
            LERROR("Could not load xml file " << path);
        }

        // Check the vessel graph and perform the matching.
        std::string animal = tgt::FileSystem::baseName(file);
        try {
            vesselGraphs[animal] = performMatching(std::move(output));
        }
        catch(InvalidInputException& e) {
            LWARNING("Matching failed for " << file << ": " << e.msg_);
        }
    }


    std::string labels [] = {"MPA", "RPA", "LPA"};
    const size_t numSamples = numSamples_.get();

    std::function<std::vector<float>(const std::vector<tgt::vec3>&)> outputFunc;
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

    std::map<std::string, std::unique_ptr<PlotData>> output;
    for(const auto& member : ensemble->getMembers()) {
        std::unique_ptr<PlotData> plotData(new PlotData(1, MatchingResult::NUM * numSamples));
        plotData->setColumnLabel(0, "Time [s]");
        for (size_t i = 0; i < MatchingResult::NUM; i++) {
            for (size_t j = 0; j < numSamples; j++) {
                plotData->setColumnLabel(1 + numSamples * i + j, labels[i] + "_" + std::to_string(j));
            }
        }

        output[member.getName()] = std::move(plotData);
    }

    return FlowCenterlineAnalysisInput {
            std::move(ensemble),
            std::move(vesselGraphs),
            std::move(output),
            std::move(outputFunc),
            numSamples,
            transformSamples_.get(),
            neighborSampleOverlap_.get(),
            exportCurvePath_.get()
    };
}

FlowCenterlineAnalysisOutput FlowCenterlineAnalysis::compute(FlowCenterlineAnalysisInput input, ProgressReporter& progressReporter) const {

    auto ensemble = std::move(input.ensemble);
    auto vesselGraphs = std::move(input.vesselGraphs);
    std::map<std::string, std::unique_ptr<PlotData>> output = std::move(input.output);
    auto outputFunc = std::move(input.outputFunc);
    size_t numSamples = input.numSamples;
    float overlap = input.neighborSampleOverlap;

    progressReporter.setProgress(0.0f);

    size_t minNumTimeSteps = ensemble->getMinNumTimeSteps();

    for(const EnsembleMember& member : ensemble->getMembers()) {

        if(vesselGraphs.find(member.getName()) == vesselGraphs.end()) {
            continue;
        }

        const MatchingResult& matching = vesselGraphs[member.getName()];
        VesselGraph* vesselGraph = matching.vesselGraph.get();

        for (const auto& field: ensemble->getCommonFieldNames()) {

            PlotData data(*output[member.getName()].get());

            std::vector<TimeStep> timeSteps;
            for(size_t t=0; t<minNumTimeSteps; t++) {
                timeSteps.push_back(member.getTimeSteps()[t].createSubset({1, field}));
            }

            EnsembleDataset tmp;
            tmp.addMember(EnsembleMember(member.getName(), member.getColor(), timeSteps));
            auto volumes = tmp.getVolumes();

            bool timeSeries = isTimeSeries(volumes);
            if (timeSeries) {
                data.setColumnLabel(0, "Time [s]");
            } else {
                data.setColumnLabel(0, "Time Step");
            }

            for (size_t t = 0; t < volumes.size(); t++) {
                const VolumeBase* volume = volumes.at(t);
                std::vector<PlotCellValue> values;
                if (timeSeries) {
                    values.emplace_back(PlotCellValue(volume->getTimestep()));
                } else {
                    values.emplace_back(PlotCellValue(t + 1));
                }

                for (size_t i = 0; i < MatchingResult::NUM; i++) {

                    const auto& node = vesselGraph->getNode(matching.nodeMapping[i]);
                    const auto& edge = vesselGraph->getEdge(matching.edgeMapping[i]);

                    const auto& voxels = edge.getVoxels();
                    size_t numVoxels = voxels.size();

                    std::function<size_t(size_t)> index;
                    if (edge.getNode1().getID() == node.getID()) {
                        index = [](size_t i) { return i; };
                    } else {
                        index = [numVoxels](size_t i) { return numVoxels - 1 - i; };
                    }

                    //// Parametrize edge.
                    float length = edge.getLength();

                    // Calculate the length of a single sample.
                    float spacing = length / (numSamples + 1);
                    float lengthPerSample = (length / numSamples) * (1.0f + 2.0f * overlap);

                    for (size_t j = 0; j < numSamples; j++) {

                        size_t centerlinePosition = j * spacing;

                        size_t mid = std::min<size_t>(centerlinePosition, numVoxels - 1);
                        size_t num = 2; // Number of reference nodes in both directions.

                        size_t frontIdx = mid > num ? (mid - num) : 0;
                        size_t backIdx = std::min(mid + num, numVoxels - 1);

                        const VesselSkeletonVoxel* ref = &edge.getVoxels().at(index(mid));
                        const VesselSkeletonVoxel* front = &edge.getVoxels().at(index(frontIdx));
                        const VesselSkeletonVoxel* back = &edge.getVoxels().at(index(backIdx));

                        // Calculate average radius.
                        float radius = 0.0f;
                        for (size_t k = frontIdx; k <= backIdx; k++) {
                            radius += edge.getVoxels().at(index(k)).avgDistToSurface_;
                        }
                        radius /= (backIdx - frontIdx + 1);

                        tgt::vec3 center = ref->pos_;
                        tgt::vec3 normal = tgt::normalize(back->pos_ - front->pos_);

                        std::vector<tgt::vec3> samples = utils::sampleCylinder(volume, center, normal, radius,
                                                                               lengthPerSample,
                                                                               input.transformSamples);

                        // Generate output.
                        std::vector<float> output = outputFunc(samples);

                        // Add to plot data.
                        for (size_t k = 0; k < output.size(); k++) {
                            values.emplace_back(PlotCellValue(output[k]));
                        }
                    }
                }

                data.insert(values);

                progressReporter.setProgress((t + 1.0f) / volumes.size());
            }

            // If we only have a single row, copy and shift to render the line visible.
            if (data.getRowsCount() == 1) {
                std::vector<PlotCellValue> copy = data.getRow(0).getCells();
                copy[0] = PlotCellValue(copy[0].getValue() + 1.0); // Add one second.
                data.insert(copy);
            }

            std::string path = input.basePath + "/" + field + "/";
            tgt::FileSystem::createDirectoryRecursive(path);
            exportVelocityCurve(path + member.getName() + ".csv", &data);
        }
    }

    progressReporter.setProgress(1.0f);

    // TODO: only outputs first dataset.
    std::unique_ptr<PlotData> representative = std::move(output.begin()->second);
    output.erase(output.begin());
    return FlowCenterlineAnalysisOutput {
            std::move(representative)
    };
}

void FlowCenterlineAnalysis::processComputeOutput(FlowCenterlineAnalysisOutput output) {
    outport_.setData(output.plotData.release(), true);
}

void FlowCenterlineAnalysis::exportVelocityCurve(const std::string& path, PlotData* data) const {

    std::ofstream lineStream(path.c_str());
    if (lineStream.fail()) {
        LWARNING("CSV file could not be opened");
        return;
    }

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

bool FlowCenterlineAnalysis::isTimeSeries(const std::vector<const VolumeBase*>& list) {
    std::set<float> timestamps;
    for(size_t t = 0; t < list.size(); t++) {
        if(!list[t]->hasMetaData(VolumeBase::META_DATA_NAME_TIMESTEP)) {
            return false;
        }
        else {
            timestamps.insert(list[t]->getTimestep());
        }
    }

    return timestamps.size() == list.size();
}

}
