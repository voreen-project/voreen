/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2016 University of Muenster, Germany.                        *
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

#include "vesselgraphglobalstats.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "custommodules/bigdataimageprocessing/util/csvwriter.h"

namespace voreen {

const std::string VesselGraphGlobalStats::loggerCat_("voreen.vesselgraphstatplotter");


#define ADD_EDGE_PROPERTY(Type, name) \
    activeEdgeProperty_.addOption("edge" #name, "Edge " #name, EdgeProperty::Type); \

VesselGraphGlobalStats::VesselGraphGlobalStats()
    : Processor()
    , graphInport_(Port::INPORT, "graph.input", "Graph Input")
    , segmentExportFilePath_("segmentExportFilePath", "Segment Export File Path", "Export File Path", "", "*.csv", FileDialogProperty::SAVE_FILE)
    , nodeExportFilePath_("nodeExportFilePath", "Node Export File Path", "Export File Path", "", "*.csv", FileDialogProperty::SAVE_FILE)
    , numNodes_("numNodes", "Number of Nodes", 0, 0, std::numeric_limits<int>::max())
    , numEdges_("numEdges", "Number of Edges", 0, 0, std::numeric_limits<int>::max())
    , exportButton_("exportButton_", "Export")
    , autoExport_("autoExport", "Auto Export", false)
    , exportForced_(false)
{
    addPort(graphInport_);

    addProperty(segmentExportFilePath_);
    addProperty(nodeExportFilePath_);
    addProperty(autoExport_);
    addProperty(exportButton_);
        ON_CHANGE_LAMBDA(exportButton_, [this] () {
                exportForced_ = true;
                });

    addProperty(numNodes_);
        numNodes_.setReadOnlyFlag(true);
    addProperty(numEdges_);
        numEdges_.setReadOnlyFlag(true);
}

VesselGraphGlobalStats::~VesselGraphGlobalStats() {
}

void VesselGraphGlobalStats::process() {
    tgtAssert(graphInport_.isReady(), "inport not ready");
    if(graphInport_.hasChanged()) {
        adaptToNewInput();
    }
    const VesselGraph* graph = graphInport_.getData();
    if(!graph) {
        return;
    }


    numNodes_.set(graph->getNodes().size());
    numEdges_.set(graph->getEdges().size());

    // Write to file
    if(exportForced_ || autoExport_.get()) {
        exportForced_ = false;
        exportToFile(*graph);
    }
}
void VesselGraphGlobalStats::exportToFile(const VesselGraph& graph) {
    if(segmentExportFilePath_.get().empty()) {
        LERROR("No segment export file path");
        return;
    } else {
        try {
            CSVWriter<size_t, size_t, size_t, float, float, float, float, float, float, float, float, float, float, float, float, float, size_t, size_t, size_t, bool> writer(segmentExportFilePath_.get());
            writer.writeHeader("id", "node1id", "node2id", "length", "distance", "curveness", "volume", "avgCrossSection", "minRadiusAvg", "minRadiusStd", "avgRadiusAvg", "avgRadiusStd", "maxRadiusAvg", "maxRadiusStd", "roundnessAvg", "roundnessStd", "node1_degree", "node2_degree", "num_voxels", "hasNodeAtSampleBorder");
            for(const VesselGraphEdge& edge : graph.getEdges()) {
                writer.write(
                        edge.getID(),
                        edge.getNodeID1(),
                        edge.getNodeID2(),
                        edge.getLength(),
                        edge.getDistance(),
                        edge.getCurveness(),
                        edge.getVolume(),
                        edge.getAvgCrossSection(),
                        edge.getMinRadiusAvg(),
                        edge.getMinRadiusStdDeviation(),
                        edge.getAvgRadiusAvg(),
                        edge.getAvgRadiusStdDeviation(),
                        edge.getMaxRadiusAvg(),
                        edge.getMaxRadiusStdDeviation(),
                        edge.getRoundnessAvg(),
                        edge.getRoundnessStdDeviation(),
                        std::min(edge.getNode1().getDegree(), edge.getNode2().getDegree()),
                        std::max(edge.getNode1().getDegree(), edge.getNode2().getDegree()),
                        edge.getVoxels().size(),
                        edge.getNode1().isAtSampleBorder_ || edge.getNode2().isAtSampleBorder_
                        );
            }
            LINFO("Writing edge stats: " << segmentExportFilePath_.get());
        } catch(tgt::IOException& e) {
            LERROR("Error opening file " << segmentExportFilePath_.get() << ": " << e.what());
        }
    }

    if(nodeExportFilePath_.get().empty()) {
        LERROR("No node export file path");
        return;
    } else {
        try {
            CSVWriter<size_t, float, float, float, int, bool> writer(nodeExportFilePath_.get());
            writer.writeHeader("id", "pos_x", "pos_y", "pos_z", "degree", "isAtSampleBorder");
            for(const VesselGraphNode& node : graph.getNodes()) {
                writer.write(
                        node.getID(),
                        node.pos_.x,
                        node.pos_.y,
                        node.pos_.z,
                        node.getDegree(),
                        node.isAtSampleBorder_
                        );
            }
            LINFO("Writing node stats: " << nodeExportFilePath_.get());
        } catch(tgt::IOException& e) {
            LERROR("Error opening file " << nodeExportFilePath_.get() << ": " << e.what());
        }
    }
}

void VesselGraphGlobalStats::initialize() {
    Processor::initialize();
}

void VesselGraphGlobalStats::deinitialize() {
    Processor::deinitialize();
}
bool VesselGraphGlobalStats::isEndProcessor() const {
    return true;
}
void VesselGraphGlobalStats::adaptToNewInput() {
}
} // namespace voreen
