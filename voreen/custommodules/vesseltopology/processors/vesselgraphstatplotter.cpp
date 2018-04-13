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

#include "vesselgraphstatplotter.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "custommodules/bigdataimageprocessing/util/csvwriter.h"

#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"
#include "modules/plotting/datastructures/plotcell.h"

namespace voreen {

const std::string VesselGraphStatPlotter::loggerCat_("voreen.vesselgraphstatplotter");

VesselGraphStatPlotter::VesselGraphStatPlotter()
    : Processor()
    , graphInport_(Port::INPORT, "graph.input", "Graph Input")
    , plotOutport_(Port::OUTPORT, "plot.output", "Plot Data Output", false,  Processor::VALID)
    , exportFilePath_("exportFilePath", "Export File Path", "Export File Path", "", "*.csv", FileDialogProperty::SAVE_FILE)
    , exportButton_("exportButton", "Export")
    , activeEdgeID_("activeEdgeID", "Active Edge ID", -1, -1, std::numeric_limits<int>::max()-1)
    , length_("length", "Length", 0, 0, std::numeric_limits<float>::max())
    , distance_("distance", "Distance", 0, 0, std::numeric_limits<float>::max())
    , curveness_("curveness", "Curveness", 0, -1, std::numeric_limits<float>::max())
    , minRadiusAvg_("minRadiusAvg", "Min Radius Mean", 0, 0, std::numeric_limits<float>::max())
    , minRadiusStdDeviation_("minRadiusStdDeviation", "Min Radius Std Deviation", 0, 0, std::numeric_limits<float>::max())
    , avgRadiusAvg_("avgRadiusAvg", "Avg Radius Mean", 0, 0, std::numeric_limits<float>::max())
    , avgRadiusStdDeviation_("avgRadiusStdDeviation", "Avg Radius Std Deviation", 0, 0, std::numeric_limits<float>::max())
    , maxRadiusAvg_("maxRadiusAvg", "Max Radius Mean", 0, 0, std::numeric_limits<float>::max())
    , maxRadiusStdDeviation_("maxRadiusStdDeviation", "Max Radius Std Deviation", 0, 0, std::numeric_limits<float>::max())
    , roundnessAvg_("roundnessAvg", "Roundness Mean", 0, 0, std::numeric_limits<float>::max())
    , roundnessStdDeviation_("roundnessStdDeviation", "Roundness Std Deviation", 0, 0, std::numeric_limits<float>::max())
    , avgCrossSection_("avgCrossSection", "Average Cross Section", 0, 0, std::numeric_limits<float>::max())
    , volume_("volume", "Volume", 0, 0, std::numeric_limits<float>::max())
    , numEdgesLeftNode_("numEdgesLeftNode", "Num Edges Left Node", 0, 0, std::numeric_limits<int>::max()-1)
    , numEdgesRightNode_("numEdgesRightNode", "Num Edges Right Node", 0, 0, std::numeric_limits<int>::max()-1)
    , numSkeletonVoxels_("numSkeletonVoxels", "Num Skeleton Voxels", 0, 0, std::numeric_limits<int>::max()-1)
    , exportForced_(false)
{
    addPort(graphInport_);
    addPort(plotOutport_);
    addProperty(activeEdgeID_);

    addProperty(exportFilePath_);
    addProperty(exportButton_);
        ON_CHANGE_LAMBDA(exportButton_, [this] () {
                exportForced_ = true;
                });

    const int decimals = 5;
    addProperty(length_);
        length_.setReadOnlyFlag(true);
        length_.setNumDecimals(decimals);
    addProperty(distance_);
        distance_.setReadOnlyFlag(true);
        distance_.setNumDecimals(decimals);
    addProperty(curveness_);
        curveness_.setReadOnlyFlag(true);
        curveness_.setNumDecimals(decimals);
    addProperty(minRadiusAvg_);
        minRadiusAvg_.setReadOnlyFlag(true);
        minRadiusAvg_.setNumDecimals(decimals);
    addProperty(minRadiusStdDeviation_);
        minRadiusStdDeviation_.setReadOnlyFlag(true);
        minRadiusStdDeviation_.setNumDecimals(decimals);
    addProperty(avgRadiusAvg_);
        avgRadiusAvg_.setReadOnlyFlag(true);
        avgRadiusAvg_.setNumDecimals(decimals);
    addProperty(avgRadiusStdDeviation_);
        avgRadiusStdDeviation_.setReadOnlyFlag(true);
        avgRadiusStdDeviation_.setNumDecimals(decimals);
    addProperty(maxRadiusAvg_);
        maxRadiusAvg_.setReadOnlyFlag(true);
        maxRadiusAvg_.setNumDecimals(decimals);
    addProperty(maxRadiusStdDeviation_);
        maxRadiusStdDeviation_.setReadOnlyFlag(true);
        maxRadiusStdDeviation_.setNumDecimals(decimals);
    addProperty(roundnessAvg_);
        roundnessAvg_.setReadOnlyFlag(true);
        roundnessAvg_.setNumDecimals(decimals);
    addProperty(roundnessStdDeviation_);
        roundnessStdDeviation_.setReadOnlyFlag(true);
        roundnessStdDeviation_.setNumDecimals(decimals);
    addProperty(avgCrossSection_);
        avgCrossSection_.setReadOnlyFlag(true);
        avgCrossSection_.setNumDecimals(decimals);
    addProperty(volume_);
        volume_.setReadOnlyFlag(true);
        volume_.setNumDecimals(decimals);
    addProperty(numEdgesLeftNode_);
        numEdgesLeftNode_.setReadOnlyFlag(true);
    addProperty(numEdgesRightNode_);
        numEdgesRightNode_.setReadOnlyFlag(true);
    addProperty(numSkeletonVoxels_);
        numSkeletonVoxels_.setReadOnlyFlag(true);

}

VesselGraphStatPlotter::~VesselGraphStatPlotter() {
}

void VesselGraphStatPlotter::process() {
    tgtAssert(graphInport_.isReady(), "inport not ready");
    if(graphInport_.hasChanged()) {
        adaptToNewInput();
    }
    const VesselGraph* graph = graphInport_.getData();
    if(!graph) {
        return;
    }

    if(activeEdgeID_.get() == -1) {
        return;
    }


    const VesselGraphEdge& edge = graph->getEdges().at(activeEdgeID_.get());

    length_.set(edge.getLength());
    distance_.set(edge.getDistance());
    curveness_.set(edge.getCurveness());
    minRadiusAvg_.set(edge.getMinRadiusAvg());
    minRadiusStdDeviation_.set(edge.getMinRadiusStdDeviation());
    avgRadiusAvg_.set(edge.getAvgRadiusAvg());
    avgRadiusStdDeviation_.set(edge.getAvgRadiusStdDeviation());
    maxRadiusAvg_.set(edge.getMaxRadiusAvg());
    maxRadiusStdDeviation_.set(edge.getMaxRadiusStdDeviation());
    roundnessAvg_.set(edge.getRoundnessAvg());
    roundnessStdDeviation_.set(edge.getRoundnessStdDeviation());
    avgCrossSection_.set(edge.getAvgCrossSection());
    volume_.set(edge.getVolume());
    numEdgesLeftNode_.set(edge.getNode1().getDegree());
    numEdgesRightNode_.set(edge.getNode2().getDegree());
    numSkeletonVoxels_.set(edge.getVoxels().size());

    // Write to file
    if(exportForced_) {
        exportForced_ = false;
        exportToFile(edge);
    }

    // Write to plot port
    std::unique_ptr<PlotData> plotData(new PlotData(0, 5));
    plotData->setColumnLabel(0, "Voxel Number");
    plotData->setColumnLabel(1, "MinRadius");
    plotData->setColumnLabel(2, "AvgRadius");
    plotData->setColumnLabel(3, "MaxRadius");
    plotData->setColumnLabel(4, "Roundness");

    int voxelNumber=0;
    for(const VesselSkeletonVoxel& voxel : edge.getVoxels()) {
        if(voxel.hasValidData()) { // Skip invalid voxels...
            std::vector<PlotCellValue> v(5);
            v.at(0) = PlotCellValue(static_cast<plot_t>(voxelNumber));
            v.at(1) = PlotCellValue(static_cast<plot_t>(voxel.minDistToSurface_));
            v.at(2) = PlotCellValue(static_cast<plot_t>(voxel.avgDistToSurface_));
            v.at(3) = PlotCellValue(static_cast<plot_t>(voxel.maxDistToSurface_));
            v.at(4) = PlotCellValue(static_cast<plot_t>(voxel.roundness()));
            bool inserted = plotData->insert(v);
            if (!inserted) {
                LERROR("Could not insert data into plot");
                break;
            }
        }
        ++voxelNumber;
    }
    plotOutport_.setData(plotData.release());

}
void VesselGraphStatPlotter::exportToFile(const VesselGraphEdge& edge) {
    if(exportFilePath_.get().empty()) {
        LERROR("No export file path");
        return;
    }
    try {
        CSVWriter<float, float, float> writer(exportFilePath_.get());
        writer.writeHeader("minDist", "avgDist", "maxdist");
        for(const VesselSkeletonVoxel& voxel : edge.getVoxels()) {
            writer.write(voxel.minDistToSurface_, voxel.avgDistToSurface_, voxel.maxDistToSurface_);
        }
    } catch(tgt::IOException& e) {
        LERROR("Error opening file " << exportFilePath_.get() << ": " << e.what());
        return;
    }
}
void VesselGraphStatPlotter::adaptToNewInput() {
    const VesselGraph* graph = graphInport_.getData();
    if(!graph) {
        return;
    }
    activeEdgeID_.setMaxValue(graph->getEdges().size() - 1);
}
} // namespace voreen
