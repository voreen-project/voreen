/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "flowindicatordetection.h"

namespace voreen {

const std::string FlowIndicatorDetection::loggerCat_("voreen.flowreen.FlowIndicatorDetection");

FlowIndicatorDetection::FlowIndicatorDetection()
    : Processor()
    , vesselGraphPort_(Port::INPORT, "vesselgraph.inport", "Vessel Graph (Optional)")
    , flowParametrizationPort_(Port::OUTPORT, "flowParametrization.outport", "Flow Parametrization")
    , ensembleName_("ensembleName", "Ensemble Name", "test_ensemble")
    , simulationTime_("simulationTime", "Simulation Time (s)", 2.0f, 0.1f, 10.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution (ms)", 3.1f, 1.0f, 30.0f)
    , spatialResolution_("spatialResolution", "Spatial Resolution", 128, 32, 1024)
    , flowDirection_("flowDirection", "Flow Direction")
    , radius_("radius", "Radius", 1.0f, 0.0f, 10.0f)
    , flowIndicatorTable_("flowIndicators", "Flow Indicators", 4)
{
    addPort(vesselGraphPort_);
    ON_CHANGE(vesselGraphPort_, FlowIndicatorDetection, onVesselGraphChange);
    addPort(flowParametrizationPort_);

    addProperty(ensembleName_);
        ensembleName_.setGroupID("ensemble");
    addProperty(simulationTime_);
        simulationTime_.setGroupID("ensemble");
    addProperty(temporalResolution_);
        temporalResolution_.setGroupID("ensemble");
    addProperty(spatialResolution_);
        spatialResolution_.setGroupID("ensemble");
    setPropertyGroupGuiName("ensemble", "Ensemble");

    addProperty(flowDirection_);
        flowDirection_.addOption("none", "NONE", FlowDirection::NONE);
        flowDirection_.addOption("in", "IN", FlowDirection::IN);
        flowDirection_.addOption("out", "OUT", FlowDirection::OUT);
        flowDirection_.setGroupID("indicator");
        ON_CHANGE(flowDirection_, FlowIndicatorDetection, onConfigChange);
    //addProperty(radius_);
        radius_.setGroupID("indicator");
        ON_CHANGE(radius_, FlowIndicatorDetection, onConfigChange);
    setPropertyGroupGuiName("indicator", "Indicator");

    addProperty(flowIndicatorTable_);
    flowIndicatorTable_.setColumnLabel(0, "Dir.");
    flowIndicatorTable_.setColumnLabel(1, "Center");
    flowIndicatorTable_.setColumnLabel(2, "Normal");
    flowIndicatorTable_.setColumnLabel(3, "Radius");
    ON_CHANGE(flowIndicatorTable_, FlowIndicatorDetection, onSelectionChange);
}

void FlowIndicatorDetection::adjustPropertiesToInput() {

    if(!vesselGraphPort_.hasData())
        return;

    radius_.setMaxValue(tgt::length(vesselGraphPort_.getData()->getBounds().diagonal() / 2.0f));
}

void FlowIndicatorDetection::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("flowIndicators", flowIndicators_);
}

void FlowIndicatorDetection::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.optionalDeserialize("flowIndicators", flowIndicators_, flowIndicators_);
}

void FlowIndicatorDetection::process() {

    FlowParametrizationList* flowParametrizationList = new FlowParametrizationList(ensembleName_.get());
    flowParametrizationList->setSimulationTime(simulationTime_.get());
    flowParametrizationList->setTemporalResolution(temporalResolution_.get());
    flowParametrizationList->setSpatialResolution(spatialResolution_.get());

    for(const FlowIndicator& indicator : flowIndicators_) {
        // NONE means invalid or not being selected for output.
        if(indicator.direction_ != NONE) {
            flowParametrizationList->addFlowIndicator(indicator);
        }
    }

    flowParametrizationPort_.setData(flowParametrizationList);
}

void FlowIndicatorDetection::onSelectionChange() {
    flowDirection_.setReadOnlyFlag(flowIndicatorTable_.getSelectedRowIndex() < 0);
    //setPropertyGroupVisible("indicator", flowIndicatorTable_.getSelectedRowIndex() >= 0);
}

void FlowIndicatorDetection::onConfigChange() {
    if(flowIndicatorTable_.getNumRows() > 0 &&
       flowIndicatorTable_.getSelectedRowIndex() >= 0 &&
       flowIndicatorTable_.getSelectedRowIndex() < static_cast<int>(flowIndicators_.size()) ) {

        FlowIndicator& indicator = flowIndicators_[flowIndicatorTable_.getSelectedRowIndex()];
        indicator.direction_ = flowDirection_.getValue();
        //indicator.radius_ = radius_.get(); // Estimate is quite accurate.

        buildTable();
    }
}

void FlowIndicatorDetection::onVesselGraphChange() {

    flowIndicators_.clear();

    const VesselGraph* vesselGraph = vesselGraphPort_.getData();
    if(!vesselGraph) {
        flowParametrizationPort_.clear();
        buildTable();
        return;
    }

    for(const VesselGraphNode& node : vesselGraph->getNodes()) {
        // Look for end-nodes.
        if(node.getDegree() == 1) {

            const VesselGraphEdge& edge = node.getEdges().back().get();

            const VesselSkeletonVoxel& end = edge.getVoxels().back();
            size_t n = std::min<size_t>(5, edge.getVoxels().size());
            const VesselSkeletonVoxel& ref = edge.getVoxels().at(edge.getVoxels().size() - n);

            FlowIndicator indicator;
            indicator.center_ = end.pos_;
            indicator.normal_ = tgt::normalize(end.pos_ - ref.pos_);
            indicator.radius_ = end.avgDistToSurface_;
            indicator.direction_ = FlowDirection::NONE;
            flowIndicators_.push_back(indicator);
        }
    }

    buildTable();
}

void FlowIndicatorDetection::buildTable() {
    flowIndicatorTable_.reset();

    for(const FlowIndicator& indicator : flowIndicators_) {
        std::vector<std::string> row(4);
        row[0] = indicator.direction_ == FlowDirection::IN ? "IN" : (indicator.direction_ == FlowDirection::OUT
                                                                     ? "OUT" : "NONE");
        row[1] = "(" + std::to_string(indicator.center_.x) + ", " + std::to_string(indicator.center_.y) + ", " +
                 std::to_string(indicator.center_.z) + ")";
        row[2] = "(" + std::to_string(indicator.normal_.x) + ", " + std::to_string(indicator.normal_.y) + ", " +
                 std::to_string(indicator.normal_.z) + ")";
        row[3] = std::to_string(indicator.radius_);
        flowIndicatorTable_.addRow(row);
    }
}

}   // namespace
