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

#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

const std::string FlowIndicatorDetection::loggerCat_("voreen.flowsimulation.FlowIndicatorDetection");

FlowIndicatorDetection::FlowIndicatorDetection()
    : Processor()
    , vesselGraphPort_(Port::INPORT, "vesselgraph.inport", "Vessel Graph (Optional)")
    , volumePort_(Port::INPORT, "volume.inport", "Velocity Data Port (Optional)")
    , flowParametrizationPort_(Port::OUTPORT, "flowParametrization.outport", "Flow Parametrization")
    , ensembleName_("ensembleName", "Ensemble Name", "test_ensemble")
    , simulationTime_("simulationTime", "Simulation Time (s)", 2.0f, 0.1f, 20.0f)
    , numTimeSteps_("numTimeSteps", "Num. Output Time Steps", 50, 1, 1000)
    , outputResolution_("outputResolution", "Max. Output Resolution", 128, 32, 1024)
    , flowFeatures_("flowFeatures", "Flow Features")
    , flowDirection_("flowDirection", "Flow Direction")
    , startPhaseFunction_("startPhaseFunction", "Start Phase Function")
    , startPhaseDuration_("startPhaseDuration", "Start Phase Duration (s)", 0.0f, 0.0f, 20.0f)
    , radius_("radius", "Radius", 1.0f, 0.0f, 10.0f)
    , flowIndicatorTable_("flowIndicators", "Flow Indicators", 4)
    , firstRefNode_("firstRefNode", "First Ref. Nodes", 0, 0, 20)
    , numRefNodes_("numRefNodes", "Num. Ref. Nodes", 3, 1, 10)
    , angleThreshold_("angleThreshold", "Angle Threshold", 15, 0, 90)
    , triggertBySelection_(false)
{
    addPort(vesselGraphPort_);
    ON_CHANGE(vesselGraphPort_, FlowIndicatorDetection, onInputChange);
    addPort(volumePort_);
    volumePort_.addCondition(new PortConditionVolumeType3xFloat());
    ON_CHANGE(volumePort_, FlowIndicatorDetection, onInputChange);
    addPort(flowParametrizationPort_);

    addProperty(ensembleName_);
        ensembleName_.setGroupID("ensemble");
    addProperty(simulationTime_);
        simulationTime_.setGroupID("ensemble");
        ON_CHANGE_LAMBDA(simulationTime_, [this] {
            startPhaseDuration_.setMaxValue(simulationTime_.get());
        });
    addProperty(numTimeSteps_);
        numTimeSteps_.setGroupID("ensemble");
    addProperty(outputResolution_);
        outputResolution_.setGroupID("ensemble");

    addProperty(flowFeatures_);
    addFeature("Velocity", FT_VELOCITY);
    addFeature("Magnitude", FT_MAGNITUDE);
    addFeature("Pressure", FT_PRESSURE);
    addFeature("Wall Shear Stress", FT_WALLSHEARSTRESS);
    flowFeatures_.addInstance("Velocity"); // Default selection.
    flowFeatures_.setGroupID("ensemble");
    setPropertyGroupGuiName("ensemble", "Ensemble");

    addProperty(flowDirection_);
        flowDirection_.addOption("none", "NONE", FlowDirection::FD_NONE);
        flowDirection_.addOption("in", "IN", FlowDirection::FD_IN);
        flowDirection_.addOption("out", "OUT", FlowDirection::FD_OUT);
        flowDirection_.setGroupID("indicator");
        ON_CHANGE(flowDirection_, FlowIndicatorDetection, onConfigChange);
    addProperty(startPhaseFunction_);
        startPhaseFunction_.addOption("none", "NONE", FlowFunction::FF_NONE); // get's selected automatically
        startPhaseFunction_.addOption("constant", "CONSTANT", FlowFunction ::FF_CONSTANT);
        startPhaseFunction_.addOption("sinus", "SINUS", FlowFunction::FF_SINUS);
        startPhaseFunction_.setGroupID("indicator");
        ON_CHANGE(startPhaseFunction_, FlowIndicatorDetection, onConfigChange);
    addProperty(startPhaseDuration_);
        startPhaseDuration_.setTracking(false);
        startPhaseDuration_.setGroupID("indicator");
        ON_CHANGE(startPhaseDuration_, FlowIndicatorDetection, onConfigChange);
    addProperty(radius_);
        radius_.setTracking(false);
        radius_.setGroupID("indicator");
        ON_CHANGE(radius_, FlowIndicatorDetection, onConfigChange);
    setPropertyGroupGuiName("indicator", "Indicator");

    addProperty(flowIndicatorTable_);
    flowIndicatorTable_.setColumnLabel(0, "Dir.");
    flowIndicatorTable_.setColumnLabel(1, "St. Ph. Fun.");
    flowIndicatorTable_.setColumnLabel(2, "St. Ph. Dur.");
    flowIndicatorTable_.setColumnLabel(3, "Radius");
    //flowIndicatorTable_.setColumnLabel(4, "Center");
    //flowIndicatorTable_.setColumnLabel(5, "Normal");
    ON_CHANGE(flowIndicatorTable_, FlowIndicatorDetection, onSelectionChange);

    addProperty(firstRefNode_);
    ON_CHANGE(firstRefNode_, FlowIndicatorDetection, onInputChange);
    addProperty(numRefNodes_);
    ON_CHANGE(numRefNodes_, FlowIndicatorDetection, onInputChange);
    addProperty(angleThreshold_);
    ON_CHANGE(angleThreshold_, FlowIndicatorDetection, onInputChange);
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

bool FlowIndicatorDetection::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized");
        return false;
    }

    // Both inports are optional

    return true;
}

void FlowIndicatorDetection::process() {

    FlowParametrizationList* flowParametrizationList = new FlowParametrizationList(ensembleName_.get());
    flowParametrizationList->setSimulationTime(simulationTime_.get());
    flowParametrizationList->setNumTimeSteps(numTimeSteps_.get());
    flowParametrizationList->setOutputResolution(outputResolution_.get());

    int flowFeatures = FT_NONE;
    for(const InteractiveListProperty::Instance& instance : flowFeatures_.getInstances()) {
        flowFeatures |=  flowFeatureIds_[instance.itemId_];
    }
    flowParametrizationList->setFlowFeatures(flowFeatures);

    for(const FlowIndicator& indicator : flowIndicators_) {
        // NONE means invalid or not being selected for output.
        if(indicator.direction_ != FD_NONE) {
            flowParametrizationList->addFlowIndicator(indicator);
        }
    }

    flowParametrizationPort_.setData(flowParametrizationList);
}

void FlowIndicatorDetection::addFeature(const std::string &name, int id) {
    flowFeatures_.addItem(name);
    flowFeatureIds_.push_back(id);
}

void FlowIndicatorDetection::onSelectionChange() {
    bool validSelection = flowIndicatorTable_.getNumRows() > 0 && flowIndicatorTable_.getSelectedRowIndex() >= 0;
    if(validSelection) {
        triggertBySelection_ = true;
        size_t index = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
        flowDirection_.selectByValue(flowIndicators_.at(index).direction_);
        startPhaseFunction_.selectByValue(flowIndicators_.at(index).startPhaseFunction_);
        startPhaseDuration_.set(flowIndicators_.at(index).startPhaseDuration_);
        radius_.set(flowIndicators_.at(index).radius_);
        triggertBySelection_ = false;
    }

    flowDirection_.setReadOnlyFlag(!validSelection);
    startPhaseFunction_.setReadOnlyFlag(!validSelection);
    startPhaseDuration_.setReadOnlyFlag(!validSelection);
    radius_.setReadOnlyFlag(!validSelection);

    //setPropertyGroupVisible("indicator", flowIndicatorTable_.getSelectedRowIndex() >= 0);
}

void FlowIndicatorDetection::onConfigChange() {

    // Ignore calls being triggert during selection callback.
    if(triggertBySelection_)
        return;

    if(flowIndicatorTable_.getNumRows() > 0 &&
       flowIndicatorTable_.getSelectedRowIndex() >= 0 &&
       flowIndicatorTable_.getSelectedRowIndex() < static_cast<int>(flowIndicators_.size()) ) {

        FlowIndicator& indicator = flowIndicators_[flowIndicatorTable_.getSelectedRowIndex()];
        indicator.direction_ = flowDirection_.getValue();
        indicator.startPhaseFunction_ = startPhaseFunction_.getValue();
        indicator.startPhaseDuration_ = startPhaseDuration_.get();
        indicator.radius_ = radius_.get(); // Estimate is quite accurate.

        buildTable();
    }
}

void FlowIndicatorDetection::onInputChange() {

    flowIndicators_.clear();

    const VesselGraph* vesselGraph = vesselGraphPort_.getData();
    if(!vesselGraph) {
        flowParametrizationPort_.clear();
        buildTable();
        return;
    }

    const VolumeBase* volume = volumePort_.getData();

    for(const VesselGraphNode& node : vesselGraph->getNodes()) {
        // Look for end-nodes.
        if(node.getDegree() == 1) {

            const VesselGraphEdge& edge = node.getEdges().back().get();
            size_t numVoxels = edge.getVoxels().size();
            if (numVoxels == 0) {
                //tgtAssert(false, "No voxels assigned to edge");
                continue;
            }

            size_t mid = std::min<size_t>(firstRefNode_.get(), numVoxels-1);
            size_t num = static_cast<size_t>(numRefNodes_.get());

            size_t frontIdx = mid > num ? (mid - num) : 0;
            size_t backIdx  = std::min(mid + num, numVoxels - 1);

            std::function<size_t (size_t)> index;
            if(edge.getNode1().getID() == node.getID()) {
                index = [](size_t i) { return i; };
            }
            else {
                index = [numVoxels](size_t i) { return numVoxels - 1 - i;  };
            }

            const VesselSkeletonVoxel* ref   = &edge.getVoxels().at(index(mid));
            const VesselSkeletonVoxel* front = &edge.getVoxels().at(index(frontIdx));
            const VesselSkeletonVoxel* back  = &edge.getVoxels().at(index(backIdx));

            // Calculate average radius.
            float radius = 0.0f;
            for (size_t i = frontIdx; i <= backIdx; i++) {
                radius += edge.getVoxels().at(index(i)).avgDistToSurface_;
            }
            radius /= (backIdx - frontIdx + 1);

            FlowIndicator indicator;
            indicator.center_ = ref->pos_;
            indicator.normal_ = tgt::normalize(back->pos_ - front->pos_);
            indicator.radius_ = radius;
            indicator.direction_ = FlowDirection::FD_NONE;
            // Define default values here:
            indicator.startPhaseFunction_ = FlowFunction::FF_SINUS;
            indicator.startPhaseDuration_ = simulationTime_.get() * 0.25f;

            // Estimate flow direction based on underlying velocities.
            if (volume) {

                tgt::vec3 velocity = tgt::vec3::zero;
                if (auto velocities = dynamic_cast<const VolumeRAM_3xFloat*>(volume->getRepresentation<VolumeRAM>())) {
                    tgt::svec3 voxel = volume->getWorldToVoxelMatrix() * indicator.center_;
                    velocity = velocities->voxel(voxel);
                }

                if (velocity != tgt::vec3::zero) {
                    velocity = tgt::normalize(velocity);
                    float threshold = tgt::deg2rad(static_cast<float>(angleThreshold_.get()));
                    float angle = std::acos(tgt::dot(velocity, indicator.normal_) /
                        (tgt::length(velocity) * tgt::length(indicator.normal_)));
                    if (angle < threshold) {
                        indicator.direction_ = FlowDirection::FD_IN;
                    }
                    else if (tgt::PIf - angle < threshold) {
                        indicator.direction_ = FlowDirection::FD_OUT;
                    }
                }
            }

            flowIndicators_.push_back(indicator);
        }
    }

    buildTable();
}

void FlowIndicatorDetection::buildTable() {
    int selectedIndex = flowIndicatorTable_.getSelectedRowIndex();
    flowIndicatorTable_.reset();

    for(const FlowIndicator& indicator : flowIndicators_) {
        std::vector<std::string> row(flowIndicatorTable_.getNumColumns());
        row[0] = indicator.direction_ == FlowDirection::FD_IN ? "IN" :
                (indicator.direction_ == FlowDirection::FD_OUT ? "OUT" : "NONE");
        row[1] = indicator.startPhaseFunction_ == FlowFunction::FF_CONSTANT ? "CONSTANT" :
                (indicator.startPhaseFunction_ == FlowFunction::FF_SINUS ? "SINUS" : "NONE");
        row[2] = std::to_string(indicator.startPhaseDuration_);
        row[3] = std::to_string(indicator.radius_);
        /*
        row[4] = "(" + std::to_string(indicator.center_.x) + ", " + std::to_string(indicator.center_.y) + ", " +
                 std::to_string(indicator.center_.z) + ")";
        row[5] = "(" + std::to_string(indicator.normal_.x) + ", " + std::to_string(indicator.normal_.y) + ", " +
                 std::to_string(indicator.normal_.z) + ")";
        */
        flowIndicatorTable_.addRow(row);
    }

    if(selectedIndex < static_cast<int>(flowIndicatorTable_.getNumRows())) {
        flowIndicatorTable_.setSelectedRowIndex(selectedIndex);
    }
}

}   // namespace
