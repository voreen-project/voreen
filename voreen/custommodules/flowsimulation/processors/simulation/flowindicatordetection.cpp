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
#include "../../utils/utils.h"

namespace voreen {

const std::string FlowIndicatorDetection::loggerCat_("voreen.flowsimulation.FlowIndicatorDetection");

FlowIndicatorDetection::FlowIndicatorDetection()
    : Processor()
    , parameterInport_(Port::INPORT, "floParametrization.inport", "Flow Parametrization Input")
    , vesselGraphPort_(Port::INPORT, "vesselgraph.inport", "Vessel Graph")
    , volumePort_(Port::INPORT, "volume.inport", "Velocity Data Port (Optional)")
    , parameterOutport(Port::OUTPORT, "flowParametrization.outport", "Flow Parametrization Output")
    , indicatorType_("flowType", "Flow Type")
    , flowProfile_("flowProfile", "Flow Profile")
    , startPhaseFunction_("startPhaseFunction", "Start Phase Function")
    , startPhaseDuration_("startPhaseDuration", "Start Phase Duration (s)", 0.0f, 0.0f, 20.0f)
    , radius_("radius", "Radius", 1.0f, 0.0f, 10.0f)
    , targetVelocity_("targetVelocity", "Target Velocity", 0.0f, 0.0f, 1000.0f)
    , voxelId_("voxelId", "Voxel ID", 0, 0, 0)
    , cloneFlowIndicator_("cloneFlowIndicator", "Clone Flow Indicator")
    , removeFlowIndicator_("removeFlowIndicator", "Remove Flow Indicator")
    , flowIndicatorTable_("flowIndicators", "Flow Indicators", 3)
    , angleThreshold_("angleThreshold", "Angle Threshold", 15, 0, 90, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
    , triggertBySelection_(false)
{
    addPort(parameterInport_);
    addPort(vesselGraphPort_);
    ON_CHANGE(vesselGraphPort_, FlowIndicatorDetection, onInputChange);
    addPort(volumePort_);
    volumePort_.addCondition(new PortConditionVolumeChannelCount(3));
    ON_CHANGE(volumePort_, FlowIndicatorDetection, onInputChange);
    addPort(parameterOutport);

    addProperty(indicatorType_);
        indicatorType_.addOption("candidate", "Candidate", FlowIndicatorType::FIT_CANDIDATE);
        indicatorType_.addOption("generator", "Flow Generator", FlowIndicatorType::FIT_GENERATOR);
        indicatorType_.addOption("pressure", "Pressure Boundary", FlowIndicatorType::FIT_PRESSURE);
        indicatorType_.addOption("measure", "Flux Measure", FlowIndicatorType::FIT_MEASURE);
        indicatorType_.setGroupID("indicator");
        ON_CHANGE(indicatorType_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(flowProfile_);
        flowProfile_.addOption("none", "None", FlowProfile::FP_NONE);
        flowProfile_.addOption("poiseuille", "Poiseuille", FlowProfile::FP_POISEUILLE);
        flowProfile_.addOption("powerlaw", "Powerlaw", FlowProfile::FP_POWERLAW);
        //flowProfile_.addOption("constant", "constant", FlowProfile::FP_CONSTANT); // Not really useful.
        flowProfile_.setGroupID("indicator");
        ON_CHANGE(flowProfile_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(startPhaseFunction_);
        startPhaseFunction_.addOption("none", "None", FlowStartPhase::FSP_NONE); // get's selected automatically
        startPhaseFunction_.addOption("constant", "Constant", FlowStartPhase ::FSP_CONSTANT);
        startPhaseFunction_.addOption("sinus", "Sinus", FlowStartPhase::FSP_SINUS);
        startPhaseFunction_.setGroupID("indicator");
        ON_CHANGE(startPhaseFunction_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(startPhaseDuration_);
        startPhaseDuration_.setTracking(false);
        startPhaseDuration_.setGroupID("indicator");
        ON_CHANGE(startPhaseDuration_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(radius_);
        radius_.setTracking(false);
        radius_.setGroupID("indicator");
        ON_CHANGE(radius_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(targetVelocity_);
        targetVelocity_.setTracking(false);
        targetVelocity_.setGroupID("indicator");
        ON_CHANGE(targetVelocity_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(voxelId_);
        voxelId_.setTracking(false);
        voxelId_.setGroupID("indicator");
        ON_CHANGE(voxelId_, FlowIndicatorDetection, onIndicatorPositionChange);
    addProperty(cloneFlowIndicator_);
        cloneFlowIndicator_.setGroupID("indicator");
        ON_CHANGE(cloneFlowIndicator_, FlowIndicatorDetection, onCloneFlowIndicator);
    addProperty(removeFlowIndicator_);
        removeFlowIndicator_.setGroupID("indicator");
        ON_CHANGE(removeFlowIndicator_, FlowIndicatorDetection, onRemoveFlowIndicator);
    setPropertyGroupGuiName("indicator", "Indicator");

    addProperty(flowIndicatorTable_);
    flowIndicatorTable_.setColumnLabel(0, "ID");
    flowIndicatorTable_.setColumnLabel(1, "Type");
    flowIndicatorTable_.setColumnLabel(2, "Profile");
    //flowIndicatorTable_.setColumnLabel(3, "St. Ph. Fun.");
    //flowIndicatorTable_.setColumnLabel(4, "St. Ph. Dur.");
    //flowIndicatorTable_.setColumnLabel(5, "Radius");
    ON_CHANGE(flowIndicatorTable_, FlowIndicatorDetection, onIndicatorSelectionChange);

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

    if(!parameterInport_.isReady()) {
        setNotReadyErrorMessage("No input parametrization");
        return false;
    }

    if(!vesselGraphPort_.isReady()) {
        setNotReadyErrorMessage("No input vessel graph");
        return false;
    }

    // Reference volume is optional.

    return true;
}

void FlowIndicatorDetection::process() {

    FlowParameterSetEnsemble* flowParameterSetEnsemble = new FlowParameterSetEnsemble(*parameterInport_.getData());

    for(const FlowIndicator& indicator : flowIndicators_) {
        // Candidates wont be added until an automatic or manual selection was performed.
        if(indicator.type_ != FIT_CANDIDATE) {
            flowParameterSetEnsemble->addFlowIndicator(indicator);
        }
    }

    parameterOutport.setData(flowParameterSetEnsemble);
}

void FlowIndicatorDetection::onIndicatorSelectionChange() {
    bool validSelection = flowIndicatorTable_.getNumRows() > 0 && flowIndicatorTable_.getSelectedRowIndex() >= 0;
    if(validSelection) {
        triggertBySelection_ = true;
        size_t index = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
        indicatorType_.selectByValue(flowIndicators_.at(index).type_);
        flowProfile_.selectByValue(flowIndicators_.at(index).flowProfile_);
        startPhaseFunction_.selectByValue(flowIndicators_.at(index).startPhaseFunction_);
        startPhaseDuration_.set(flowIndicators_.at(index).startPhaseDuration_);
        radius_.set(flowIndicators_.at(index).radius_);
        targetVelocity_.set(flowIndicators_.at(index).targetVelocity_);

        const IndicatorSettings& settings = indicatorSettings_.at(index);
        const VesselGraph* vesselGraph = vesselGraphPort_.getData();
        voxelId_.setMaxValue(vesselGraph->getEdge(settings.edgeId).getVoxels().size());
        voxelId_.set(settings.voxelId);
        triggertBySelection_ = false;
    }

    indicatorType_.setReadOnlyFlag(!validSelection);
    flowProfile_.setReadOnlyFlag(!validSelection);
    startPhaseFunction_.setReadOnlyFlag(!validSelection);
    startPhaseDuration_.setReadOnlyFlag(!validSelection);
    radius_.setReadOnlyFlag(!validSelection);
    targetVelocity_.setReadOnlyFlag(!validSelection);
    voxelId_.setReadOnlyFlag(!validSelection);
    cloneFlowIndicator_.setReadOnlyFlag(!validSelection);
    removeFlowIndicator_.setReadOnlyFlag(!validSelection);

    //setPropertyGroupVisible("indicator", flowIndicatorTable_.getSelectedRowIndex() >= 0);
}

void FlowIndicatorDetection::onIndicatorConfigChange() {

    // Ignore calls being triggered during selection callback.
    if(triggertBySelection_)
        return;

    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {

        FlowIndicator& indicator = flowIndicators_[indicatorIdx];
        indicator.type_ = indicatorType_.getValue();
        indicator.flowProfile_ = flowProfile_.getValue();
        indicator.startPhaseFunction_ = startPhaseFunction_.getValue();
        indicator.startPhaseDuration_ = startPhaseDuration_.get();
        indicator.radius_ = radius_.get();
        indicator.targetVelocity_ = targetVelocity_.get();

        buildTable();
    }
}

void FlowIndicatorDetection::onIndicatorPositionChange() {

    // Ignore calls being triggered during selection callback.
    if(triggertBySelection_)
        return;

    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {

        FlowIndicator& indicator = flowIndicators_[indicatorIdx];
        IndicatorSettings& settings = indicatorSettings_[indicatorIdx];
        settings.voxelId = voxelId_.get();
        updateIndicator(indicator, settings);

        buildTable();
    }
}

void FlowIndicatorDetection::onInputChange() {

    // TODO: Dont ALWAYS discard deserialized values.
    flowIndicators_.clear();
    indicatorSettings_.clear();

    const FlowParameterSetEnsemble* flowParameterSetEnsemble = parameterInport_.getData();
    const VesselGraph* vesselGraph = vesselGraphPort_.getData();
    if(!flowParameterSetEnsemble || !vesselGraph) {
        parameterOutport.clear();
        buildTable();
        return;
    }

    for(const VesselGraphNode& node : vesselGraph->getNodes()) {
        // Look for end-nodes.
        if(node.getDegree() == 1) {

            IndicatorSettings settings;
            settings.nodeId = node.getID();
            settings.edgeId = node.getEdges().back().get().getID();
            settings.voxelId = 0;

            // Default is candidate. However, we set initial values as if it was a flow generator,
            // since it might be classified as one later.
            FlowIndicator indicator;
            indicator.flowProfile_ = FlowProfile::FP_POISEUILLE;
            indicator.startPhaseFunction_ = FlowStartPhase::FSP_SINUS;
            indicator.startPhaseDuration_ = flowParameterSetEnsemble->getSimulationTime() * 0.25f;
            updateIndicator(indicator, settings);

            flowIndicators_.push_back(indicator);
            indicatorSettings_.push_back(settings);
        }
    }

    buildTable();
}

FlowIndicatorType FlowIndicatorDetection::estimateType(const FlowIndicator& indicator, const tgt::vec3& velocity) const {

    // Do not overwrite already set indicator types!
    if(indicator.type_ != FlowIndicatorType::FIT_CANDIDATE) {
        return indicator.type_;
    }

    if (velocity != tgt::vec3::zero) {
        tgt::vec3 v = tgt::normalize(velocity);
        float threshold = tgt::deg2rad(static_cast<float>(angleThreshold_.get()));
        float angle = std::acos(tgt::dot(v, indicator.normal_));
        if (angle < threshold) {
            return FlowIndicatorType::FIT_GENERATOR;
        }
        else if (tgt::PIf - angle < threshold) {
            return FlowIndicatorType::FIT_PRESSURE;
        }
    }

    return FlowIndicatorType::FIT_CANDIDATE;
}

void FlowIndicatorDetection::onCloneFlowIndicator() {
    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {
        FlowIndicator clonedIndicator = flowIndicators_[indicatorIdx];
        clonedIndicator.type_ = FIT_CANDIDATE; // This now becomes a candidate.
        flowIndicators_.push_back(clonedIndicator);
        indicatorSettings_.push_back(indicatorSettings_[indicatorIdx]);

        buildTable();
    }
}

void FlowIndicatorDetection::onRemoveFlowIndicator() {
    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {
        flowIndicators_.erase(flowIndicators_.begin() + indicatorIdx);
        indicatorSettings_.erase(indicatorSettings_.begin() + indicatorIdx);

        buildTable();
    }
}

void FlowIndicatorDetection::updateIndicator(FlowIndicator& indicator, IndicatorSettings& settings) {
    const VesselGraph* vesselGraph = vesselGraphPort_.getData();
    const VesselGraphNode& node = vesselGraph->getNode(settings.nodeId);
    const VesselGraphEdge& edge = vesselGraph->getEdge(settings.edgeId);

    size_t numVoxels = edge.getVoxels().size();
    if (numVoxels == 0) {
        //tgtAssert(false, "No voxels assigned to edge");
        return;
    }

    size_t mid = std::min<size_t>(voxelId_.get(), numVoxels - 1);
    size_t num = 2; // Number of reference nodes.

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

    indicator.center_ = ref->pos_;
    indicator.normal_ = tgt::normalize(back->pos_ - front->pos_);
    indicator.radius_ = radius;

    tgt::vec3 velocity = utils::sampleDisk(volumePort_.getData(), indicator.center_, indicator.normal_, indicator.radius_);
    indicator.type_ = estimateType(indicator, velocity);
    indicator.targetVelocity_ = tgt::length(velocity);
}

void FlowIndicatorDetection::buildTable() {
    int selectedIndex = flowIndicatorTable_.getSelectedRowIndex();
    flowIndicatorTable_.reset();

    for(const FlowIndicator& indicator : flowIndicators_) {
        std::vector<std::string> row(flowIndicatorTable_.getNumColumns());
        row[0] = std::to_string(flowIndicatorTable_.getNumRows() + 1); // Very simple ID.
        row[1] = indicator.type_ == FlowIndicatorType::FIT_GENERATOR ? "Generator" :
                 (indicator.type_ == FlowIndicatorType::FIT_PRESSURE ? "Pressure" :
                  (indicator.type_ == FlowIndicatorType::FIT_MEASURE ? "Measure" : "Candidate"));
        row[2] = indicator.flowProfile_ == FlowProfile ::FP_POISEUILLE ? "Poiseuille" :
                 (indicator.flowProfile_ == FlowProfile::FP_POWERLAW ? "Powerlaw" :
                  (indicator.flowProfile_ == FlowProfile::FP_CONSTANT ? "Constant" : "None"));
        /*
        row[3] = indicator.startPhaseFunction_ == FlowStartPhase::FSP_CONSTANT ? "Constant" :
                 (indicator.startPhaseFunction_ == FlowStartPhase::FSP_SINUS ? "Sinus" : "None");
        row[4] = std::to_string(indicator.startPhaseDuration_);
        row[5] = std::to_string(indicator.radius_);
        row[6] = "(" + std::to_string(indicator.center_.x) + ", " + std::to_string(indicator.center_.y) + ", " +
                 std::to_string(indicator.center_.z) + ")";
        row[7] = "(" + std::to_string(indicator.normal_.x) + ", " + std::to_string(indicator.normal_.y) + ", " +
                 std::to_string(indicator.normal_.z) + ")";
        */
        flowIndicatorTable_.addRow(row);
    }

    if(selectedIndex < static_cast<int>(flowIndicatorTable_.getNumRows())) {
        flowIndicatorTable_.setSelectedRowIndex(selectedIndex);
    }
}

}   // namespace
