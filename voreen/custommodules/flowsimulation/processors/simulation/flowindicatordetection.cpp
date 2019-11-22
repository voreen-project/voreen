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
#include "voreen/core/utils/hashing.h"
#include "../../utils/utils.h"

namespace {

/**
 * This function returns a hash for the specified vessel graph.
 * The current implementation hashes the serialized string.
 */
std::string hashVesselGraph(const voreen::VesselGraph& vesselGraph) {
    std::stringstream stream;
    voreen::JsonSerializer serializer;
    voreen::Serializer s(serializer);
    s.serialize("vesselGraph", vesselGraph);
    serializer.write(stream);
    return voreen::VoreenHash::getHash(stream.str());
}

}

namespace voreen {

const std::string FlowIndicatorDetection::loggerCat_("voreen.flowsimulation.FlowIndicatorDetection");

FlowIndicatorDetection::FlowIndicatorDetection()
    : Processor()
    , parameterInport_(Port::INPORT, "floParametrization.inport", "Flow Parametrization Input")
    , vesselGraphPort_(Port::INPORT, "vesselgraph.inport", "Vessel Graph")
    , volumePort_(Port::INPORT, "volume.inport", "Velocity Data Port (Optional)")
    , parameterOutport_(Port::OUTPORT, "flowParametrization.outport", "Flow Parametrization Output")
    , indicatorType_("flowType", "Flow Type")
    , flowProfile_("flowProfile", "Flow Profile")
    , startPhaseFunction_("startPhaseFunction", "Start Phase Function")
    , startPhaseDuration_("startPhaseDuration", "Start Phase Duration (s)", 0.0f, 0.0f, 20.0f)
    , radius_("radius", "Radius (mm)", 1.0f, 0.0f, 10.0f)
    , targetVelocity_("targetVelocity", "Target Velocity (mm/s)", 0.0f, 0.0f, 1000.0f)
    , centerlinePosition_("position", "Position", 0, 0, 0)
    , invertDirection_("invertDirection", "Invert Direction", false)
    , cloneFlowIndicator_("cloneFlowIndicator", "Clone Flow Indicator")
    , removeFlowIndicator_("removeFlowIndicator", "Remove Flow Indicator")
    , flowIndicatorTable_("flowIndicatorTable", "Flow Indicators", 2, Processor::VALID)
    , resetFlowIndicators_("resetFlowIndicators", "Reset Flow Indicators")
    , angleThreshold_("angleThreshold", "Angle Threshold", 15, 0, 90, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
    , triggertBySelection_(false)
{
    addPort(parameterInport_);
    addPort(vesselGraphPort_);
        ON_CHANGE_LAMBDA(vesselGraphPort_, [this] {
            detectFlowIndicators(false);
        });
    addPort(volumePort_);
    volumePort_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(parameterOutport_);

    addProperty(indicatorType_);
        indicatorType_.addOption("candidate", "Candidate", FlowIndicatorType::FIT_CANDIDATE);
        indicatorType_.addOption("generator", "Flow Generator", FlowIndicatorType::FIT_GENERATOR);
        indicatorType_.addOption("pressure", "Pressure Boundary", FlowIndicatorType::FIT_PRESSURE);
        indicatorType_.addOption("measure", "Flux Measure", FlowIndicatorType::FIT_MEASURE);
        indicatorType_.setGroupID("indicator");
        ON_CHANGE(indicatorType_, FlowIndicatorDetection, onIndicatorConfigChange);
        ON_CHANGE(indicatorType_, FlowIndicatorDetection, buildTable);
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
        startPhaseDuration_.setGroupID("indicator");
        ON_CHANGE(startPhaseDuration_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(radius_);
        radius_.setGroupID("indicator");
        ON_CHANGE(radius_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(targetVelocity_);
        targetVelocity_.setGroupID("indicator");
        ON_CHANGE(targetVelocity_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(centerlinePosition_);
        centerlinePosition_.setGroupID("indicator");
        ON_CHANGE(centerlinePosition_, FlowIndicatorDetection, onIndicatorSettingsChange);
    addProperty(invertDirection_);
        invertDirection_.setGroupID("indicator");
        ON_CHANGE(invertDirection_, FlowIndicatorDetection, onIndicatorSettingsChange);
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
        //flowIndicatorTable_.setColumnLabel(2, "Profile");
        //flowIndicatorTable_.setColumnLabel(3, "St. Ph. Fun.");
        //flowIndicatorTable_.setColumnLabel(4, "St. Ph. Dur.");
        //flowIndicatorTable_.setColumnLabel(5, "Radius");
        ON_CHANGE(flowIndicatorTable_, FlowIndicatorDetection, updateIndicatorUI);
    addProperty(resetFlowIndicators_);
        ON_CHANGE_LAMBDA(resetFlowIndicators_, [this] {
            detectFlowIndicators(true);
        });
    addProperty(angleThreshold_);

    // Setup property visibility.
    updateIndicatorUI();
}

void FlowIndicatorDetection::adjustPropertiesToInput() {
    if(!vesselGraphPort_.hasData())
        return;

    radius_.setMaxValue(tgt::length(vesselGraphPort_.getData()->getBounds().diagonal() / 2.0f));
}

void FlowIndicatorDetection::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("flowIndicators", flowIndicators_);
    s.serializeBinaryBlob("flowIndicatorSettings", flowIndicatorSettings_);
    s.serialize("vesselGraphHash", vesselGraphHash_);
}

void FlowIndicatorDetection::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.deserialize("flowIndicators", flowIndicators_);
    s.deserializeBinaryBlob("flowIndicatorSettings", flowIndicatorSettings_);
    s.deserialize("vesselGraphHash", vesselGraphHash_);

    buildTable();
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

    // If we have deserialized flow indicators, use them!
    if(!flowIndicators_.empty()) {
        return true;
    }

    if(!vesselGraphPort_.isReady()) {
        setNotReadyErrorMessage("No input vessel graph");
        return false;
    }

    // Reference volume is optional.

    return true;
}

void FlowIndicatorDetection::process() {

    // Here, we just set the output according to our currently set indicators.
    FlowParameterSetEnsemble* flowParameterSetEnsemble = new FlowParameterSetEnsemble(*parameterInport_.getData());

    for(const FlowIndicator& indicator : flowIndicators_) {
        // Candidates won't be added until an automatic or manual selection was performed.
        if(indicator.type_ != FIT_CANDIDATE) {
            flowParameterSetEnsemble->addFlowIndicator(indicator);
        }
    }

    parameterOutport_.setData(flowParameterSetEnsemble);
}

void FlowIndicatorDetection::updateIndicatorUI() {
    bool validSelection = flowIndicatorTable_.getNumRows() > 0 && flowIndicatorTable_.getSelectedRowIndex() >= 0;
    if(validSelection) {
        triggertBySelection_ = true;

        size_t index = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
        const FlowIndicator& indicator = flowIndicators_.at(index);

        indicatorType_.selectByValue(indicator.type_);
        flowProfile_.selectByValue(indicator.flowProfile_);
        startPhaseFunction_.selectByValue(indicator.startPhaseFunction_);
        startPhaseDuration_.set(indicator.startPhaseDuration_);
        radius_.set(indicator.radius_);
        targetVelocity_.set(indicator.targetVelocity_);

        if(vesselGraphPort_.hasData()) {
            const FlowIndicatorSettings& settings = flowIndicatorSettings_.at(index);
            const VesselGraph* vesselGraph = vesselGraphPort_.getData();
            centerlinePosition_.setMaxValue(vesselGraph->getEdge(settings.edgeId_).getVoxels().size());
            centerlinePosition_.set(settings.centerlinePosition_);
            invertDirection_.set(settings.invertDirection_);
        }

        triggertBySelection_ = false;
    }

    indicatorType_.setReadOnlyFlag(!validSelection);
    flowProfile_.setReadOnlyFlag(!validSelection);
    startPhaseFunction_.setReadOnlyFlag(!validSelection);
    startPhaseDuration_.setReadOnlyFlag(!validSelection);
    radius_.setReadOnlyFlag(!validSelection);
    targetVelocity_.setReadOnlyFlag(!validSelection);

    bool settingsEditable = validSelection && vesselGraphPort_.hasData();
    centerlinePosition_.setReadOnlyFlag(!settingsEditable);
    invertDirection_.setReadOnlyFlag(!settingsEditable);
    cloneFlowIndicator_.setReadOnlyFlag(!settingsEditable);
    removeFlowIndicator_.setReadOnlyFlag(!settingsEditable);
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
    }
}

void FlowIndicatorDetection::onIndicatorSettingsChange() {

    // Ignore calls being triggered during selection callback.
    if(triggertBySelection_)
        return;

    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {

        FlowIndicatorSettings& settings = flowIndicatorSettings_[indicatorIdx];
        settings.centerlinePosition_ = centerlinePosition_.get();
        settings.invertDirection_ = invertDirection_.get();

        // Reinitialize indicator, but keep type and id!
        FlowIndicator& indicator = flowIndicators_[indicatorIdx];
        FlowIndicatorType type = indicator.type_;
        int id = indicator.id_;
        indicator = initializeIndicator(settings);

        // Restore config.
        indicator.type_ = type;
        indicator.id_ = id;

        // Update UI.
        updateIndicatorUI();
    }
}

void FlowIndicatorDetection::detectFlowIndicators(bool forced) {

    const VesselGraph* vesselGraph = vesselGraphPort_.getData();
    if(!vesselGraph) {
        return;
    }

    // If input didn't change, we can skip the computation.
    std::string vesselGraphHash = hashVesselGraph(*vesselGraph);
    if(!forced && vesselGraphHash == vesselGraphHash_) {
        updateIndicatorUI();
        return;
    }

    // Discard old settings and recalculate data.
    vesselGraphHash_ = vesselGraphHash;
    flowIndicators_.clear();
    flowIndicatorSettings_.clear();

    for(const VesselGraphNode& node : vesselGraph->getNodes()) {
        // Look for end-nodes.
        if(node.getDegree() == 1) {

            // Generate initial settings.
            FlowIndicatorSettings settings;
            settings.nodeId_ = node.getID();
            settings.edgeId_ = node.getEdges().back().get().getID();
            settings.centerlinePosition_ = 0;
            settings.invertDirection_ = false;
            flowIndicatorSettings_.push_back(settings);

            // Initialize indicator according to those settings.
            FlowIndicator indicator = initializeIndicator(settings);
            if(indicator.type_ != FIT_INVALID) {
                flowIndicators_.push_back(indicator);
            }
        }
    }

    buildTable();
}

FlowIndicatorType FlowIndicatorDetection::estimateType(const FlowIndicator& indicator, const tgt::vec3& velocity) const {
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
        const FlowIndicatorSettings& settings = flowIndicatorSettings_.at(indicatorIdx);
        FlowIndicator clonedIndicator = initializeIndicator(settings);
        clonedIndicator.type_ = FIT_MEASURE; // Change type to MEASURE.
        flowIndicators_.push_back(clonedIndicator);
        flowIndicatorSettings_.push_back(settings);

        buildTable();
    }
}

void FlowIndicatorDetection::onRemoveFlowIndicator() {
    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {
        flowIndicators_.erase(flowIndicators_.begin() + indicatorIdx);
        flowIndicatorSettings_.erase(flowIndicatorSettings_.begin() + indicatorIdx);

        buildTable();
    }
}

FlowIndicator FlowIndicatorDetection::initializeIndicator(const FlowIndicatorSettings& settings) {

    const VesselGraph* vesselGraph = vesselGraphPort_.getData();
    const VesselGraphNode& node = vesselGraph->getNode(settings.nodeId_);
    const VesselGraphEdge& edge = vesselGraph->getEdge(settings.edgeId_);

    size_t numVoxels = edge.getVoxels().size();
    if(numVoxels == 0) {
        return FlowIndicator(); // Invalid by default.
    }

    size_t mid = std::min<size_t>(settings.centerlinePosition_, numVoxels - 1);
    size_t num = 2; // Number of reference nodes in both directions.

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

    // Default is candidate. However, we set initial values as if it was a flow generator,
    // since it might be classified as one later.
    FlowIndicator indicator;
    indicator.id_ = flowIndicators_.size() + FlowParameterSetEnsemble::getFlowIndicatorIdOffset(); // TODO: The should be determined by FlowParameterSetEnsemble.
    indicator.type_ = FlowIndicatorType::FIT_CANDIDATE;
    indicator.flowProfile_ = FlowProfile::FP_POISEUILLE;
    indicator.startPhaseFunction_ = FlowStartPhase::FSP_SINUS;
    indicator.startPhaseDuration_ = 0.2f;

    indicator.center_ = ref->pos_;
    indicator.normal_ = tgt::normalize(back->pos_ - front->pos_);
    indicator.radius_ = radius;

    if(settings.invertDirection_) {
        indicator.normal_ *= -1.0f;
    }

    // Estimate velocity and type.
    tgt::vec3 velocity = utils::sampleDisk(volumePort_.getData(), indicator.center_, indicator.normal_, indicator.radius_);
    indicator.type_ = estimateType(indicator, velocity);
    indicator.targetVelocity_ = tgt::length(velocity);

    return indicator;
}

void FlowIndicatorDetection::buildTable() {
    int selectedIndex = flowIndicatorTable_.getSelectedRowIndex();
    flowIndicatorTable_.reset();

    for(const FlowIndicator& indicator : flowIndicators_) {
        std::vector<std::string> row(flowIndicatorTable_.getNumColumns());
        row[0] = std::to_string(indicator.id_);
        row[1] = indicator.type_ == FlowIndicatorType::FIT_GENERATOR ? "Generator" :
                 (indicator.type_ == FlowIndicatorType::FIT_PRESSURE ? "Pressure" :
                  (indicator.type_ == FlowIndicatorType::FIT_MEASURE ? "Measure" : "Candidate"));
        /*
        row[2] = indicator.flowProfile_ == FlowProfile ::FP_POISEUILLE ? "Poiseuille" :
                 (indicator.flowProfile_ == FlowProfile::FP_POWERLAW ? "Powerlaw" :
                  (indicator.flowProfile_ == FlowProfile::FP_CONSTANT ? "Constant" : "None"));
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
