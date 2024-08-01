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

#include "flowindicatordetection.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/utils/hashing.h"
#include "../../utils/serializationhelper.h"
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

FlowIndicatorDetection::FlowIndicatorSettings::FlowIndicatorSettings()
    : FlowIndicatorSettings(VGNodeID::INVALID, VGEdgeID::INVALID)
{
}

FlowIndicatorDetection::FlowIndicatorSettings::FlowIndicatorSettings(VGNodeID nodeId, VGEdgeID edgeId)
    : nodeId_(nodeId)
    , edgeId_(edgeId)
    , centerlinePosition_(0)
    , relativeRadiusCorrection_(1.0f)
    , invertDirection_(false)
    , forceAxisAlignment_(false)
    , velocityCurveType_("sinus")
    , velocityCurveDuration_(0.25f)
    , targetVelocity_(0.0f)
    , velocityCurveFile_()
    , velocityCurvePeriodic_(false)
    , velocityCurveScale_(1.0f)
{
}

void FlowIndicatorDetection::FlowIndicatorSettings::serialize(Serializer& s) const {
    s.serialize("nodeId", nodeId_.raw());
    s.serialize("edgeId", edgeId_.raw());
    s.serialize("centerlinePosition", centerlinePosition_);
    s.serialize("relativeRadiusCorrection", relativeRadiusCorrection_);
    s.serialize("invertDirection", invertDirection_);
    s.serialize("forceAxisAlignment", forceAxisAlignment_);
    s.serialize("velocityCurveType", velocityCurveType_);
    s.serialize("velocityCurveDuration", velocityCurveDuration_);
    s.serialize("targetVelocity", targetVelocity_);
    s.serialize("velocityCurveFile", velocityCurveFile_);
    s.serialize("velocityCurvePeriodic", velocityCurvePeriodic_);
}

void FlowIndicatorDetection::FlowIndicatorSettings::deserialize(Deserializer& s) {
    uint32_t nodeId = -1;
    s.deserialize("nodeId", nodeId);
    nodeId_ = VGNodeID(nodeId);
    uint32_t edgeId = -1;
    s.deserialize("edgeId", edgeId);
    edgeId_ = VGEdgeID(edgeId);
    s.deserialize("centerlinePosition", centerlinePosition_);
    s.optionalDeserialize("relativeRadiusCorrection", relativeRadiusCorrection_, 1.0f);
    s.deserialize("invertDirection", invertDirection_);
    s.deserialize("forceAxisAlignment", forceAxisAlignment_);
    s.deserialize("velocityCurveType", velocityCurveType_);
    s.deserialize("velocityCurveDuration", velocityCurveDuration_);
    s.deserialize("targetVelocity", targetVelocity_);
    s.deserialize("velocityCurveFile", velocityCurveFile_);
    s.deserialize("velocityCurvePeriodic", velocityCurvePeriodic_);
}


const std::string FlowIndicatorDetection::loggerCat_("voreen.flowsimulation.FlowIndicatorDetection");

FlowIndicatorDetection::FlowIndicatorDetection()
    : Processor()
    , parameterInport_(Port::INPORT, "floParametrization.inport", "Flow Parametrization Input")
    , vesselGraphPort_(Port::INPORT, "vesselgraph.inport", "Vessel Graph")
    , volumePort_(Port::INPORT, "volume.inport", "Velocity Data Port (Optional)")
    , parameterOutport_(Port::OUTPORT, "flowParametrization.outport", "Flow Parametrization Output")
    , flowIndicatorTable_("flowIndicatorTable", "Flow Indicators", 3, Processor::VALID)
    , cloneFlowIndicator_("cloneFlowIndicator", "Clone Flow Indicator")
    , removeFlowIndicator_("removeFlowIndicator", "Remove Flow Indicator")
    , resetFlowIndicators_("resetFlowIndicators", "Reset Flow Indicators")
    , swapVelocityAndPressureBoundaries_("swapVelocityAndPressureBoundaries", "Swap Velocity and Pressure Boundaries")
    , resetTransformation_("resetTransformation", "Reset Grid Transformation")
    , angleThreshold_("angleThreshold", "Angle Threshold", 15, 0, 90, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
    , indicatorName_("indicatorName", "Name")
    , indicatorColor_("indicatorColor", "Color", tgt::vec4(1.0f))
    , centerlinePosition_("position", "Position", 0, 0, 0)
    , radius_("radius", "Radius (mm)", 1.0f, 0.0f, 10.0f)
    , length_("length", "Length (mm)", 0.0f, 0.0f, 10.0f)
    , relativeRadiusCorrection_("relativeRadiusCorrection", "Relative Radius Correction", 1.0f, 0.1f, 2.0f, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , invertDirection_("invertDirection", "Invert Direction", false)
    , forceAxisAlignment_("forceAxisAlignment", "Force Axis Alignment", false)
    , adaptTransformation_("adaptTransformation", "Align Grid to Indicator")
    , transformationMatrix_("transformationMatrix", "Transformation Matrix (Linking)", tgt::mat4::identity, tgt::mat4(-99999), tgt::mat4(99999))
    , indicatorType_("flowType", "Flow Type")
    , flowProfile_("flowProfile", "Flow Profile")
    , velocityCurveType_("velocityCurveType", "Velocity Curve Type")
    , velocityCurveDuration_("velocityCurveDuration", "Velocity Curve Duration (s)", 0.0f, 0.0f, 20.0f)
    , targetVelocity_("targetVelocity", "Target Velocity (m/s)", 0.0f, 0.0f, 10.0f)
    , velocityCurveFile_("velocityCurveFile", "Velocity Curve File", "Velocity Curve File", "", "*.csv", FileDialogProperty::OPEN_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT, FileDialogProperty::ALWAYS_OFF)
    , velocityCurvePeriodicity_("velocityCurvePeriodicity", "Repeat Velocity Curve periodically", false)
    , velocityCurveScale_("velocityCurveScale", "Scale Velocity Magnitude", 1.0f, 0.01f, 10.0f)
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

    addProperty(flowIndicatorTable_);
        flowIndicatorTable_.setColumnLabel(0, "ID");
        flowIndicatorTable_.setColumnLabel(1, "Name");
        flowIndicatorTable_.setColumnLabel(2, "Type");
        flowIndicatorTable_.setGroupID("detection");
        ON_CHANGE(flowIndicatorTable_, FlowIndicatorDetection, updateIndicatorUI);
    addProperty(cloneFlowIndicator_);
        cloneFlowIndicator_.setGroupID("detection");
        ON_CHANGE(cloneFlowIndicator_, FlowIndicatorDetection, onCloneFlowIndicator);
    addProperty(removeFlowIndicator_);
        removeFlowIndicator_.setGroupID("detection");
        ON_CHANGE(removeFlowIndicator_, FlowIndicatorDetection, onRemoveFlowIndicator);
    addProperty(resetFlowIndicators_);
        resetFlowIndicators_.setGroupID("detection");
        ON_CHANGE_LAMBDA(resetFlowIndicators_, [this] {
            detectFlowIndicators(true);
        });
    addProperty(swapVelocityAndPressureBoundaries_);
        swapVelocityAndPressureBoundaries_.setGroupID("detection");
        ON_CHANGE(swapVelocityAndPressureBoundaries_, FlowIndicatorDetection, swapVelocityAndPressureBoundaries);
    addProperty(resetTransformation_);
        resetTransformation_.setGroupID("detection");
        ON_CHANGE(resetTransformation_, FlowIndicatorDetection, resetTransformation);
    addProperty(angleThreshold_);
        angleThreshold_.setGroupID("detection");
    setPropertyGroupGuiName("detection", "Detection");

    addProperty(indicatorName_);
        indicatorName_.setInstantUpdate(false);
        indicatorName_.setGroupID("indicator");
        ON_CHANGE_LAMBDA(indicatorName_, [this] { onIndicatorConfigChange(false); });
        ON_CHANGE(indicatorName_, FlowIndicatorDetection, buildTable);
    addProperty(indicatorColor_);
        indicatorColor_.setGroupID("indicator");
        ON_CHANGE_LAMBDA(indicatorColor_, [this] { onIndicatorConfigChange(false); });
    addProperty(centerlinePosition_);
        centerlinePosition_.setGroupID("indicator");
        ON_CHANGE_LAMBDA(centerlinePosition_, [this] { onIndicatorConfigChange(true); });
    addProperty(invertDirection_);
        invertDirection_.setGroupID("indicator");
        ON_CHANGE_LAMBDA(invertDirection_, [this] { onIndicatorConfigChange(true); });
    addProperty(forceAxisAlignment_);
        forceAxisAlignment_.setGroupID("indicator");
        ON_CHANGE_LAMBDA(forceAxisAlignment_, [this] { onIndicatorConfigChange(true); });
    addProperty(transformationMatrix_);
        //transformationMatrix_.setGroupID("indicator");
        transformationMatrix_.setVisibleFlag(false);
    addProperty(adaptTransformation_);
        adaptTransformation_.setGroupID("indicator");
        ON_CHANGE_LAMBDA(adaptTransformation_, [this] { adaptTransformation(); });
    addProperty(radius_);
        radius_.setGroupID("indicator");
        ON_CHANGE(radius_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(length_);
        length_.setGroupID("indicator");
        ON_CHANGE(length_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(relativeRadiusCorrection_);
        relativeRadiusCorrection_.setGroupID("indicator");
        ON_CHANGE_LAMBDA(relativeRadiusCorrection_, [this] { onIndicatorConfigChange(true); });
    addProperty(indicatorType_);
        indicatorType_.addOption("candidate", "Candidate", FlowIndicatorType::FIT_CANDIDATE);
        indicatorType_.addOption("velocity", "Velocity Boundary", FlowIndicatorType::FIT_VELOCITY);
        indicatorType_.addOption("pressure", "Pressure Boundary", FlowIndicatorType::FIT_PRESSURE);
        indicatorType_.addOption("measure", "Flux Measure", FlowIndicatorType::FIT_MEASURE);
        indicatorType_.setGroupID("indicator");
        ON_CHANGE(indicatorType_, FlowIndicatorDetection, onIndicatorConfigChange);
        ON_CHANGE(indicatorType_, FlowIndicatorDetection, buildTable);
    setPropertyGroupGuiName("indicator", "Indicator");

    addProperty(flowProfile_);
        flowProfile_.addOption("poiseuille", "Poiseuille", FlowProfile::FP_POISEUILLE);
        flowProfile_.addOption("powerlaw", "Powerlaw", FlowProfile::FP_POWERLAW);
        //flowProfile_.addOption("constant", "constant", FlowProfile::FP_CONSTANT); // Not really useful.
        flowProfile_.addOption("volume", "Volume", FlowProfile::FP_VOLUME);
        flowProfile_.setGroupID("velocity");
        ON_CHANGE(flowProfile_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(velocityCurveType_);
        velocityCurveType_.addOption("constant", "Constant"); // get's selected automatically
        velocityCurveType_.addOption("linear", "Linear");
        velocityCurveType_.addOption("sinus", "Sinus");
        velocityCurveType_.addOption("heartBeat", "Heart Beat");
        velocityCurveType_.addOption("custom", "Custom");
        velocityCurveType_.setGroupID("velocity");
        ON_CHANGE(velocityCurveType_, FlowIndicatorDetection, onIndicatorConfigChange);
        ON_CHANGE_LAMBDA(velocityCurveType_, [this] { velocityCurveFile_.set(""); });
    addProperty(velocityCurveFile_);
        velocityCurveFile_.setGroupID("velocity");
        ON_CHANGE(velocityCurveFile_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(velocityCurveDuration_);
        velocityCurveDuration_.setGroupID("velocity");
        ON_CHANGE(velocityCurveDuration_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(targetVelocity_);
        targetVelocity_.setNumDecimals(3);
        targetVelocity_.setGroupID("velocity");
        ON_CHANGE(targetVelocity_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(velocityCurvePeriodicity_);
        velocityCurvePeriodicity_.setGroupID("velocity");
        ON_CHANGE(velocityCurvePeriodicity_, FlowIndicatorDetection, onIndicatorConfigChange);
    addProperty(velocityCurveScale_);
        velocityCurveScale_.setGroupID("velocity");
        ON_CHANGE(velocityCurveScale_, FlowIndicatorDetection, onIndicatorConfigChange);
    setPropertyGroupGuiName("velocity", "Velocity Boundary");

    // Setup property visibility.
    updateIndicatorUI();
}

void FlowIndicatorDetection::adjustPropertiesToInput() {
    if(!vesselGraphPort_.hasData())
        return;

    radius_.setMaxValue(tgt::length(vesselGraphPort_.getData()->getBounds().diagonal() / 2.0f));
    length_.setMaxValue(tgt::length(vesselGraphPort_.getData()->getBounds().diagonal() / 2.0f));
}

void FlowIndicatorDetection::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("vesselGraphHash", vesselGraphHash_);
    serializeVector<FlowIndicatorSerializable, FlowIndicator>(s, "flowIndicators", flowIndicators_);
    s.serialize("flowIndicatorSettings", flowIndicatorSettings_);
}

void FlowIndicatorDetection::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.deserialize("vesselGraphHash", vesselGraphHash_);
    deserializeVector<FlowIndicatorSerializable, FlowIndicator>(s, "flowIndicators", flowIndicators_);
    try {
        s.deserialize("flowIndicatorSettings", flowIndicatorSettings_);
    } catch(SerializationException&) {
        s.removeLastError();
        vesselGraphHash_.clear();
        flowIndicators_.clear();
    }

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
    FlowSimulationConfig* config = new FlowSimulationConfig(*parameterInport_.getData());

    // Set transformation matrix.
    if (config->getTransformationMatrix() != tgt::mat4::identity &&
        transformationMatrix_.get() != config->getTransformationMatrix()) {
        LWARNING("Overriding transformation matrix");
    }
    config->setTransformationMatrix(transformationMatrix_.get());

    for(const FlowIndicator& indicator : flowIndicators_) {
        // Candidates won't be added until an automatic or manual selection was performed.
        if(indicator.type_ != FIT_CANDIDATE) {
            config->addFlowIndicator(indicator);
        }
    }

    parameterOutport_.setData(config);
}

void FlowIndicatorDetection::updateIndicatorUI() {
    bool validSelection = flowIndicatorTable_.getNumRows() > 0 && flowIndicatorTable_.getSelectedRowIndex() >= 0;
    bool isRoleSwapped = false;
    if(validSelection) {
        triggertBySelection_ = true;

        auto index = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
        const FlowIndicator& indicator = flowIndicators_.at(index);
        const FlowIndicatorSettings& settings = flowIndicatorSettings_.at(index);

        if(vesselGraphPort_.hasData()) {
            const VesselGraph* vesselGraph = vesselGraphPort_.getData();
            centerlinePosition_.setMinValue(0);
            centerlinePosition_.setMaxValue(vesselGraph->getEdge(settings.edgeId_).getVoxels().size());
        }
        else {
            centerlinePosition_.setMinValue(settings.centerlinePosition_);
            centerlinePosition_.setMaxValue(settings.centerlinePosition_);
        }

        indicatorName_.set(indicator.name_);
        indicatorColor_.set(indicator.color_);
        centerlinePosition_.set(settings.centerlinePosition_);
        radius_.set(indicator.radius_);
        relativeRadiusCorrection_.set(settings.relativeRadiusCorrection_);
        length_.set(indicator.length_);
        invertDirection_.set(settings.invertDirection_);
        forceAxisAlignment_.set(settings.forceAxisAlignment_);
        indicatorType_.selectByValue(indicator.type_);
        isRoleSwapped = indicator.roleSwapped_;

        flowProfile_.selectByValue(indicator.flowProfile_);
        velocityCurveType_.select(settings.velocityCurveType_);
        velocityCurveDuration_.set(settings.velocityCurveDuration_);
        targetVelocity_.set(settings.targetVelocity_);
        velocityCurveFile_.set(settings.velocityCurveFile_);
        velocityCurvePeriodicity_.set(settings.velocityCurvePeriodic_);
        velocityCurveScale_.set(settings.velocityCurveScale_);

        triggertBySelection_ = false;
    }

    bool settingsEditable = validSelection && vesselGraphPort_.hasData();
    cloneFlowIndicator_.setReadOnlyFlag(!settingsEditable);
    removeFlowIndicator_.setReadOnlyFlag(!settingsEditable);
    centerlinePosition_.setReadOnlyFlag(!settingsEditable);

    indicatorName_.setReadOnlyFlag(!validSelection);
    indicatorColor_.setReadOnlyFlag(!validSelection);
    radius_.setReadOnlyFlag(!validSelection);
    relativeRadiusCorrection_.setReadOnlyFlag(!validSelection);
    length_.setReadOnlyFlag(!validSelection);
    invertDirection_.setReadOnlyFlag(!validSelection);
    forceAxisAlignment_.setReadOnlyFlag(!validSelection);
    adaptTransformation_.setReadOnlyFlag(!validSelection);
    indicatorType_.setReadOnlyFlag(!validSelection);

    bool isVelocityBoundary = validSelection && indicatorType_.getValue() == FIT_VELOCITY;
    setPropertyGroupVisible("velocity", isVelocityBoundary);
    flowProfile_.setReadOnlyFlag(!isVelocityBoundary);
    velocityCurveType_.setReadOnlyFlag(!validSelection);

    bool isVolumeCondition = flowProfile_.getValue() == FP_VOLUME;
    bool isSimpleCurveType = isVelocityBoundary && velocityCurveType_.get() != "custom" && velocityCurveType_.get() != "heartBeat";

    velocityCurveDuration_.setReadOnlyFlag(!isSimpleCurveType);
    targetVelocity_.setReadOnlyFlag(!isSimpleCurveType || isVolumeCondition);
    velocityCurveFile_.setReadOnlyFlag(velocityCurveType_.get() != "custom");
    velocityCurvePeriodicity_.setReadOnlyFlag(!isVelocityBoundary);
    velocityCurveScale_.setReadOnlyFlag(!isVelocityBoundary || !isRoleSwapped);
}

void FlowIndicatorDetection::onIndicatorConfigChange(bool needReinitialization) {

    // Ignore calls being triggered during selection callback.
    if(triggertBySelection_)
        return;

    // Handle selection.
    auto indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {

        FlowIndicatorSettings& settings = flowIndicatorSettings_[indicatorIdx];
        settings.centerlinePosition_ = centerlinePosition_.get();
        settings.relativeRadiusCorrection_ = relativeRadiusCorrection_.get();
        settings.invertDirection_ = invertDirection_.get();
        settings.forceAxisAlignment_ = forceAxisAlignment_.get();
        settings.velocityCurveType_ = velocityCurveType_.get();
        settings.velocityCurveDuration_ = velocityCurveDuration_.get();
        settings.targetVelocity_ = targetVelocity_.get();
        settings.velocityCurveFile_ = velocityCurveFile_.get();
        settings.velocityCurvePeriodic_ = velocityCurvePeriodicity_.get();
        settings.velocityCurveScale_ = velocityCurveScale_.get();

        FlowIndicator& indicator = flowIndicators_[indicatorIdx];
        indicator.type_ = indicatorType_.getValue();
        indicator.name_ = indicatorName_.get();
        indicator.color_ = indicatorColor_.get();
        indicator.length_ = length_.get();

        if(needReinitialization) {
            // initialization will overwrite some fields that we don't want to lose.
            // Hence, we back them up.
            const auto backup = indicator;

            indicator = initializeIndicator(settings);
            indicator.type_ = backup.type_;
            indicator.id_ = backup.id_;
            indicator.name_ = backup.name_;
            indicator.color_ = backup.color_;
            indicator.length_ = backup.length_;
            indicator.roleSwapped_ = backup.roleSwapped_;
        }
        else {
            indicator.radius_ = radius_.get();
            indicator.flowProfile_ = flowProfile_.getValue();
            indicator.velocityCurve_ = createCurveFromSettings(settings);
        }

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
            FlowIndicatorSettings settings(node.getID(), node.getEdges().back().get().getID());

            // Initialize indicator according to those settings.
            FlowIndicator indicator = initializeIndicator(settings);

            if(indicator.type_ != FIT_INVALID) {
                flowIndicators_.push_back(indicator);
                flowIndicatorSettings_.push_back(settings);
            }
        }
    }

    buildTable();
    updateIndicatorUI();
}

FlowIndicatorType FlowIndicatorDetection::estimateType(const FlowIndicator& indicator, const tgt::vec3& velocity) const {
    if (velocity != tgt::vec3::zero) {
        tgt::vec3 v = tgt::normalize(velocity);
        float threshold = tgt::deg2rad(static_cast<float>(angleThreshold_.get()));
        float angle = std::acos(tgt::dot(v, indicator.normal_));
        if (angle < threshold) {
            return FlowIndicatorType::FIT_VELOCITY;
        }
        else if (tgt::PIf - angle < threshold) {
            return FlowIndicatorType::FIT_PRESSURE;
        }
    }

    return FlowIndicatorType::FIT_CANDIDATE;
}

VelocityCurve FlowIndicatorDetection::createCurveFromSettings(FlowIndicatorSettings& settings) {
    VelocityCurve velocityCurve;

    if (settings.velocityCurveType_ == "constant") {
        velocityCurve = VelocityCurve::createConstantCurve(settings.targetVelocity_);
    } else if (settings.velocityCurveType_ == "linear") {
        velocityCurve = VelocityCurve::createLinearCurve(settings.velocityCurveDuration_,
                                                                    settings.targetVelocity_);
    } else if (settings.velocityCurveType_ == "sinus") {
        velocityCurve = VelocityCurve::createSinusoidalCurve(settings.velocityCurveDuration_,
                                                                        settings.targetVelocity_);
    } else if (settings.velocityCurveType_ == "heartBeat") {
        velocityCurve = VelocityCurve::createHumanHeartBeat();
    } else if (settings.velocityCurveType_ == "custom") {
        if(!settings.velocityCurveFile_.empty()) {
            try {
                velocityCurve = VelocityCurve::createFromCSV(settings.velocityCurveFile_);
                settings.targetVelocity_ = velocityCurve.getMaxVelocity();
            }
            catch (VoreenException& e) {
                VoreenApplication::app()->showMessageBox("Failed loading Curve", e.what());
                settings.velocityCurveFile_.clear();
            }
        }
    } else {
        tgtAssert(false, "Unhandled velocity curve");
        return velocityCurve;
    }

    velocityCurve.setPeriodic(settings.velocityCurvePeriodic_);
    velocityCurve.setScale(settings.velocityCurveScale_);

    return velocityCurve;
}

void FlowIndicatorDetection::onCloneFlowIndicator() {
    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {
        FlowIndicatorSettings& settings = flowIndicatorSettings_[indicatorIdx];
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

FlowIndicator FlowIndicatorDetection::initializeIndicator(FlowIndicatorSettings& settings) {

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

    // Default is candidate. However, we set initial values as if it was a velocity inlet,
    // since it might be classified as one later on.
    FlowIndicator indicator;
    indicator.id_ = flowIndicators_.size() + FlowSimulationConfig::getFlowIndicatorIdOffset();
    indicator.name_ = "Indicator " + std::to_string(indicator.id_);
    indicator.type_ = FlowIndicatorType::FIT_CANDIDATE;
    indicator.flowProfile_ = FlowProfile::FP_POISEUILLE;

    indicator.center_ = ref->pos_;
    indicator.normal_ = tgt::normalize(back->pos_ - front->pos_);
    indicator.radius_ = settings.relativeRadiusCorrection_ * radius;

    if(settings.forceAxisAlignment_) {
        float angleX = std::acos(std::abs(indicator.normal_.x));
        float angleY = std::acos(std::abs(indicator.normal_.y));
        float angleZ = std::acos(std::abs(indicator.normal_.z));

        if(angleX < angleY) {
            if(angleX < angleZ) {
                float sign = indicator.normal_.x < 0.0f ? -1.0f : 1.0f;
                indicator.normal_ = tgt::vec3(sign, 0.0f, 0.0f);
            }
            else {
                float sign = indicator.normal_.z < 0.0f ? -1.0f : 1.0f;
                indicator.normal_ = tgt::vec3(0.0f, 0.0f, sign);
            }
        }
        else {
            if(angleY < angleZ) {
                float sign = indicator.normal_.y < 0.0f ? -1.0f : 1.0f;
                indicator.normal_ = tgt::vec3(0.0f, sign, 0.0f);
            }
            else {
                float sign = indicator.normal_.z < 0.0f ? -1.0f : 1.0f;
                indicator.normal_ = tgt::vec3(0.0f, 0.0f, sign);
            }
        }
    }

    if(settings.invertDirection_) {
        indicator.normal_ *= -1.0f;
    }

    // Estimate velocity(direction) and therefore type.
    std::vector<tgt::vec3> samples = utils::sampleDisk(volumePort_.getData(), indicator.center_, indicator.normal_, indicator.radius_);

    tgt::vec3 accumDirection = tgt::vec3::zero;
    float maxMagnitudeSq = 0.0f;
    for(const auto& sample : samples) {
        accumDirection += sample;
        maxMagnitudeSq = std::max(maxMagnitudeSq, tgt::lengthSq(sample));
    }
    indicator.type_ = estimateType(indicator, accumDirection);

    // Pressure outlets point outwards.
    if (indicator.type_ == FIT_PRESSURE) {
        indicator.normal_ *= -1.0f;
    }

    // Setup velocity curve.
    settings.targetVelocity_ = std::sqrt(maxMagnitudeSq);
    indicator.velocityCurve_ = createCurveFromSettings(settings);

    return indicator;
}

void FlowIndicatorDetection::buildTable() {
    int selectedIndex = flowIndicatorTable_.getSelectedRowIndex();
    flowIndicatorTable_.reset();

    for(const FlowIndicator& indicator : flowIndicators_) {
        std::vector<std::string> row(flowIndicatorTable_.getNumColumns());
        row[0] = std::to_string(indicator.id_);
        row[1] = indicator.name_;
        row[2] = indicator.type_ == FlowIndicatorType::FIT_VELOCITY ? "Velocity" :
                 (indicator.type_ == FlowIndicatorType::FIT_PRESSURE ? "Pressure" :
                  (indicator.type_ == FlowIndicatorType::FIT_MEASURE ? "Measure" : "Candidate"));
        flowIndicatorTable_.addRow(row);
    }

    selectedIndex = std::min(selectedIndex, flowIndicatorTable_.getNumRows()-1);
    flowIndicatorTable_.setSelectedRowIndex(selectedIndex);
}

void FlowIndicatorDetection::swapVelocityAndPressureBoundaries() {

    std::vector<FlowIndicator> pressureBoundaries;
    std::vector<FlowIndicator> velocityBoundaries;
    for (const FlowIndicator& indicator: flowIndicators_) {
        if(indicator.type_ == FIT_PRESSURE) {
            pressureBoundaries.push_back(indicator);
        }
        else if(indicator.type_ == FIT_VELOCITY)  {
            velocityBoundaries.push_back(indicator);
        }
    }

    if(velocityBoundaries.size() != 1 && pressureBoundaries.size() != 1) {
        LERROR("Number of Velocity : Pressure Boundaries must be 1:n or n:1");
        return;
    }

    for (FlowIndicator& indicator: flowIndicators_) {
        if (indicator.type_ == FIT_PRESSURE) {
            indicator.type_ = FIT_VELOCITY;
            indicator.roleSwapped_ = !indicator.roleSwapped_;

            if(indicator.roleSwapped_ && velocityBoundaries.size() == 1) {
                auto refIndicator = velocityBoundaries.front();

                indicator.velocityCurve_ = refIndicator.velocityCurve_;

                // TODO: very bad heuristic. Should use a different one here or none at all.
                float areaRatio = std::sqrt(indicator.radius_ / refIndicator.radius_);
                indicator.velocityCurve_.setScale(areaRatio);
            }
        }
        else if (indicator.type_ == FIT_VELOCITY) {
            indicator.type_ = FIT_PRESSURE;
            indicator.roleSwapped_ = !indicator.roleSwapped_;
        }
        else {
            indicator.roleSwapped_ = false;
        }
    }

    buildTable();
    updateIndicatorUI();
}

void FlowIndicatorDetection::adaptTransformation() {
    size_t indicatorIdx = static_cast<size_t>(flowIndicatorTable_.getSelectedRowIndex());
    if(flowIndicatorTable_.getNumRows() > 0 && indicatorIdx < flowIndicators_.size()) {
        auto& indicator = flowIndicators_.at(indicatorIdx);

        auto rotation = utils::createTransformationMatrix(indicator.center_, indicator.normal_).getRotationalPart();
        tgt::mat4 transformationMatrix;
        rotation.invert(transformationMatrix);
        transformationMatrix_.set(transformationMatrix);
    }
}

void FlowIndicatorDetection::resetTransformation() {
    transformationMatrix_.set(tgt::mat4::identity);
}

}   // namespace
