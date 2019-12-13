/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "flowparameters.h"

#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"

#include <sstream>

namespace voreen {

FlowIndicator::FlowIndicator()
    : type_(FIT_INVALID)
    , id_(-1)
    , center_(tgt::vec3::zero)
    , normal_(tgt::vec3::zero)
    , radius_(0.0f)
    , flowProfile_(FP_NONE)
    , startPhaseFunction_(FSP_NONE)
    , startPhaseDuration_(0.0f)
    , targetVelocity_(0.0f)
    , selected_(false)
{
}

void FlowIndicator::serialize(Serializer& s) const {
    s.serialize("type_", type_); // "type" is a reserved xml keyword.
    s.serialize("id_", id_); // "id" is a reserved xml keyword.
    s.serialize("center", center_);
    s.serialize("normal", normal_);
    s.serialize("radius", radius_);
    s.serialize("flowProfile", flowProfile_);
    s.serialize("startPhaseFunction", startPhaseFunction_);
    s.serialize("startPhaseDuration", startPhaseDuration_);
    s.serialize("targetVelocity", targetVelocity_);
}

void FlowIndicator::deserialize(Deserializer& s) {
    int type = FIT_CANDIDATE;
    s.deserialize("type_", type);
    type_ = static_cast<FlowIndicatorType>(type);
    s.deserialize("id_", id_);
    s.deserialize("center", center_);
    s.deserialize("normal", normal_);
    s.deserialize("radius", radius_);
    int profile = FP_POISEUILLE;
    s.deserialize("flowProfile", profile);
    flowProfile_ = static_cast<FlowProfile>(profile);
    int function = FSP_NONE;
    s.deserialize("startPhaseFunction", function);
    startPhaseFunction_ = static_cast<FlowStartPhase>(function);
    s.deserialize("startPhaseDuration", startPhaseDuration_);
    s.deserialize("targetVelocity", targetVelocity_);
}

FlowParameterSet::FlowParameterSet()
    : FlowParameterSet("")
{
}

FlowParameterSet::FlowParameterSet(const std::string& name)
    : name_(name)
    , spatialResolution_(0)
    , temporalResolution_(0.0f)
    , characteristicLength_(0.0f)
    , characteristicVelocity_(0.0f)
    , viscosity_(0.0f)
    , density_(0.0f)
    , smagorinskyConstant_(0.0f)
    , bouzidi_(false)
{
}

const std::string& FlowParameterSet::getName() const {
    return name_;
}

int FlowParameterSet::getSpatialResolution() const {
    return spatialResolution_;
}

void FlowParameterSet::setSpatialResolution(int spatialResolution) {
    spatialResolution_ = spatialResolution;
}

float FlowParameterSet::getTemporalResolution() const {
    return temporalResolution_;
}

void FlowParameterSet::setTemporalResolution(float temporalResolution) {
    temporalResolution_ = temporalResolution;
}

float FlowParameterSet::getCharacteristicLength() const {
    return characteristicLength_;
}

void FlowParameterSet::setCharacteristicLength(float characteristicLength) {
    characteristicLength_ = characteristicLength;
}

float FlowParameterSet::getCharacteristicVelocity() const {
    return characteristicVelocity_;
}

void FlowParameterSet::setCharacteristicVelocity(float characteristicVelocity) {
    characteristicVelocity_ = characteristicVelocity;
}

float FlowParameterSet::getViscosity() const {
    return viscosity_;
}

void FlowParameterSet::setViscosity(float viscosity) {
    viscosity_ = viscosity;
}

float FlowParameterSet::getDensity() const {
    return density_;
}

void FlowParameterSet::setDensity(float density) {
    density_ = density;
}

float FlowParameterSet::getSmagorinskyConstant() const {
    return smagorinskyConstant_;
}

void FlowParameterSet::setSmagorinskyConstant(float smagorinskyConstant) {
    smagorinskyConstant_ = smagorinskyConstant;
}

bool FlowParameterSet::getBouzidi() const {
    return bouzidi_;
}

void FlowParameterSet::setBouzidi(bool bouzidi) {
    bouzidi_ = bouzidi;
}

void FlowParameterSet::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("spatialResolution", spatialResolution_);
    s.serialize("temporalResolution", temporalResolution_);
    s.serialize("characteristicLength", characteristicLength_);
    s.serialize("characteristicVelocity", characteristicVelocity_);
    s.serialize("viscosity", viscosity_);
    s.serialize("density", density_);
    s.serialize("smagorinskyConstant", smagorinskyConstant_);
    s.serialize("bouzidi", bouzidi_);
}

void FlowParameterSet::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("spatialResolution", spatialResolution_);
    s.deserialize("temporalResolution", temporalResolution_);
    s.deserialize("characteristicLength", characteristicLength_);
    s.deserialize("characteristicVelocity", characteristicVelocity_);
    s.deserialize("viscosity", viscosity_);
    s.deserialize("density", density_);
    s.deserialize("smagorinskyConstant", smagorinskyConstant_);
    s.deserialize("bouzidi", bouzidi_);
}


const int FlowParameterSetEnsemble::VERSION = 1;
const int FlowParameterSetEnsemble::FLOW_INDICATOR_ID_OFFSET = 3;
const size_t FlowParameterSetEnsemble::ALL_PARAMETER_SETS = static_cast<size_t>(-1);

int FlowParameterSetEnsemble::getFlowIndicatorIdOffset() {
    // TODO: The offset is determined by the number of reserved material IDs of the
    //  simulation framework. OpenLB has an default empty (0) and fluid(1) and wall (2) material.
    //  Thus, the offset is set to 3. This should be set explicitly by something like:
    //  SimulationFramework::RESERVED_MATERIAL_IDS
    return FLOW_INDICATOR_ID_OFFSET;
}

FlowParameterSetEnsemble::FlowParameterSetEnsemble(const std::string& name)
    : name_(name)
    , simulationTime_(0.0f)
    , numTimeSteps_(0)
    , outputResolution_(0)
    , flowFeatures_(FF_NONE)
{
}

FlowParameterSetEnsemble::FlowParameterSetEnsemble(const FlowParameterSetEnsemble& origin)
    : name_(origin.name_)
    , simulationTime_(origin.simulationTime_)
    , numTimeSteps_(origin.numTimeSteps_)
    , outputResolution_(origin.outputResolution_)
    , flowFeatures_(origin.flowFeatures_)
    , flowIndicators_(origin.flowIndicators_)
    , flowParameterSets(origin.flowParameterSets)
{
}

const std::string& FlowParameterSetEnsemble::getName() const {
    return name_;
}

float FlowParameterSetEnsemble::getSimulationTime() const {
    return simulationTime_;
}

void FlowParameterSetEnsemble::setSimulationTime(float simulationTime) {
    notifyPendingDataInvalidation();
    simulationTime_ = simulationTime;
}

int FlowParameterSetEnsemble::getNumTimeSteps() const {
    return numTimeSteps_;
}

void FlowParameterSetEnsemble::setNumTimeSteps(int numTimeSteps) {
    notifyPendingDataInvalidation();
    numTimeSteps_ = numTimeSteps;
}

int FlowParameterSetEnsemble::getOutputResolution() const {
    return outputResolution_;
}

void FlowParameterSetEnsemble::setOutputResolution(int outputResolution) {
    notifyPendingDataInvalidation();
    outputResolution_ = outputResolution;
}

int FlowParameterSetEnsemble::getFlowFeatures() const {
    return flowFeatures_;
}

void FlowParameterSetEnsemble::setFlowFeatures(int flowFeatures) {
    notifyPendingDataInvalidation();
    flowFeatures_ = flowFeatures;
}

void FlowParameterSetEnsemble::addFlowIndicator(const FlowIndicator& flowIndicator) {
    notifyPendingDataInvalidation();

    // Copy the indicator and set its id.
    FlowIndicator indicator = flowIndicator;
    indicator.id_ = generateIndicatorId();
    flowIndicators_.push_back(indicator);
}

const std::vector<FlowIndicator>& FlowParameterSetEnsemble::getFlowIndicators() const {
    return flowIndicators_;
}

void FlowParameterSetEnsemble::addFlowParameterSet(const voreen::FlowParameterSet& parameters) {
    notifyPendingDataInvalidation();
    flowParameterSets.push_back(parameters);
}

const std::vector<FlowParameterSet>& FlowParameterSetEnsemble::getFlowParameterSets() const {
    return flowParameterSets;
}

bool FlowParameterSetEnsemble::empty() const {
    return flowParameterSets.empty();
}
size_t FlowParameterSetEnsemble::size() const {
    return flowParameterSets.size();
}
const FlowParameterSet& FlowParameterSetEnsemble::at(size_t index) const {
    return flowParameterSets.at(index);
}

std::string FlowParameterSetEnsemble::toJSONString(size_t param) const {
    std::stringstream stream;
    JsonSerializer json;
    Serializer s(json);
    serializeInternal(s, param);
    json.write(stream, true, false);
    return stream.str();
}

std::string FlowParameterSetEnsemble::toXMLString(size_t param) const {
    std::stringstream stream;
    XmlSerializer xml;
    Serializer s(xml);
    serializeInternal(s, param);
    xml.write(stream);
    return stream.str();
}

void FlowParameterSetEnsemble::serialize(Serializer& s) const {
    serializeInternal(s, ALL_PARAMETER_SETS);
}

void FlowParameterSetEnsemble::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("simulationTime", simulationTime_);
    s.deserialize("numTimeSteps", numTimeSteps_);
    s.deserialize("outputResolution", outputResolution_);
    s.deserialize("flowFeatures", flowFeatures_);
    s.deserialize("flowIndicators", flowIndicators_);
    s.deserialize("flowParametrizations", flowParameterSets);
}

void FlowParameterSetEnsemble::serializeInternal(Serializer& s, size_t param) const {
    s.serialize("name", name_);
    s.serialize("simulationTime", simulationTime_);
    s.serialize("numTimeSteps", numTimeSteps_);
    s.serialize("outputResolution", outputResolution_);
    s.serialize("flowFeatures", flowFeatures_);
    s.serialize("flowIndicators", flowIndicators_);
    if(param == ALL_PARAMETER_SETS) {
        s.serialize("flowParametrizations", flowParameterSets);
    }
    else {
        s.serialize("flowParameters", flowParameterSets[param]);
    }
}

int FlowParameterSetEnsemble::generateIndicatorId() const {
    return static_cast<int>(flowIndicators_.size()) + getFlowIndicatorIdOffset();
}

}   // namespace
