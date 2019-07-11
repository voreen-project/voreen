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

#include "flowparameters.h"

#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"

#include <sstream>

namespace voreen {

FlowIndicator::FlowIndicator()
    : direction_(FD_NONE)
    , startPhaseFunction_(FF_NONE)
    , startPhaseDuration_(0.0f)
    , center_(tgt::vec3::zero)
    , normal_(tgt::vec3::zero)
    , radius_(0.0f)
{
}

void FlowIndicator::serialize(Serializer& s) const {
    s.serialize("direction", direction_);
    s.serialize("startPhaseFunction", startPhaseFunction_);
    s.serialize("startPhaseDuration", startPhaseDuration_);
    s.serialize("center", center_);
    s.serialize("normal", normal_);
    s.serialize("radius", radius_);
}

void FlowIndicator::deserialize(Deserializer& s) {
    int direction = FD_NONE;
    s.deserialize("direction", direction);
    direction_ = static_cast<FlowDirection>(direction);
    int function = FF_NONE;
    s.deserialize("function", function);
    startPhaseFunction_ = static_cast<FlowFunction>(function);
    s.deserialize("startPhaseDuration", startPhaseDuration_);
    s.deserialize("center", center_);
    s.deserialize("normal", normal_);
    s.deserialize("radius", radius_);
}

FlowParameters::FlowParameters()
    : FlowParameters("")
{
}

FlowParameters::FlowParameters(const std::string& name)
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

const std::string& FlowParameters::getName() const {
    return name_;
}

int FlowParameters::getSpatialResolution() const {
    return spatialResolution_;
}

void FlowParameters::setSpatialResolution(int spatialResolution) {
    spatialResolution_ = spatialResolution;
}

float FlowParameters::getTemporalResolution() const {
    return temporalResolution_;
}

void FlowParameters::setTemporalResolution(float temporalResolution) {
    temporalResolution_ = temporalResolution;
}

float FlowParameters::getCharacteristicLength() const {
    return characteristicLength_;
}

void FlowParameters::setCharacteristicLength(float characteristicLength) {
    characteristicLength_ = characteristicLength;
}

float FlowParameters::getCharacteristicVelocity() const {
    return characteristicVelocity_;
}

void FlowParameters::setCharacteristicVelocity(float characteristicVelocity) {
    characteristicVelocity_ = characteristicVelocity;
}

float FlowParameters::getViscosity() const {
    return viscosity_;
}

void FlowParameters::setViscosity(float viscosity) {
    viscosity_ = viscosity;
}

float FlowParameters::getDensity() const {
    return density_;
}

void FlowParameters::setDensity(float density) {
    density_ = density;
}

float FlowParameters::getSmagorinskyConstant() const {
    return smagorinskyConstant_;
}

void FlowParameters::setSmagorinskyConstant(float smagorinskyConstant) {
    smagorinskyConstant_ = smagorinskyConstant;
}

bool FlowParameters::getBouzidi() const {
    return bouzidi_;
}

void FlowParameters::setBouzidi(bool bouzidi) {
    bouzidi_ = bouzidi;
}

void FlowParameters::serialize(Serializer& s) const {
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

void FlowParameters::deserialize(Deserializer& s) {
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


const int FlowParametrizationList::VERSION = 1;
const size_t FlowParametrizationList::ALL_PARAMETRIZATIONS = static_cast<size_t>(-1);

FlowParametrizationList::FlowParametrizationList(const std::string& name)
    : name_(name)
    , simulationTime_(0.0f)
    , numTimeSteps_(0)
    , outputResolution_(0)
{
}

FlowParametrizationList::FlowParametrizationList(const FlowParametrizationList& origin)
    : name_(origin.name_)
    , simulationTime_(origin.simulationTime_)
    , numTimeSteps_(origin.numTimeSteps_)
    , outputResolution_(origin.outputResolution_)
    , flowIndicators_(origin.flowIndicators_)
    , flowParametrizations_(origin.flowParametrizations_)
{
}

const std::string& FlowParametrizationList::getName() const {
    return name_;
}

float FlowParametrizationList::getSimulationTime() const {
    return simulationTime_;
}

void FlowParametrizationList::setSimulationTime(float simulationTime) {
    notifyPendingDataInvalidation();
    simulationTime_ = simulationTime;
}

int FlowParametrizationList::getNumTimeSteps() const {
    return numTimeSteps_;
}

void FlowParametrizationList::setNumTimeSteps(int numTimeSteps) {
    notifyPendingDataInvalidation();
    numTimeSteps_ = numTimeSteps;
}

int FlowParametrizationList::getOutputResolution() const {
    return outputResolution_;
}

void FlowParametrizationList::setOutputResolution(int outputResolution) {
    notifyPendingDataInvalidation();
    outputResolution_ = outputResolution;
}

void FlowParametrizationList::setStartPhaseFunction(FlowFunction startPhaseFunction) {
    notifyPendingDataInvalidation();
    for(FlowIndicator& flowIndicator : flowIndicators_) {
        if(flowIndicator.direction_ == FD_IN) {
            flowIndicator.startPhaseFunction_ = startPhaseFunction;
        }
    }
}

void FlowParametrizationList::addFlowIndicator(const FlowIndicator& flowIndicator) {
    notifyPendingDataInvalidation();
    flowIndicators_.push_back(flowIndicator);
}

const std::vector<FlowIndicator>& FlowParametrizationList::getFlowIndicators() const {
    return flowIndicators_;
}

void FlowParametrizationList::addFlowParameters(const voreen::FlowParameters& flowParameters) {
    notifyPendingDataInvalidation();
    flowParametrizations_.push_back(flowParameters);
}

const std::vector<FlowParameters>& FlowParametrizationList::getFlowParametrizations() const {
    return flowParametrizations_;
}

bool FlowParametrizationList::empty() const {
    return flowParametrizations_.empty();
}
size_t FlowParametrizationList::size() const {
    return flowParametrizations_.size();
}
const FlowParameters& FlowParametrizationList::at(size_t index) const {
    return flowParametrizations_.at(index);
}

std::string FlowParametrizationList::toCSVString(size_t param) const {
    std::stringstream output;

    output << VERSION;
    output << ", " << simulationTime_;
    output << ", " << numTimeSteps_;
    output << ", " << outputResolution_;

    if(param == ALL_PARAMETRIZATIONS) {
        output << ", " << flowParametrizations_.size();
        for (const FlowParameters& flowParameters : flowParametrizations_) {
            output << ", " << flowParameters.getSpatialResolution();
            output << ", " << flowParameters.getTemporalResolution();
            output << ", " << flowParameters.getCharacteristicLength();
            output << ", " << flowParameters.getCharacteristicVelocity();
            output << ", " << flowParameters.getViscosity();
            output << ", " << flowParameters.getDensity();
            output << ", " << flowParameters.getSmagorinskyConstant();
            output << ", " << flowParameters.getBouzidi();
        }
    }
    else {
        output << ", " << 1;
        output << ", " << flowParametrizations_[param].getSpatialResolution();
        output << ", " << flowParametrizations_[param].getTemporalResolution();
        output << ", " << flowParametrizations_[param].getCharacteristicLength();
        output << ", " << flowParametrizations_[param].getCharacteristicVelocity();
        output << ", " << flowParametrizations_[param].getViscosity();
        output << ", " << flowParametrizations_[param].getDensity();
        output << ", " << flowParametrizations_[param].getSmagorinskyConstant();
        output << ", " << flowParametrizations_[param].getBouzidi();
    }

    output << ", " << flowIndicators_.size();
    for(const FlowIndicator& flowIndicator : flowIndicators_) {
        output << ", " << flowIndicator.direction_;
        output << ", " << flowIndicator.startPhaseFunction_;
        output << ", " << flowIndicator.startPhaseDuration_;
        output << ", " << flowIndicator.center_.x << ", " << flowIndicator.center_.y << ", " << flowIndicator.center_.z;
        output << ", " << flowIndicator.normal_.x << ", " << flowIndicator.normal_.y << ", " << flowIndicator.normal_.z;
        output << ", " << flowIndicator.radius_;
    }

    return output.str();
}

std::string FlowParametrizationList::toJSONString(size_t param) const {
    std::stringstream stream;
    JsonSerializer json;
    Serializer s(json);
    serializeInternal(s, param);
    json.write(stream, true, false);
    return stream.str();
}

std::string FlowParametrizationList::toXMLString(size_t param) const {
    std::stringstream stream;
    XmlSerializer xml;
    Serializer s(xml);
    serializeInternal(s, param);
    xml.write(stream);
    return stream.str();
}

void FlowParametrizationList::serialize(Serializer& s) const {
    serializeInternal(s, ALL_PARAMETRIZATIONS);
}

void FlowParametrizationList::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("simulationTime", simulationTime_);
    s.deserialize("numTimeSteps", numTimeSteps_);
    s.deserialize("outputResolution", outputResolution_);
    s.deserialize("flowIndicators", flowIndicators_);
    s.deserialize("flowParametrizations", flowParametrizations_);
}

void FlowParametrizationList::serializeInternal(Serializer& s, size_t param) const {
    s.serialize("name", name_);
    s.serialize("simulationTime", simulationTime_);
    s.serialize("numTimeSteps", numTimeSteps_);
    s.serialize("outputResolution", outputResolution_);
    s.serialize("flowIndicators", flowIndicators_);
    if(param == ALL_PARAMETRIZATIONS) {
        s.serialize("flowParametrizations", flowParametrizations_);
    }
    else {
        s.serialize("flowParameters", flowParametrizations_[param]);
    }
}

}   // namespace
