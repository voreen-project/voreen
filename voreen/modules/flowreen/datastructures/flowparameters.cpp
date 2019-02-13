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

FlowParameters::FlowParameters(const std::string& name)
    : name_(name)
    , characteristicLength_(0.0f)
    , viscosity_(0.0f)
    , density_(0.0f)
    , bouzidi_(false)
{
}

const std::string& FlowParameters::getName() const {
    return name_;
}

float FlowParameters::getCharacteristicLength() const {
    return characteristicLength_;
}

void FlowParameters::setCharacteristicLength(float characteristicLength) {
    characteristicLength_ = characteristicLength;
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

bool FlowParameters::getBouzidi() const {
    return bouzidi_;
}

void FlowParameters::setBouzidi(bool bouzidi) {
    bouzidi_ = bouzidi;
}

void FlowParameters::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("characteristicLength", characteristicLength_);
    s.serialize("viscosity", viscosity_);
    s.serialize("density", density_);
    s.serialize("bouzidi", bouzidi_);
}

void FlowParameters::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("characteristicLength", characteristicLength_);
    s.deserialize("viscosity", viscosity_);
    s.deserialize("density", density_);
    s.deserialize("bouzidi", bouzidi_);
}


const int FlowParametrizationList::VERSION = 1;

FlowParametrizationList::FlowParametrizationList(const std::string& name)
    : name_(name)
    , simulationTime_(0.0f)
    , temporalResolution_(0.0f)
{
}

FlowParametrizationList::FlowParametrizationList(const FlowParametrizationList& origin)
    : name_(origin.name_)
    , simulationTime_(origin.simulationTime_)
    , temporalResolution_(origin.temporalResolution_)
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

float FlowParametrizationList::getTemporalResolution() const {
    return temporalResolution_;
}

void FlowParametrizationList::setTemporalResolution(float temporalResolution) {
    notifyPendingDataInvalidation();
    temporalResolution_ = temporalResolution;
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

std::string FlowParametrizationList::toCSVString() const {
    std::stringstream output;

    output << VERSION;
    output << ", " << simulationTime_;
    output << ", " << temporalResolution_;

    output << ", " << flowParametrizations_.size();
    for(const FlowParameters& flowParameters : flowParametrizations_) {
        output << ", " << flowParameters.getCharacteristicLength();
        output << ", " << flowParameters.getViscosity();
        output << ", "<< flowParameters.getDensity();
        output << ", "<< flowParameters.getBouzidi();
    }

    output << ", " << flowIndicators_.size();
    for(const FlowIndicator& flowIndicator : flowIndicators_) {
        output << ", " << flowIndicator.direction_;
        output << ", " << flowIndicator.center_.x << ", " << flowIndicator.center_.y << ", " << flowIndicator.center_.z;
        output << ", " << flowIndicator.normal_.x << ", " << flowIndicator.normal_.y << ", " << flowIndicator.normal_.z;
        output << ", " << flowIndicator.radius_;
    }

    return output.str();
}

std::string FlowParametrizationList::toJSONString() const {
    std::stringstream stream;
    JsonSerializer json;
    Serializer s(json);
    serialize(s);
    json.write(stream, true, false);
    return stream.str();
}

void FlowParametrizationList::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("simulationTime", simulationTime_);
    s.serialize("temporalResolution", temporalResolution_);
    s.serializeBinaryBlob("flowIndicators", flowIndicators_);
    s.serialize("flowParametrizations", flowParametrizations_);
}

void FlowParametrizationList::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("simulationTime", simulationTime_);
    s.deserialize("temporalResolution", temporalResolution_);
    s.deserializeBinaryBlob("flowIndicators", flowIndicators_);
    s.deserialize("flowParametrizations", flowParametrizations_,
                  XmlSerializationConstants::ITEMNODE,
                  std::function<FlowParameters()>([]{ return FlowParameters(""); }));
}

}   // namespace
