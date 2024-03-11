/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "flowsimulationconfig.h"

#include "../utils/serializationhelper.h"

namespace voreen {

const int FlowSimulationConfig::VERSION = 1;
const size_t FlowSimulationConfig::ALL_PARAMETER_SETS = static_cast<size_t>(-1);

FlowSimulationConfig::FlowSimulationConfig(const std::string& name)
    : name_(name)
    , simulationTime_(0.0f)
    , numTimeSteps_(0)
    , outputResolution_(0)
    , flowFeatures_(FF_NONE)
    , transformation_(tgt::mat4::identity)
    , geometryIsMesh_(true)
{
}

FlowSimulationConfig::FlowSimulationConfig(const FlowSimulationConfig& origin)
    : name_(origin.name_)
    , simulationTime_(origin.simulationTime_)
    , numTimeSteps_(origin.numTimeSteps_)
    , outputResolution_(origin.outputResolution_)
    , outputFileFormat_(origin.outputFileFormat_)
    , flowFeatures_(origin.flowFeatures_)
    , flowIndicators_(origin.flowIndicators_)
    , flowParameters_(origin.flowParameters_)
    , transformation_(origin.transformation_)
    , geometryIsMesh_(true)
{
}

const std::string& FlowSimulationConfig::getName() const {
    return name_;
}

float FlowSimulationConfig::getSimulationTime() const {
    return simulationTime_;
}

void FlowSimulationConfig::setSimulationTime(float simulationTime) {
    notifyPendingDataInvalidation();
    simulationTime_ = simulationTime;
}

int FlowSimulationConfig::getNumTimeSteps() const {
    return numTimeSteps_;
}

void FlowSimulationConfig::setNumTimeSteps(int numTimeSteps) {
    notifyPendingDataInvalidation();
    numTimeSteps_ = numTimeSteps;
}

int FlowSimulationConfig::getOutputResolution() const {
    return outputResolution_;
}

void FlowSimulationConfig::setOutputResolution(int outputResolution) {
    notifyPendingDataInvalidation();
    outputResolution_ = outputResolution;
}

const std::string& FlowSimulationConfig::getOutputFileFormat() const {
    return outputFileFormat_;
}

void FlowSimulationConfig::setOutputFileFormat(const std::string& format) {
    notifyPendingDataInvalidation();
    outputFileFormat_ = format;
}

int FlowSimulationConfig::getFlowFeatures() const {
    return flowFeatures_;
}

void FlowSimulationConfig::setFlowFeatures(int flowFeatures) {
    notifyPendingDataInvalidation();
    flowFeatures_ = flowFeatures;
}

tgt::mat4 FlowSimulationConfig::getTransformationMatrix() const {
    return transformation_;
}

tgt::mat4 FlowSimulationConfig::getInvertedTransformationMatrix() const {
    tgt::mat4 inverted = tgt::mat4::identity;
    if (!transformation_.invert(inverted)) {
        LERRORC("FlowSimulationConfig", "Could not invert transformation matrix");
    }
    return inverted;
}

void FlowSimulationConfig::setTransformationMatrix(tgt::mat4 transformation) {
    transformation_ = transformation;
}

void FlowSimulationConfig::setGeometryFiles(const std::map<float, std::string>& geometryFiles, bool isMesh) {
    notifyPendingDataInvalidation();
    geometryFiles_ = geometryFiles;
    geometryIsMesh_ = isMesh;
}
const std::map<float, std::string>& FlowSimulationConfig::getGeometryFiles() const {
    return geometryFiles_;
}

void FlowSimulationConfig::setGeometryIsMesh(bool isMesh) {
    notifyPendingDataInvalidation();
    geometryIsMesh_ = isMesh;
}
bool FlowSimulationConfig::isGeometryMesh() const {
    return geometryIsMesh_;
}

void FlowSimulationConfig::setMeasuredDataFiles(const std::map<float, std::string>& measuredDataFiles) {
    notifyPendingDataInvalidation();
    measuredDataFiles_ = measuredDataFiles;
}
const std::map<float, std::string>& FlowSimulationConfig::getMeasuredDataFiles() const {
    return measuredDataFiles_;
}

void FlowSimulationConfig::addFlowIndicator(const FlowIndicator& flowIndicator) {
    notifyPendingDataInvalidation();

    // Copy the indicator and set its id.
    FlowIndicator indicator = flowIndicator;
    indicator.id_ = generateIndicatorId();
    flowIndicators_.push_back(indicator);
}

std::vector<FlowIndicator> FlowSimulationConfig::getFlowIndicators(bool transformed) const {
    std::vector<FlowIndicator> indicators = flowIndicators_;
    if (transformed) {
        auto normalMatrix = transpose(getInvertedTransformationMatrix());
        for (auto& indicator : indicators) {
            indicator.center_ = transformation_ * indicator.center_;
            indicator.normal_ = normalMatrix * indicator.normal_;
        }
    }
    return indicators;
}

void FlowSimulationConfig::addFlowParameterSet(const Parameters& parameters) {
    notifyPendingDataInvalidation();
    flowParameters_.push_back(parameters);
}

const std::vector<Parameters>& FlowSimulationConfig::getFlowParameterSets() const {
    return flowParameters_;
}

bool FlowSimulationConfig::empty() const {
    return flowParameters_.empty();
}
size_t FlowSimulationConfig::size() const {
    return flowParameters_.size();
}
const Parameters& FlowSimulationConfig::at(size_t index) const {
    return flowParameters_.at(index);
}

std::string FlowSimulationConfig::toJSONString(size_t param) const {
    std::stringstream stream;
    JsonSerializer json;
    Serializer s(json);
    serializeInternal(s, param);
    json.write(stream, true, false);
    return stream.str();
}

std::string FlowSimulationConfig::toXMLString(size_t param) const {
    std::stringstream stream;
    XmlSerializer xml;
    Serializer s(xml);
    serializeInternal(s, param);
    xml.write(stream);
    return stream.str();
}

void FlowSimulationConfig::serialize(Serializer& s) const {
    serializeInternal(s, ALL_PARAMETER_SETS);
}

void FlowSimulationConfig::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("simulationTime", simulationTime_);
    s.deserialize("numTimeSteps", numTimeSteps_);
    s.deserialize("outputResolution", outputResolution_);
    s.deserialize("outputFileFormat", outputFileFormat_);
    s.deserialize("flowFeatures", flowFeatures_);
    s.deserialize("transformation", transformation_);
    s.optionalDeserialize("geometryFiles", geometryFiles_, decltype(geometryFiles_)());
    s.optionalDeserialize("geometryIsMesh", geometryIsMesh_, true);
    s.optionalDeserialize("measuredDataFiles", measuredDataFiles_, decltype(measuredDataFiles_)());
    deserializeVector<FlowIndicatorSerializable, FlowIndicator>(s, "flowIndicators", flowIndicators_);
    deserializeVector<ParametersSerializable, Parameters>(s, "flowParameters", flowParameters_);
}

void FlowSimulationConfig::serializeInternal(Serializer& s, size_t param) const {
    s.serialize("name", name_);
    s.serialize("simulationTime", simulationTime_);
    s.serialize("numTimeSteps", numTimeSteps_);
    s.serialize("outputResolution", outputResolution_);
    s.serialize("outputFileFormat", outputFileFormat_);
    s.serialize("flowFeatures", flowFeatures_);
    s.serialize("transformation", transformation_);
    s.serialize("geometryFiles", geometryFiles_);
    s.serialize("geometryIsMesh", geometryIsMesh_);
    s.serialize("measuredDataFiles", measuredDataFiles_);
    serializeVector<FlowIndicatorSerializable, FlowIndicator>(s, "flowIndicators", flowIndicators_);
    if(param == ALL_PARAMETER_SETS) {
        serializeVector<ParametersSerializable, Parameters>(s, "flowParametrizations", flowParameters_);
    }
    else {
        s.serialize("flowParameters", ParametersSerializable(flowParameters_[param]));
    }
}

int FlowSimulationConfig::generateIndicatorId() const {
    return static_cast<int>(flowIndicators_.size()) + getFlowIndicatorIdOffset();
}

}   // namespace
