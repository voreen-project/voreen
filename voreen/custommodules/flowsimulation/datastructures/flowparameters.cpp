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

#include "flowparameters.h"

#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"

#include <sstream>
#include <fstream>

namespace voreen {

VelocityCurve::VelocityCurve()
    : periodic_(false)
{
    peakVelocities_[0.0f] = 0.0f;
}

float VelocityCurve::operator()(float t) const {

    if(isPeriodic()) {
        float begin = peakVelocities_.begin()->first;
        float end = peakVelocities_.rbegin()->first;
        t = std::fmod(t - begin, end - begin);
    }
    else {
        if(t < peakVelocities_.begin()->first) {
            return peakVelocities_.begin()->second;
        }

        if(t > peakVelocities_.rbegin()->first) {
            return peakVelocities_.rbegin()->second;
        }
    }

    struct Comparator {
        bool operator()(const std::pair<float, float>& p, float value) {
            return p.first < value;
        }
    };

    auto upper = std::lower_bound(peakVelocities_. begin(), peakVelocities_.end(), t, Comparator());
    auto lower = upper++;

    float a = (t - lower->first) / (upper->first - lower->first);

    return (1.0f - a) * lower->second + a * upper->second;
}

float& VelocityCurve::operator[](float t) {
    return peakVelocities_[t];
}

void VelocityCurve::setPeriodic(bool enabled) {
    periodic_ = enabled;
}

bool VelocityCurve::isPeriodic() const {
    return periodic_;
}

float VelocityCurve::getMinVelocity() const {
    float min = peakVelocities_.begin()->second;
    for(auto iter = ++peakVelocities_.begin(); iter != peakVelocities_.end(); iter++) {
        min = std::min(iter->second, min);
    }
    return min;
}

float VelocityCurve::getMaxVelocity() const {
    float max = peakVelocities_.begin()->second;
    for(auto iter = ++peakVelocities_.begin(); iter != peakVelocities_.end(); iter++) {
        max = std::max(iter->second, max);
    }
    return max;
}

VelocityCurve VelocityCurve::createConstantCurve(float value) {
    VelocityCurve curve;
    curve[0.0f] = value;
    return curve;
}

VelocityCurve VelocityCurve::createLinearCurve(float duration, float maxValue) {
    VelocityCurve curve;
    curve[duration] = maxValue;
    return curve;
}

VelocityCurve VelocityCurve::createSinusoidalCurve(float duration, float maxValue, int steps) {
    tgtAssert(steps > 0, "invalid steps");

    VelocityCurve curve;
    for(int i=0; i<=steps; i++) {
        float ts = i * duration / steps;
        float value = maxValue * (std::sin(-tgt::PIf * 0.5f + i * tgt::PIf / steps) + 1.0f) * 0.5f;
        curve[ts] = value;
    }

    return curve;
}

VelocityCurve VelocityCurve::createHumanHeartBeat() {
    VelocityCurve curve;

    curve[0.040f] =  1200.0f / 20000.0f * 1200.0f;
    curve[0.075f] =  7000.0f / 20000.0f * 1200.0f;
    curve[0.100f] = 17000.0f / 20000.0f * 1200.0f;
    curve[0.125f] = 20000.0f / 20000.0f * 1200.0f;
    curve[0.150f] = 20000.0f / 20000.0f * 1200.0f;
    curve[0.250f] = 10000.0f / 20000.0f * 1200.0f;
    curve[0.280f] =  7500.0f / 20000.0f * 1200.0f;
    curve[0.320f] =  4200.0f / 20000.0f * 1200.0f;
    curve[0.350f] =  1000.0f / 20000.0f * 1200.0f;
    curve[0.380f] =     0.0f / 20000.0f * 1200.0f;
    curve[0.390f] =     0.0f / 20000.0f * 1200.0f;
    curve[0.420f] =   800.0f / 20000.0f * 1200.0f;
    curve[0.455f] =  1200.0f / 20000.0f * 1200.0f;
    curve[0.495f] =  2500.0f / 20000.0f * 1200.0f;
    curve[0.525f] =  2000.0f / 20000.0f * 1200.0f;
    curve[0.560f] =  1200.0f / 20000.0f * 1200.0f;
    curve[0.595f] =   300.0f / 20000.0f * 1200.0f;
    curve[0.660f] =   600.0f / 20000.0f * 1200.0f;
    curve[0.700f] =   300.0f / 20000.0f * 1200.0f;

    return curve;
}

VelocityCurve VelocityCurve::createFromCSV(const std::string& file) {
    VelocityCurve curve;
    curve.peakVelocities_.clear();

    std::ifstream lineStream(file.c_str());
    if (lineStream.fail()) {
        throw VoreenException("CSV file could not be opened");
    }

    std::string line;
    while(std::getline(lineStream, line)) {

        std::vector<float> cellValues;

        std::stringstream cellStream;
        std::string cell;
        while(std::getline(cellStream, cell, ',')) {
            float value;
            std::stringstream helper(cell);
            helper >> value;

            if (!helper.good()) {
                throw VoreenException("Invalid value: " + cell);
            }

            cellValues.push_back(value);
        }

        if(cellValues.size() != 2) {
            throw VoreenException("Invalid tupel size");
        }

        curve[cellValues[0]] = cellValues[1];
    }

    if(curve.peakVelocities_.empty()) {
        throw VoreenException("Empty curve");
    }

    return curve;
}

void VelocityCurve::serialize(Serializer& s) const {
    s.serialize("peakVelocities", peakVelocities_);
    s.serialize("periodic", periodic_);
}
void VelocityCurve::deserialize(Deserializer& s) {
    s.deserialize("peakVelocities", peakVelocities_);
    s.deserialize("periodic", periodic_);
}


FlowIndicator::FlowIndicator()
    : type_(FIT_INVALID)
    , id_(-1)
    , center_(tgt::vec3::zero)
    , normal_(tgt::vec3::zero)
    , radius_(0.0f)
    , flowProfile_(FP_NONE)
{
}

void FlowIndicator::serialize(Serializer& s) const {
    s.serialize("type_", type_); // "type" is a reserved xml keyword.
    s.serialize("id_", id_); // "id" is a reserved xml keyword.
    s.serialize("name", name_);
    s.serialize("center", center_);
    s.serialize("normal", normal_);
    s.serialize("radius", radius_);
    s.serialize("flowProfile", flowProfile_);
    s.serialize("velocityCurve", velocityCurve_);
}

void FlowIndicator::deserialize(Deserializer& s) {
    int type = FIT_CANDIDATE;
    s.deserialize("type_", type);
    type_ = static_cast<FlowIndicatorType>(type);
    s.deserialize("id_", id_);
    s.optionalDeserialize("name", name_, "Indicator " + std::to_string(id_));
    s.deserialize("center", center_);
    s.deserialize("normal", normal_);
    s.deserialize("radius", radius_);
    int profile = FP_POISEUILLE;
    s.deserialize("flowProfile", profile);
    flowProfile_ = static_cast<FlowProfile>(profile);
    s.deserialize("velocityCurve", velocityCurve_);
}

FlowParameterSet::FlowParameterSet()
    : FlowParameterSet("")
{
}

FlowParameterSet::FlowParameterSet(const std::string& name)
    : name_(name)
    , spatialResolution_(0)
    , relaxationTime_(0.0f)
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

float FlowParameterSet::getRelaxationTime() const {
    return relaxationTime_;
}

void FlowParameterSet::setRelaxationTime(float relaxationTime) {
    relaxationTime_ = relaxationTime;
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

float FlowParameterSet::getReynoldsNumber() const {
    return characteristicVelocity_ * characteristicLength_ / (viscosity_ / density_);
}

bool FlowParameterSet::isValid() const {
    float dx = (characteristicLength_ / spatialResolution_);
    float convertVelocity = 3.0f / (relaxationTime_ - 0.5f) * viscosity_ / density_ / dx;
    float uLatticeMax = characteristicVelocity_ / convertVelocity;
    return uLatticeMax < 0.4f; // Intrinsic property of LBM.
}

void FlowParameterSet::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("spatialResolution", spatialResolution_);
    s.serialize("relaxationTime", relaxationTime_);
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
    s.deserialize("relaxationTime", relaxationTime_);
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
    , outputFileFormat_(origin.outputFileFormat_)
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

const std::string& FlowParameterSetEnsemble::getOutputFileFormat() const {
    return outputFileFormat_;
}

void FlowParameterSetEnsemble::setOutputFileFormat(const std::string& format) {
    notifyPendingDataInvalidation();
    outputFileFormat_ = format;
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
    s.deserialize("outputFileFormat", outputFileFormat_);
    s.deserialize("flowFeatures", flowFeatures_);
    s.deserialize("flowIndicators", flowIndicators_);
    s.deserialize("flowParametrizations", flowParameterSets);
}

void FlowParameterSetEnsemble::serializeInternal(Serializer& s, size_t param) const {
    s.serialize("name", name_);
    s.serialize("simulationTime", simulationTime_);
    s.serialize("numTimeSteps", numTimeSteps_);
    s.serialize("outputResolution", outputResolution_);
    s.serialize("outputFileFormat", outputFileFormat_);
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
