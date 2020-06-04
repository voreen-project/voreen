/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "vorticityfilterproperties.h"
#include "../volumefiltering/slicereader.h"

namespace voreen {

VorticityFilterProperties::VorticityFilterProperties()
    : gradientType_(getId("gradientType"), "Gradient Type", Processor::INVALID_RESULT)
    , samplingStrategyType_(getId("samplingStrategyType"), "Sampling Strategy", Processor::INVALID_RESULT)
    , outsideVolumeValue_(getId("outsideVolumeValue"), "Outside Volume Value", 0, 0, 1, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
{
    gradientType_.addOption("centralDifferences", "Central Differences", GradientType::VOG_CENTRAL_DIFFERENCE);
    gradientType_.addOption("linearRegression", "Linear Regression", GradientType::VOG_LINEAR_REGRESSION);
    gradientType_.addOption("sobel", "Sobel", GradientType::VOG_SOBEL);

    samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
    samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
    samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
    ON_CHANGE_LAMBDA(samplingStrategyType_, [this]() {
        outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
    });

    // Update property state.
    samplingStrategyType_.invalidate();

    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);

    // Add properties to list.
    addProperties();
}

std::string VorticityFilterProperties::getVolumeFilterName() const {
    return "Vorticity";
}

void VorticityFilterProperties::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    const auto& mm = input.estimateMinMax();

    outsideVolumeValue_.setMinValue(mm.x);
    outsideVolumeValue_.setMaxValue(mm.y);
}

VolumeFilter* VorticityFilterProperties::getVolumeFilter(const SliceReaderMetaData& inputmetadata, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);
    return new VorticityFilter(
            settings.gradientType_,
            inputmetadata.getSpacing(),
            SamplingStrategy<tgt::vec3>(settings.samplingStrategyType_, tgt::vec3(settings.outsideVolumeValue_))
    );
}
void VorticityFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    gradientType_.selectByValue(settings.gradientType_);
    samplingStrategyType_.selectByValue(settings.samplingStrategyType_);
    outsideVolumeValue_.set(settings.outsideVolumeValue_);
}
void VorticityFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.gradientType_ = gradientType_.getValue();
    settings.samplingStrategyType_ = samplingStrategyType_.getValue();
    settings.outsideVolumeValue_ = outsideVolumeValue_.get();
}
void VorticityFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void VorticityFilterProperties::addProperties() {
    properties_.push_back(&gradientType_);
    properties_.push_back(&samplingStrategyType_);
    properties_.push_back(&outsideVolumeValue_);
}
void VorticityFilterProperties::serialize(Serializer& s) const {
    s.serialize(getId("instanceSettings"), instanceSettings_);
}
void VorticityFilterProperties::deserialize(Deserializer& s) {
    s.deserialize(getId("instanceSettings"), instanceSettings_);
}
std::vector<int> VorticityFilterProperties::getStoredInstances() const {
    std::vector<int> output;
    for(auto& kv : instanceSettings_) {
        if(kv.first != DEFAULT_SETTINGS) {
            output.push_back(kv.first);
        }
    }
    return output;
}

void VorticityFilterProperties::Settings::serialize(Serializer& s) const {
    s.serialize("gradientType", gradientType_);
    s.serialize("samplingStrategyType", samplingStrategyType_);
    s.serialize("outsideVolumeValue", outsideVolumeValue_);
    }
    void VorticityFilterProperties::Settings::deserialize(Deserializer& s) {
        int gradientType = 0;
        s.deserialize("gradientType", gradientType);
        gradientType_ = static_cast<GradientType>(gradientType);
        int samplingStrategyType = 0;
        s.deserialize("samplingStrategyType", samplingStrategyType);
        samplingStrategyType_ = static_cast<SamplingStrategyType>(samplingStrategyType);
        s.deserialize("outsideVolumeValue", outsideVolumeValue_);
    }

}
