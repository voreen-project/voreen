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

#include "thresholdingfilterproperties.h"
#include "../volumefiltering/slicereader.h"

namespace voreen {

ThresholdingFilterProperties::ThresholdingFilterProperties()
    : thresholdValue_(getId("thresholdValue"), "Threshold Value", 0, 0, 1, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
    , replacementValue_(getId("replacementValue"), "Replacement Value", 0, 0, 1, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
    , thresholdingStrategyType_(getId("thresholdingStrategyType"), "Thresholding Strategy", Processor::INVALID_RESULT)
{
    thresholdingStrategyType_.addOption("lower", "Lower", ThresholdingStrategyType::LOWER_T);
    thresholdingStrategyType_.addOption("upper", "Upper", ThresholdingStrategyType::UPPER_T);

    // Update property state.
    thresholdingStrategyType_.invalidate();

    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);

    // Add properties to list.
    addProperties();
}

std::string ThresholdingFilterProperties::getVolumeFilterName() const {
    return "Thresholding Filter";
}

void ThresholdingFilterProperties::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    const auto& mm = input.estimateMinMax();

    thresholdValue_.setMinValue(mm.x);
    thresholdValue_.setMaxValue(mm.y);
    replacementValue_.setMinValue(mm.x);
    replacementValue_.setMaxValue(mm.y);
}

VolumeFilter* ThresholdingFilterProperties::getVolumeFilter(const SliceReaderMetaData& inputmetadata, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);
    RealWorldMapping rwm = inputmetadata.getRealworldMapping();

    // Currently, only 1D thresholding is supported.
    return new ThresholdingFilter1D(
            rwm.realWorldToNormalized(settings.thresholdValue_),
            rwm.realWorldToNormalized(settings.replacementValue_),
            settings.thresholdingStrategyType_
    );
}
void ThresholdingFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    thresholdValue_.set(settings.thresholdValue_);
    replacementValue_.set(settings.replacementValue_);
    thresholdingStrategyType_.selectByValue(settings.thresholdingStrategyType_);
}
void ThresholdingFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.thresholdValue_ = thresholdValue_.get();
    settings.replacementValue_ = replacementValue_.get();
    settings.thresholdingStrategyType_ = thresholdingStrategyType_.getValue();
}
void ThresholdingFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void ThresholdingFilterProperties::addProperties() {
    properties_.push_back(&thresholdValue_);
    properties_.push_back(&replacementValue_);
    properties_.push_back(&thresholdingStrategyType_);
}
void ThresholdingFilterProperties::serialize(Serializer& s) const {
    s.serialize(getId("instanceSettings"), instanceSettings_);
}
void ThresholdingFilterProperties::deserialize(Deserializer& s) {
    try {
        s.deserialize(getId("instanceSettings"), instanceSettings_);
    }
    catch (SerializationException&) {
        s.removeLastError();
        LERROR("You need to reconfigure " << getVolumeFilterName() << " instances of " << ( properties_[0]->getOwner() ? properties_[0]->getOwner()->getGuiName() : "VolumeFilterList"));
    }
}
std::vector<int> ThresholdingFilterProperties::getStoredInstances() const {
    std::vector<int> output;
    for(auto& kv : instanceSettings_) {
        if(kv.first != DEFAULT_SETTINGS) {
            output.push_back(kv.first);
        }
    }
    return output;
}

void ThresholdingFilterProperties::Settings::serialize(Serializer& s) const {
    s.serialize("thresholdValue", thresholdValue_);
    s.serialize("replacementValue", replacementValue_);
    s.serialize("thresholdingStrategyType", thresholdingStrategyType_);
}
void ThresholdingFilterProperties::Settings::deserialize(Deserializer& s) {
    s.deserialize("thresholdValue", thresholdValue_);
    s.deserialize("replacementValue", replacementValue_);
    int thresholdingStrategyType = 0;
    s.deserialize("thresholdingStrategyType", thresholdingStrategyType);
    thresholdingStrategyType_ = static_cast<ThresholdingStrategyType>(thresholdingStrategyType);
}

}
