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

#include "rescalefilterproperties.h"

namespace voreen {

RescaleFilterProperties::RescaleFilterProperties()
    : rescaleStrategyType_(getId("rescaleStrategyType"), "Rescale Strategy", Processor::INVALID_RESULT)
{
    rescaleStrategyType_.addOption("logarithmic", "log(x)", RESCALE_LOGARITHMIC_T);
    rescaleStrategyType_.addOption("exponential", "exp(x)", RESCALE_EXPONENTIAL_T);

    // Update property state.
    rescaleStrategyType_.invalidate();

    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);

    // Add properties to list.
    addProperties();
}

std::string RescaleFilterProperties::getVolumeFilterName() const {
    return "Rescale Filter";
}

void RescaleFilterProperties::adjustPropertiesToInput(const VolumeBase& input) {
}

VolumeFilter* RescaleFilterProperties::getVolumeFilter(const VolumeBase& volume, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);

    // Currently, only 1D rescale is supported.
    return new RescaleFilter1D(
            settings.rescaleStrategyType_,
            volume.getBaseType()
    );
}
void RescaleFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    rescaleStrategyType_.selectByValue(settings.rescaleStrategyType_);
}
void RescaleFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.rescaleStrategyType_ = rescaleStrategyType_.getValue();
}
void RescaleFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void RescaleFilterProperties::addProperties() {
    properties_.push_back(&rescaleStrategyType_);
}
void RescaleFilterProperties::serialize(Serializer& s) const {
    s.serialize(getId("instanceSettings"), instanceSettings_);
}
void RescaleFilterProperties::deserialize(Deserializer& s) {
    try {
        s.deserialize(getId("instanceSettings"), instanceSettings_);
    }
    catch (SerializationException&) {
        s.removeLastError();
        LERROR("You need to reconfigure " << getVolumeFilterName() << " instances of " << ( properties_[0]->getOwner() ? properties_[0]->getOwner()->getGuiName() : "VolumeFilterList"));
    }
}
std::vector<int> RescaleFilterProperties::getStoredInstances() const {
    std::vector<int> output;
    for(auto& kv : instanceSettings_) {
        if(kv.first != DEFAULT_SETTINGS) {
            output.push_back(kv.first);
        }
    }
    return output;
}

void RescaleFilterProperties::Settings::serialize(Serializer& s) const {
    s.serialize("thresholdValue", thresholdValue_);
    s.serialize("replacementValue", replacementValue_);
    s.serialize("rescaleStrategyType", rescaleStrategyType_);
}
void RescaleFilterProperties::Settings::deserialize(Deserializer& s) {
    s.deserialize("thresholdValue", thresholdValue_);
    s.deserialize("replacementValue", replacementValue_);
    int rescaleStrategyType = 0;
    s.deserialize("rescaleStrategyType", rescaleStrategyType);
    rescaleStrategyType_ = static_cast<RescaleStrategyType>(rescaleStrategyType);
}

}
