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

#include "binarizationfilterproperties.h"

namespace voreen {

BinarizationFilterProperties::BinarizationFilterProperties()
    : threshold_(getId("threshold"), "Binarization Threshold", 0.5f, 0.0f, 1.0f, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
{
    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);

    // Add properties to list.
    addProperties();
}

std::string BinarizationFilterProperties::getVolumeFilterName() const {
    return "Binarization";
}

void BinarizationFilterProperties::adjustPropertiesToInput(const VolumeBase& input) {
    if (!input.hasDerivedData<VolumeMinMax>()) {
        LINFO("Calculating VolumeMinMax. This may take a while...");
    }
    const VolumeMinMax* mm = input.getDerivedData<VolumeMinMax>();

    threshold_.setMinValue(mm->getMin());
    threshold_.setMaxValue(mm->getMax());
}

VolumeFilter* BinarizationFilterProperties::getVolumeFilter(const VolumeBase& volume, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);

    RealWorldMapping rwm;
    if (volume.hasMetaData(VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING)) {
        rwm = volume.getRealWorldMapping();
    }

    return new BinarizationFilter(rwm.realWorldToNormalized(settings.threshold_));
}
void BinarizationFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    threshold_.set(settings.threshold_);
}
void BinarizationFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.threshold_ = threshold_.get();
}
void BinarizationFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void BinarizationFilterProperties::addProperties() {
    properties_.push_back(&threshold_);
}
void BinarizationFilterProperties::serialize(Serializer& s) const {
    std::vector<int> names;
    std::vector<Settings> settings;
    for (const auto& pair : instanceSettings_) {
        names.push_back(pair.first);
        settings.push_back(pair.second);
    }
    s.serializeBinaryBlob(getId("names"), names);
    s.serializeBinaryBlob(getId("settings"), settings);
}
void BinarizationFilterProperties::deserialize(Deserializer& s) {
    std::vector<int> names;
    std::vector<Settings> settings;
    s.deserializeBinaryBlob(getId("names"), names);
    s.deserializeBinaryBlob(getId("settings"), settings);
    tgtAssert(names.size() == settings.size(), "number of keys and values does not match");
    for (size_t i = 0; i < names.size(); i++) {
        instanceSettings_[names[i]] = settings[i];
    }
}
std::vector<int> BinarizationFilterProperties::getStoredInstances() const {
    std::vector<int> output;
    for(auto& kv : instanceSettings_) {
        if(kv.first != DEFAULT_SETTINGS) {
            output.push_back(kv.first);
        }
    }
    return output;
}

}
