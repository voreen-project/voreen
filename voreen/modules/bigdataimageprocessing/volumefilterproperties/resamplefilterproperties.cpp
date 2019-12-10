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

#include "resamplefilterproperties.h"
#include "../volumefiltering/resamplefilter.h"

namespace voreen {

#define RESAMPLE_DIM_DEFAULT_VALUE tgt::ivec3(100000)

ResampleFilterProperties::ResampleFilterProperties()
    : dimensions_(getId("dimensions"), "Target Dimensions", RESAMPLE_DIM_DEFAULT_VALUE, tgt::ivec3(2), tgt::ivec3(100000))
{
    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);

    // Add properties to list.
    addProperties();
}

std::string ResampleFilterProperties::getVolumeFilterName() const {
    return "Resample Filter";
}

void ResampleFilterProperties::adjustPropertiesToInput(const VolumeBase& input) {
    if(dimensions_.get() == RESAMPLE_DIM_DEFAULT_VALUE) {
        dimensions_.set(input.getDimensions());
    }
}

VolumeFilter* ResampleFilterProperties::getVolumeFilter(const VolumeBase& volume, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);
    return new ResampleFilter(
        settings.dimensions_,
        volume.getBaseType(), volume.getNumChannels()
    );
}
void ResampleFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    dimensions_.set(settings.dimensions_);
}
void ResampleFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.dimensions_ = dimensions_.get();
}
void ResampleFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void ResampleFilterProperties::addProperties() {
    properties_.push_back(&dimensions_);
}
void ResampleFilterProperties::serialize(Serializer& s) const {
    std::vector<int> names;
    std::vector<Settings> settings;
    for (const auto& pair : instanceSettings_) {
        names.push_back(pair.first);
        settings.push_back(pair.second);
    }
    s.serializeBinaryBlob(getId("names"), names);
    s.serializeBinaryBlob(getId("settings"), settings);
}
void ResampleFilterProperties::deserialize(Deserializer& s) {
    std::vector<int> names;
    std::vector<Settings> settings;
    s.deserializeBinaryBlob(getId("names"), names);
    s.deserializeBinaryBlob(getId("settings"), settings);
    tgtAssert(names.size() == settings.size(), "number of keys and values does not match");
    for (size_t i = 0; i < names.size(); i++) {
        instanceSettings_[names[i]] = settings[i];
    }
}

}
