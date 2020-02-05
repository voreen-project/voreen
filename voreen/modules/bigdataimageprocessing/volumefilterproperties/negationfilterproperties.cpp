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

#include "negationfilterproperties.h"

namespace voreen {

NegationFilterProperties::NegationFilterProperties()
    : negateX_(getId("negateX"), "Negate X")
    , negateY_(getId("negateY"), "Negate Y")
    , negateZ_(getId("negateZ"), "Negate Z")
    , negateW_(getId("negateW"), "Negate W")
{
    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);

    // Add properties to list.
    addProperties();
}

std::string NegationFilterProperties::getVolumeFilterName() const {
    return "Negation";
}

void NegationFilterProperties::adjustPropertiesToInput(const VolumeBase& volume) {
    size_t numChannels = volume.getNumChannels();
    //negateX_.setReadOnlyFlag(numChannels < 1); // Always true.
    negateY_.setReadOnlyFlag(numChannels < 2);
    negateZ_.setReadOnlyFlag(numChannels < 3);
    negateW_.setReadOnlyFlag(numChannels < 4);
}

VolumeFilter* NegationFilterProperties::getVolumeFilter(const VolumeBase& volume, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);

    tgt::bvec4 negate(settings.negateX_, settings.negateY_, settings.negateZ_, settings.negateW_);

    switch(volume.getNumChannels()) {
    case 1:
        return new NegationFilter1D(negate, volume.getBaseType());
    case 2:
        return new NegationFilter2D(negate, volume.getBaseType());
    case 3:
        return new NegationFilter3D(negate, volume.getBaseType());
    case 4:
        return new NegationFilter4D(negate, volume.getBaseType());
    default:
        return nullptr;
    }
}
void NegationFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    negateX_.set(settings.negateX_);
    negateY_.set(settings.negateY_);
    negateZ_.set(settings.negateZ_);
    negateW_.set(settings.negateW_);
}
void NegationFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.negateX_ = negateX_.get();
    settings.negateY_ = negateY_.get();
    settings.negateZ_ = negateZ_.get();
    settings.negateW_ = negateW_.get();
}
void NegationFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void NegationFilterProperties::addProperties() {
    properties_.push_back(&negateX_);
    properties_.push_back(&negateY_);
    properties_.push_back(&negateZ_);
    properties_.push_back(&negateW_);
}
void NegationFilterProperties::serialize(Serializer& s) const {
    s.serialize(getId("instanceSettings"), instanceSettings_);
}
void NegationFilterProperties::deserialize(Deserializer& s) {
    s.deserialize(getId("instanceSettings"), instanceSettings_);
}
std::vector<int> NegationFilterProperties::getStoredInstances() const {
    std::vector<int> output;
    for(auto& kv : instanceSettings_) {
        if(kv.first != DEFAULT_SETTINGS) {
            output.push_back(kv.first);
        }
    }
    return output;
}

void NegationFilterProperties::Settings::serialize(Serializer& s) const {
    s.serialize("negateX", negateX_);
    s.serialize("negateY", negateY_);
    s.serialize("negateZ", negateZ_);
    s.serialize("negateW", negateW_);
}
void NegationFilterProperties::Settings::deserialize(Deserializer& s) {
    s.deserialize("negateX", negateX_);
    s.deserialize("negateY", negateY_);
    s.deserialize("negateZ", negateZ_);
    s.deserialize("negateW", negateW_);
}

}
