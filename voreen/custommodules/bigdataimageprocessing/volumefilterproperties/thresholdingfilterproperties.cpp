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

#include "thresholdingfilterproperties.h"

namespace voreen {

ThresholdingFilterProperties::ThresholdingFilterProperties()
    : thresholdValue_(getId("thresholdValue"), "Threshold Value", 0, 0, 1)
    , binarize_(getId("binarize"), "Binarize?", false)
    , replacementValue_(getId("replacementValue"), "Replacement Value", 0, 0, 1)
    , thresholdingStrategyType_(getId("thresholdingStrategyType"), "Thresholding Strategy", ThresholdingStrategyType::LOWER_T)
{
    thresholdingStrategyType_.addOption("lower", "Lower", ThresholdingStrategyType::LOWER_T);
    thresholdingStrategyType_.addOption("upper", "Upper", ThresholdingStrategyType::UPPER_T);

    ON_CHANGE_LAMBDA(binarize_, [this]() {
        replacementValue_.setVisibleFlag(!binarize_.get());
    });

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

void ThresholdingFilterProperties::adjustPropertiesToInput(const VolumeBase& input) {
    if (!input.hasDerivedData<VolumeMinMax>()) {
        LINFO("Calculating VolumeMinMax. This may take a while...");
    }
    const VolumeMinMax* mm = input.getDerivedData<VolumeMinMax>();

    thresholdValue_.setMinValue(mm->getMin());
    thresholdValue_.setMaxValue(mm->getMax());
    replacementValue_.setMinValue(mm->getMin());
    replacementValue_.setMaxValue(mm->getMax());
}

VolumeFilter* ThresholdingFilterProperties::getVolumeFilter(const VolumeBase& volume, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);
    RealWorldMapping rwm;
    if (volume.hasMetaData("RealWorldMapping")) {
        rwm = volume.getRealWorldMapping();
    }
    if (binarize_.get()) {
        return new ThresholdingFilter(
            rwm.realWorldToNormalized(settings.thresholdValue_),
            settings.thresholdingStrategyType_
        );
    }
    else {
        return new ThresholdingFilter(
            rwm.realWorldToNormalized(settings.thresholdValue_),
            rwm.realWorldToNormalized(settings.replacementValue_),
            settings.thresholdingStrategyType_,
            volume.getBaseType()
        );
    }
}
void ThresholdingFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    thresholdValue_.set(settings.thresholdValue_);
    binarize_.set(settings.binarize_);
    replacementValue_.set(settings.replacementValue_);
    thresholdingStrategyType_.selectByValue(settings.thresholdingStrategyType_);
}
void ThresholdingFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.thresholdValue_ = thresholdValue_.get();
    settings.binarize_ = binarize_.get();
    settings.replacementValue_ = replacementValue_.get();
    settings.thresholdingStrategyType_ = thresholdingStrategyType_.getValue();
}
void ThresholdingFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void ThresholdingFilterProperties::addProperties() {
    properties_.push_back(&thresholdValue_);
    properties_.push_back(&binarize_);
    properties_.push_back(&replacementValue_);
    properties_.push_back(&thresholdingStrategyType_);
}
void ThresholdingFilterProperties::serialize(Serializer& s) const {
    std::vector<int> names;
    std::vector<Settings> settings;
    for (const auto& pair : instanceSettings_) {
        names.push_back(pair.first);
        settings.push_back(pair.second);
    }
    s.serializeBinaryBlob(getId("names"), names);
    s.serializeBinaryBlob(getId("settings"), settings);
}
void ThresholdingFilterProperties::deserialize(Deserializer& s) {
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
