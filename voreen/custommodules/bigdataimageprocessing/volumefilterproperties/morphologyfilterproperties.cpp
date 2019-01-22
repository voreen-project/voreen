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

#include "morphologyfilterproperties.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

namespace voreen {

MorphologyFilterProperties::MorphologyFilterProperties()
    : extentX_(getId("extentx"), "Extent X", 1, 1, 10)
    , extentY_(getId("extenty"), "Extent Y", 1, 1, 10)
    , extentZ_(getId("extentz"), "Extent Z", 1, 1, 10)
    , morphologyOperatorType_(getId("morphologyOperatorType"), "Operator Type", MorphologyOperatorType::DILATION_T)
    , morphologyOperatorShape_(getId("morphologyOperatorShape"), "Operator Shape", MorphologyOperatorShape::CUBE_T)
    , samplingStrategyType_(getId("samplingStrategyType"), "Sampling Strategy", SamplingStrategyType::CLAMP_T)
    , outsideVolumeValue_(getId("outsideVolumeValue"), "Outside Volume Value", 0, 0, 1)
{
    morphologyOperatorType_.addOption("dilation", "Dilation", MorphologyOperatorType::DILATION_T);
    morphologyOperatorType_.addOption("erosion", "Erosion", MorphologyOperatorType::EROSION_T);

    morphologyOperatorShape_.addOption("cube", "Cube", MorphologyOperatorShape::CUBE_T);
    morphologyOperatorShape_.addOption("sphere", "Sphere", MorphologyOperatorShape::SPHERE_T);

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

std::string MorphologyFilterProperties::getVolumeFilterName() const {
    return "Morphology Filter";
}

void MorphologyFilterProperties::adjustPropertiesToInput(const VolumeBase& input) {
    if (!input.hasDerivedData<VolumeMinMax>()) {
        LINFO("Calculating VolumeMinMax. This may take a while...");
    }
    const VolumeMinMax* mm = input.getDerivedData<VolumeMinMax>();

    outsideVolumeValue_.setMinValue(mm->getMin());
    outsideVolumeValue_.setMaxValue(mm->getMax());
}

VolumeFilter* MorphologyFilterProperties::getVolumeFilter(const VolumeBase& volume, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);
    return new MorphologyFilter(
        tgt::ivec3(settings.extentX_, settings.extentY_, settings.extentZ_),
        settings.morphologyOperatorType_,
        settings.morphologyOperatorShape_,
        SamplingStrategy<float>(settings.samplingStrategyType_, static_cast<float>(settings.outsideVolumeValue_)),
        volume.getBaseType()
    );
}
void MorphologyFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    extentX_.set(settings.extentX_);
    extentY_.set(settings.extentY_);
    extentZ_.set(settings.extentZ_);
    morphologyOperatorType_.selectByValue(settings.morphologyOperatorType_);
    morphologyOperatorShape_.selectByValue(settings.morphologyOperatorShape_);
    samplingStrategyType_.selectByValue(settings.samplingStrategyType_);
    outsideVolumeValue_.set(settings.outsideVolumeValue_);
}
void MorphologyFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.extentX_ = extentX_.get();
    settings.extentY_ = extentY_.get();
    settings.extentZ_ = extentZ_.get();
    settings.morphologyOperatorType_ = morphologyOperatorType_.getValue();
    settings.morphologyOperatorShape_ = morphologyOperatorShape_.getValue();
    settings.samplingStrategyType_ = samplingStrategyType_.getValue();
    settings.outsideVolumeValue_ = outsideVolumeValue_.get();
}
void MorphologyFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void MorphologyFilterProperties::addProperties() {
    properties_.push_back(&extentX_);
    properties_.push_back(&extentY_);
    properties_.push_back(&extentZ_);
    properties_.push_back(&morphologyOperatorType_);
    properties_.push_back(&morphologyOperatorShape_);
    properties_.push_back(&samplingStrategyType_);
    properties_.push_back(&outsideVolumeValue_);
}
void MorphologyFilterProperties::serialize(Serializer& s) const {
    std::vector<int> names;
    std::vector<Settings> settings;
    for (const auto& pair : instanceSettings_) {
        names.push_back(pair.first);
        settings.push_back(pair.second);
    }
    s.serializeBinaryBlob(getId("names"), names);
    s.serializeBinaryBlob(getId("settings"), settings);
}
void MorphologyFilterProperties::deserialize(Deserializer& s) {
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