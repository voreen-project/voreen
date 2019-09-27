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

#include "binarymedianfilterproperties.h"

namespace voreen {

BinaryMedianFilterProperties::BinaryMedianFilterProperties()
    : useUniformExtent_(getId("useUniformExtent"), "Uniform Extent", true)
    , extentX_(getId("extentx"), "Extent X", 1, 1, 100)
    , extentY_(getId("extenty"), "Extent Y", 1, 1, 100)
    , extentZ_(getId("extentz"), "Extent Z", 1, 1, 100)
    , binarizationThreshold_(getId("binarizationThreshold"), "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , samplingStrategyType_(getId("samplingStrategyType"), "Sampling Strategy", SamplingStrategyType::CLAMP_T)
    , outsideVolumeValue_(getId("outsideVolumeValue"), "Outside Volume Value", 0, 0, 1)
    , forceMedian_(getId("forceMedian"), "Force Median", true)
    , objectVoxelThreshold_(getId("objectVoxelThreshold"), "Object Voxel Threshold", 0, 0, std::numeric_limits<int>::max())
{
    samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
    samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
    samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
    ON_CHANGE_LAMBDA(samplingStrategyType_, [this]() {
        outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
    });

    ON_CHANGE(forceMedian_, BinaryMedianFilterProperties, updateObjectVoxelThreshold);

    ON_CHANGE_LAMBDA(extentX_, [this]() {
        if (useUniformExtent_.get()) {
            extentY_.set(extentX_.get());
            extentZ_.set(extentX_.get());
        }
        updateObjectVoxelThreshold();
    });
    ON_CHANGE_LAMBDA(extentY_, [this]() {
        if (useUniformExtent_.get()) {
            extentX_.set(extentY_.get());
            extentZ_.set(extentY_.get());
        }
        updateObjectVoxelThreshold();
    });
    ON_CHANGE_LAMBDA(extentZ_, [this]() {
        if (useUniformExtent_.get()) {
            extentX_.set(extentZ_.get());
            extentY_.set(extentZ_.get());
        }
        updateObjectVoxelThreshold();
    });

    // Update property state.
    samplingStrategyType_.invalidate();
    updateObjectVoxelThreshold();

    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);

    // Add properties to list.
    addProperties();
}

std::string BinaryMedianFilterProperties::getVolumeFilterName() const {
    return "Binary Median Filter";
}

void BinaryMedianFilterProperties::adjustPropertiesToInput(const VolumeBase& input) {
    if (!input.hasDerivedData<VolumeMinMax>()) {
        LINFO("Calculating VolumeMinMax. This may take a while...");
    }
    const VolumeMinMax* mm = input.getDerivedData<VolumeMinMax>();

    binarizationThreshold_.setMinValue(mm->getMin());
    binarizationThreshold_.setMaxValue(mm->getMax());
    binarizationThreshold_.adaptDecimalsToRange(2);
    outsideVolumeValue_.setMinValue(mm->getMin());
    outsideVolumeValue_.setMaxValue(mm->getMax());
}

VolumeFilter* BinaryMedianFilterProperties::getVolumeFilter(const VolumeBase& volume, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    Settings settings = instanceSettings_.at(instanceId);
    RealWorldMapping rwm;
    if (volume.hasMetaData(VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING)) {
        rwm = volume.getRealWorldMapping();
    }
    return new BinaryMedianFilter(
        tgt::ivec3(settings.extentX_, settings.extentY_, settings.extentZ_),
        rwm.realWorldToNormalized(settings.binarizationThreshold_),
        settings.objectVoxelThreshold_,
        SamplingStrategy<float>(settings.samplingStrategyType_, settings.outsideVolumeValue_)
    );
}
void BinaryMedianFilterProperties::restoreInstance(int instanceId) {
    auto iter = instanceSettings_.find(instanceId);
    if (iter == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    Settings settings = instanceSettings_[instanceId];
    useUniformExtent_.set(settings.useUniformExtent_);
    extentX_.set(settings.extentX_);
    extentY_.set(settings.extentY_);
    extentZ_.set(settings.extentZ_);
    binarizationThreshold_.set(settings.binarizationThreshold_);
    samplingStrategyType_.selectByValue(settings.samplingStrategyType_);
    outsideVolumeValue_.set(settings.outsideVolumeValue_);
    forceMedian_.set(settings.forceMedian_);
    objectVoxelThreshold_.set(settings.objectVoxelThreshold_);
}
void BinaryMedianFilterProperties::storeInstance(int instanceId) {
    Settings& settings = instanceSettings_[instanceId];
    settings.useUniformExtent_ = useUniformExtent_.get();
    settings.extentX_ = extentX_.get();
    settings.extentY_ = extentY_.get();
    settings.extentZ_ = extentZ_.get();
    settings.binarizationThreshold_ = binarizationThreshold_.get();
    settings.samplingStrategyType_ = samplingStrategyType_.getValue();
    settings.outsideVolumeValue_ = outsideVolumeValue_.get();
    settings.forceMedian_ = forceMedian_.get();
    settings.objectVoxelThreshold_ = objectVoxelThreshold_.get();
}
void BinaryMedianFilterProperties::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}
void BinaryMedianFilterProperties::addProperties() {
    properties_.push_back(&useUniformExtent_);
    properties_.push_back(&extentX_);
    properties_.push_back(&extentY_);
    properties_.push_back(&extentZ_);
    properties_.push_back(&binarizationThreshold_);
    properties_.push_back(&samplingStrategyType_);
    properties_.push_back(&outsideVolumeValue_);
    properties_.push_back(&forceMedian_);
    properties_.push_back(&objectVoxelThreshold_);
}
void BinaryMedianFilterProperties::serialize(Serializer& s) const {
    std::vector<int> names;
    std::vector<Settings> settings;
    for (const auto& pair : instanceSettings_) {
        names.push_back(pair.first);
        settings.push_back(pair.second);
    }
    s.serializeBinaryBlob(getId("names"), names);
    s.serializeBinaryBlob(getId("settings"), settings);
}
void BinaryMedianFilterProperties::deserialize(Deserializer& s) {
    std::vector<int> names;
    std::vector<Settings> settings;
    s.deserializeBinaryBlob(getId("names"), names);
    s.deserializeBinaryBlob(getId("settings"), settings);
    tgtAssert(names.size() == settings.size(), "number of keys and values does not match");
    for (size_t i = 0; i < names.size(); i++) {
        instanceSettings_[names[i]] = settings[i];
    }
}

void BinaryMedianFilterProperties::updateObjectVoxelThreshold() {
    bool medianForced = forceMedian_.get();
    objectVoxelThreshold_.setReadOnlyFlag(medianForced);
    objectVoxelThreshold_.setMaxValue((2 * extentX_.get() + 1)*(2 * extentY_.get() + 1)*(2 * extentZ_.get() + 1));
    if (medianForced) {
        objectVoxelThreshold_.set(objectVoxelThreshold_.getMaxValue() / 2);
    }
}

}
