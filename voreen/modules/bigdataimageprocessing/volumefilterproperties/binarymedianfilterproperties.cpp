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

#include "binarymedianfilterproperties.h"

namespace voreen {

BinaryMedianFilterProperties::BinaryMedianFilterProperties()
    : useUniformExtent_(getId("useUniformExtent"), "Uniform Extent", true)
    , extentX_(getId("extentx"), "Extent X", 1, 1, 100)
    , extentY_(getId("extenty"), "Extent Y", 1, 1, 100)
    , extentZ_(getId("extentz"), "Extent Z", 1, 1, 100)
    , binarizationThreshold_(getId("binarizationThreshold"), "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::DYNAMIC, Property::LOD_ADVANCED)
    , samplingStrategyType_(getId("samplingStrategyType"), "Sampling Strategy", Processor::INVALID_RESULT)
    , outsideVolumeValue_(getId("outsideVolumeValue"), "Outside Volume Value", 0, 0, 1, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
    , forceMedian_(getId("forceMedian"), "Force Median", true)
    , objectVoxelThreshold_(getId("objectVoxelThreshold"), "Object Voxel Threshold", 0, 0, std::numeric_limits<int>::max(), Processor::INVALID_RESULT, IntProperty::DYNAMIC)
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
    s.serialize(getId("instanceSettings"), instanceSettings_);
}
void BinaryMedianFilterProperties::deserialize(Deserializer& s) {
    try {
        s.deserialize(getId("instanceSettings"), instanceSettings_);
    }
    catch (SerializationException&) {
        s.removeLastError();
        LERROR("You need to reconfigure " << getVolumeFilterName() << " instances of " << ( properties_[0]->getOwner() ? properties_[0]->getOwner()->getGuiName() : "VolumeFilterList"));
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
std::vector<int> BinaryMedianFilterProperties::getStoredInstances() const {
    std::vector<int> output;
    for(auto& kv : instanceSettings_) {
        if(kv.first != DEFAULT_SETTINGS) {
            output.push_back(kv.first);
        }
    }
    return output;
}

void BinaryMedianFilterProperties::Settings::serialize(Serializer& s) const {
    s.serialize("useUniformExtent_", useUniformExtent_);
    s.serialize("extentX", extentX_);
    s.serialize("extentY", extentY_);
    s.serialize("extentZ", extentZ_);
    s.serialize("binarizationThreshold", binarizationThreshold_);
    s.serialize("samplingStrategyType", samplingStrategyType_);
    s.serialize("outsideVolumeValue", outsideVolumeValue_);
    s.serialize("forceMedian", forceMedian_);
    s.serialize("objectVoxelThreshold", objectVoxelThreshold_);
}
void BinaryMedianFilterProperties::Settings::deserialize(Deserializer& s) {
    s.deserialize("useUniformExtent_", useUniformExtent_);
    s.deserialize("extentX", extentX_);
    s.deserialize("extentY", extentY_);
    s.deserialize("extentZ", extentZ_);
    s.deserialize("binarizationThreshold", binarizationThreshold_);
    int samplingStrategyType = 0;
    s.deserialize("samplingStrategyType", samplingStrategyType);
    samplingStrategyType_ = static_cast<SamplingStrategyType>(samplingStrategyType);
    s.deserialize("outsideVolumeValue", outsideVolumeValue_);
    s.deserialize("forceMedian", forceMedian_);
    s.deserialize("objectVoxelThreshold", objectVoxelThreshold_);
}

}
