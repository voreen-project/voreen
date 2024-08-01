/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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
#include "../volumefiltering/slicereader.h"

#include "../volumefiltering/binarymedianfilter.h"

namespace voreen {

BinaryMedianFilterSettings::BinaryMedianFilterSettings()
    : useUniformExtent_(settingsId<BinaryMedianFilterSettings>("useUniformExtent"), "Uniform Extent", true)
    , extentX_(settingsId<BinaryMedianFilterSettings>("extentx"), "Extent X", 1, 0, 100)
    , extentY_(settingsId<BinaryMedianFilterSettings>("extenty"), "Extent Y", 1, 0, 100)
    , extentZ_(settingsId<BinaryMedianFilterSettings>("extentz"), "Extent Z", 1, 0, 100)
    , binarizationThreshold_(settingsId<BinaryMedianFilterSettings>("binarizationThreshold"), "Threshold", 0.5f, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC, Property::LOD_ADVANCED)
    , samplingStrategyType_(settingsId<BinaryMedianFilterSettings>("samplingStrategyType"), "Sampling Strategy", Processor::INVALID_RESULT)
    , outsideVolumeValue_(settingsId<BinaryMedianFilterSettings>("outsideVolumeValue"), "Outside Volume Value", 0, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
    , forceMedian_(settingsId<BinaryMedianFilterSettings>("forceMedian"), "Force Median", true)
    , objectVoxelThreshold_(settingsId<BinaryMedianFilterSettings>("objectVoxelThreshold"), "Object Voxel Threshold", 0, INT_MIN, INT_MAX, Processor::INVALID_RESULT, IntProperty::DYNAMIC)
{
    samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
    samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
    samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
    ON_CHANGE_LAMBDA(samplingStrategyType_, [this]() {
        outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
    });

    // Update property state.
    samplingStrategyType_.invalidate();
    updateObjectVoxelThreshold();


    ON_CHANGE(forceMedian_, BinaryMedianFilterSettings, updateObjectVoxelThreshold);

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
}

BinaryMedianFilterSettings& BinaryMedianFilterSettings::operator=(const BinaryMedianFilterSettings& other) {
    copyPropertyValue(other.useUniformExtent_, useUniformExtent_);
    copyPropertyValue(other.extentX_, extentX_);
    copyPropertyValue(other.extentY_, extentY_);
    copyPropertyValue(other.extentZ_, extentZ_);
    copyPropertyValue(other.binarizationThreshold_, binarizationThreshold_);
    copyPropertyValue(other.samplingStrategyType_, samplingStrategyType_);
    copyPropertyValue(other.outsideVolumeValue_, outsideVolumeValue_);
    copyPropertyValue(other.forceMedian_, forceMedian_);
    copyPropertyValue(other.objectVoxelThreshold_, objectVoxelThreshold_);

    return *this;
}

std::string BinaryMedianFilterSettings::getVolumeFilterName() {
    return "Binary Median Filter";
}

void BinaryMedianFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    auto mm = input.estimateMinMax();

    binarizationThreshold_.setMinValue(mm.x);
    binarizationThreshold_.setMaxValue(mm.y);
    binarizationThreshold_.adaptDecimalsToRange(2);
    outsideVolumeValue_.setMinValue(mm.x);
    outsideVolumeValue_.setMaxValue(mm.y);
}

VolumeFilter* BinaryMedianFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    return new BinaryMedianFilter(
        tgt::ivec3(extentX_.get(), extentY_.get(), extentZ_.get()),
        inputmetadata.getRealWorldMapping().realWorldToNormalized(binarizationThreshold_.get()),
        objectVoxelThreshold_.get(),
        SamplingStrategy<float>(samplingStrategyType_.getValue(), outsideVolumeValue_.get())
    );
}
void BinaryMedianFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&useUniformExtent_);
    output.push_back(&extentX_);
    output.push_back(&extentY_);
    output.push_back(&extentZ_);
    output.push_back(&binarizationThreshold_);
    output.push_back(&samplingStrategyType_);
    output.push_back(&outsideVolumeValue_);
    output.push_back(&forceMedian_);
    output.push_back(&objectVoxelThreshold_);
}
void BinaryMedianFilterSettings::serialize(Serializer& s) const {
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
void BinaryMedianFilterSettings::deserialize(Deserializer& s) {
    deserializeTemplatePropertyWithValueFallback(s, "useUniformExtent_", useUniformExtent_);
    deserializeTemplatePropertyWithValueFallback(s, "extentX", extentX_);
    deserializeTemplatePropertyWithValueFallback(s, "extentY", extentY_);
    deserializeTemplatePropertyWithValueFallback(s, "extentZ", extentZ_);
    deserializeTemplatePropertyWithValueFallback(s, "binarizationThreshold", binarizationThreshold_);
    deserializeTemplatePropertyWithValueFallback(s, "outsideVolumeValue", outsideVolumeValue_);
    deserializeTemplatePropertyWithValueFallback(s, "forceMedian", forceMedian_);
    deserializeTemplatePropertyWithValueFallback(s, "objectVoxelThreshold", objectVoxelThreshold_);

    try {
        s.deserialize("samplingStrategyType", samplingStrategyType_);
    } catch (SerializationException&) {
        s.removeLastError();
        int samplingStrategyType = 0;
        s.deserialize("samplingStrategyType", samplingStrategyType);
        samplingStrategyType_.selectByValue(static_cast<SamplingStrategyType>(samplingStrategyType));
    }
}

void BinaryMedianFilterSettings::updateObjectVoxelThreshold() {
    bool medianForced = forceMedian_.get();
    objectVoxelThreshold_.setReadOnlyFlag(medianForced);
    objectVoxelThreshold_.setMaxValue((2 * extentX_.get() + 1)*(2 * extentY_.get() + 1)*(2 * extentZ_.get() + 1));
    if (medianForced) {
        objectVoxelThreshold_.set(objectVoxelThreshold_.getMaxValue() / 2);
    }
}

}
