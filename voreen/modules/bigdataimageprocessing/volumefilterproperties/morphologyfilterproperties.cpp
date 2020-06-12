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

#include "morphologyfilterproperties.h"
#include "../volumefiltering/slicereader.h"

namespace voreen {

MorphologyFilterSettings::MorphologyFilterSettings()
    : extentX_(settingsId<MorphologyFilterSettings>("extentx"), "Extent X", 1, 0, 100)
    , extentY_(settingsId<MorphologyFilterSettings>("extenty"), "Extent Y", 1, 0, 100)
    , extentZ_(settingsId<MorphologyFilterSettings>("extentz"), "Extent Z", 1, 0, 100)
    , morphologyOperatorType_(settingsId<MorphologyFilterSettings>("morphologyOperatorType"), "Operator Type", Processor::INVALID_RESULT)
    , morphologyOperatorShape_(settingsId<MorphologyFilterSettings>("morphologyOperatorShape"), "Operator Shape", Processor::INVALID_RESULT)
    , samplingStrategyType_(settingsId<MorphologyFilterSettings>("samplingStrategyType"), "Sampling Strategy", Processor::INVALID_RESULT)
    , outsideVolumeValue_(settingsId<MorphologyFilterSettings>("outsideVolumeValue"), "Outside Volume Value", 0, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
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
}

MorphologyFilterSettings& MorphologyFilterSettings::operator=(const MorphologyFilterSettings& other) {
    copyPropertyValue(other.extentX_, extentX_);
    copyPropertyValue(other.extentY_, extentY_);
    copyPropertyValue(other.extentZ_, extentZ_);
    copyPropertyValue(other.morphologyOperatorType_, morphologyOperatorType_);
    copyPropertyValue(other.morphologyOperatorShape_, morphologyOperatorShape_);
    copyPropertyValue(other.samplingStrategyType_, samplingStrategyType_);
    copyPropertyValue(other.outsideVolumeValue_, outsideVolumeValue_);

    return *this;
}

std::string MorphologyFilterSettings::getVolumeFilterName() {
    return "Morphology Filter";
}

void MorphologyFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    const auto& mm = input.estimateMinMax();

    outsideVolumeValue_.setMinValue(mm.x);
    outsideVolumeValue_.setMaxValue(mm.y);
}

VolumeFilter* MorphologyFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    return new MorphologyFilter(
        tgt::ivec3(extentX_.get(), extentY_.get(), extentZ_.get()),
        morphologyOperatorType_.getValue(),
        morphologyOperatorShape_.getValue(),
        SamplingStrategy<float>(samplingStrategyType_.getValue(), static_cast<float>(outsideVolumeValue_.get()))
    );
}

void MorphologyFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&extentX_);
    output.push_back(&extentY_);
    output.push_back(&extentZ_);
    output.push_back(&morphologyOperatorType_);
    output.push_back(&morphologyOperatorShape_);
    output.push_back(&samplingStrategyType_);
    output.push_back(&outsideVolumeValue_);
}

void MorphologyFilterSettings::serialize(Serializer& s) const {
    s.serialize("extentX", extentX_);
    s.serialize("extentY", extentY_);
    s.serialize("extentZ", extentZ_);
    s.serialize("morphologyOperatorType", morphologyOperatorType_);
    s.serialize("morphologyOperatorShape", morphologyOperatorShape_);
    s.serialize("samplingStrategyType", samplingStrategyType_);
    s.serialize("outsideVolumeValue", outsideVolumeValue_);
}
void MorphologyFilterSettings::deserialize(Deserializer& s) {
    deserializeTemplatePropertyWithValueFallback(s, "extentX", extentX_);
    deserializeTemplatePropertyWithValueFallback(s, "extentY", extentY_);
    deserializeTemplatePropertyWithValueFallback(s, "extentZ", extentZ_);

    try {
        s.deserialize("morphologyOperatorType", morphologyOperatorType_);
    } catch (SerializationException&) {
        s.removeLastError();
        int morphologyOperatorType = 0;
        s.deserialize("morphologyOperatorType", morphologyOperatorType);
        morphologyOperatorType_.selectByValue(static_cast<MorphologyOperatorType>(morphologyOperatorType));
    }

    try {
        s.deserialize("morphologyOperatorShape", morphologyOperatorShape_);
    } catch (SerializationException&) {
        s.removeLastError();
        int morphologyOperatorShape = 0;
        s.deserialize("morphologyOperatorShape", morphologyOperatorShape);
        morphologyOperatorShape_.selectByValue(static_cast<MorphologyOperatorShape>(morphologyOperatorShape));
    }

    try {
        s.deserialize("samplingStrategyType", samplingStrategyType_);
    } catch (SerializationException&) {
        s.removeLastError();
        int samplingStrategyType = 0;
        s.deserialize("samplingStrategyType", samplingStrategyType);
        samplingStrategyType_.selectByValue(static_cast<SamplingStrategyType>(samplingStrategyType));
    }

    deserializeTemplatePropertyWithValueFallback(s, "outsideVolumeValue", outsideVolumeValue_);
}

}
