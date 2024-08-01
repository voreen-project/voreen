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

#include "gaussianfilterproperties.h"
#include "../volumefiltering/gaussianfilter.h"

namespace voreen {
GaussianFilterSettings::GaussianFilterSettings()
    : extentX_(settingsId<GaussianFilterSettings>("extentx"), "Extent X", 1, 0, 100)
    , extentY_(settingsId<GaussianFilterSettings>("extenty"), "Extent Y", 1, 0, 100)
    , extentZ_(settingsId<GaussianFilterSettings>("extentz"), "Extent Z", 1, 0, 100)
    , samplingStrategyType_(settingsId<GaussianFilterSettings>("samplingStrategyType"), "Sampling Strategy", Processor::INVALID_RESULT)
    , outsideVolumeValue_(settingsId<GaussianFilterSettings>("outsideVolumeValue"), "Outside Volume Value", 0, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
{
    samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
    samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
    samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
    ON_CHANGE_LAMBDA(samplingStrategyType_, [this]() {
        outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
    });

    // Update property state.
    samplingStrategyType_.invalidate();
}
GaussianFilterSettings& GaussianFilterSettings::operator=(const GaussianFilterSettings& other) {
    copyPropertyValue(other.extentX_, extentX_);
    copyPropertyValue(other.extentY_, extentY_);
    copyPropertyValue(other.extentZ_, extentZ_);
    copyPropertyValue(other.samplingStrategyType_, samplingStrategyType_);
    copyPropertyValue(other.outsideVolumeValue_, outsideVolumeValue_);

    return *this;
}

std::string GaussianFilterSettings::getVolumeFilterName() {
    return "Gaussian Filter";
}

void GaussianFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    const auto& mm = input.estimateMinMax();

    outsideVolumeValue_.setMinValue(mm.x);
    outsideVolumeValue_.setMaxValue(mm.y);
}
VolumeFilter* GaussianFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    return new GaussianFilter(
        tgt::ivec3(extentX_.get(), extentY_.get(), extentZ_.get()),
        SamplingStrategy<float>(samplingStrategyType_.getValue(), static_cast<float>(outsideVolumeValue_.get())),
        inputmetadata.getNumChannels()
    );
}
void GaussianFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&extentX_);
    output.push_back(&extentY_);
    output.push_back(&extentZ_);
    output.push_back(&samplingStrategyType_);
    output.push_back(&outsideVolumeValue_);
}

void GaussianFilterSettings::serialize(Serializer& s) const {
    s.serialize("extentX", extentX_);
    s.serialize("extentY", extentY_);
    s.serialize("extentZ", extentZ_);
    s.serialize("samplingStrategyType", samplingStrategyType_);
    s.serialize("outsideVolumeValue", outsideVolumeValue_);
}
void GaussianFilterSettings::deserialize(Deserializer& s) {
    deserializeTemplatePropertyWithValueFallback(s, "extentX", extentX_);
    deserializeTemplatePropertyWithValueFallback(s, "extentY", extentY_);
    deserializeTemplatePropertyWithValueFallback(s, "extentZ", extentZ_);
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
