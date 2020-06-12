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

#include "vorticityfilterproperties.h"
#include "../volumefiltering/slicereader.h"

namespace voreen {

VorticityFilterSettings::VorticityFilterSettings()
    : gradientType_(settingsId<VorticityFilterSettings>("gradientType"), "Gradient Type", Processor::INVALID_RESULT)
    , samplingStrategyType_(settingsId<VorticityFilterSettings>("samplingStrategyType"), "Sampling Strategy", Processor::INVALID_RESULT)
    , outsideVolumeValue_(settingsId<VorticityFilterSettings>("outsideVolumeValue"), "Outside Volume Value", 0, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
{
    gradientType_.addOption("centralDifferences", "Central Differences", GradientType::VOG_CENTRAL_DIFFERENCE);
    gradientType_.addOption("linearRegression", "Linear Regression", GradientType::VOG_LINEAR_REGRESSION);
    gradientType_.addOption("sobel", "Sobel", GradientType::VOG_SOBEL);

    samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
    samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
    samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
    ON_CHANGE_LAMBDA(samplingStrategyType_, [this]() {
        outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
    });

    // Update property state.
    samplingStrategyType_.invalidate();
}
VorticityFilterSettings& VorticityFilterSettings::operator=(const VorticityFilterSettings& other) {
    copyPropertyValue(other.gradientType_, gradientType_);
    copyPropertyValue(other.samplingStrategyType_, samplingStrategyType_);
    copyPropertyValue(other.outsideVolumeValue_, outsideVolumeValue_);

    return *this;
}

std::string VorticityFilterSettings::getVolumeFilterName() {
    return "Vorticity";
}

void VorticityFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    const auto& mm = input.estimateMinMax();

    outsideVolumeValue_.setMinValue(mm.x);
    outsideVolumeValue_.setMaxValue(mm.y);
}

VolumeFilter* VorticityFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    return new VorticityFilter(
            gradientType_.getValue(),
            inputmetadata.getSpacing(),
            SamplingStrategy<tgt::vec3>(samplingStrategyType_.getValue(), tgt::vec3(outsideVolumeValue_.get()))
    );
}
void VorticityFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&gradientType_);
    output.push_back(&samplingStrategyType_);
    output.push_back(&outsideVolumeValue_);
}

void VorticityFilterSettings::serialize(Serializer& s) const {
    s.serialize("gradientType", gradientType_);
    s.serialize("samplingStrategyType", samplingStrategyType_);
    s.serialize("outsideVolumeValue", outsideVolumeValue_);
    }
void VorticityFilterSettings::deserialize(Deserializer& s) {
    try {
        s.deserialize("gradientType", gradientType_);
    } catch (SerializationException&) {
        s.removeLastError();
        int gradientType = 0;
        s.deserialize("gradientType", gradientType);
        gradientType_.selectByValue(static_cast<GradientType>(gradientType));
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
