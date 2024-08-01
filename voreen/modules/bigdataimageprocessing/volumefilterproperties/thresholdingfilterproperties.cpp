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

#include "thresholdingfilterproperties.h"
#include "../volumefiltering/slicereader.h"

namespace voreen {

ThresholdingFilterSettings::ThresholdingFilterSettings()
    : thresholdValue_(settingsId<ThresholdingFilterSettings>("thresholdValue"), "Threshold Value", 0, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
    , replacementValue_(settingsId<ThresholdingFilterSettings>("replacementValue"), "Replacement Value", 0, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
    , thresholdingStrategyType_(settingsId<ThresholdingFilterSettings>("thresholdingStrategyType"), "Thresholding Strategy", Processor::INVALID_RESULT)
{
    thresholdingStrategyType_.addOption("lower", "Lower", ThresholdingStrategyType::LOWER_T);
    thresholdingStrategyType_.addOption("upper", "Upper", ThresholdingStrategyType::UPPER_T);
}
ThresholdingFilterSettings& ThresholdingFilterSettings::operator=(const ThresholdingFilterSettings& other) {
    copyPropertyValue(other.thresholdValue_, thresholdValue_);
    copyPropertyValue(other.replacementValue_, replacementValue_);
    copyPropertyValue(other.thresholdingStrategyType_, thresholdingStrategyType_);

    return *this;
}

std::string ThresholdingFilterSettings::getVolumeFilterName() {
    return "Thresholding Filter";
}

void ThresholdingFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    const auto& mm = input.estimateMinMax();

    thresholdValue_.setMinValue(mm.x);
    thresholdValue_.setMaxValue(mm.y);
    replacementValue_.setMinValue(mm.x);
    replacementValue_.setMaxValue(mm.y);
}

VolumeFilter* ThresholdingFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    RealWorldMapping rwm = inputmetadata.getRealworldMapping();

    // Currently, only 1D thresholding is supported.
    return new ThresholdingFilter1D(
            rwm.realWorldToNormalized(thresholdValue_.get()),
            rwm.realWorldToNormalized(replacementValue_.get()),
            thresholdingStrategyType_.getValue()
    );
}
void ThresholdingFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&thresholdValue_);
    output.push_back(&replacementValue_);
    output.push_back(&thresholdingStrategyType_);
}
void ThresholdingFilterSettings::serialize(Serializer& s) const {
    s.serialize("thresholdValue", thresholdValue_);
    s.serialize("replacementValue", replacementValue_);
    s.serialize("thresholdingStrategyType", thresholdingStrategyType_);
}
void ThresholdingFilterSettings::deserialize(Deserializer& s) {
    deserializeTemplatePropertyWithValueFallback(s, "thresholdValue", thresholdValue_);
    deserializeTemplatePropertyWithValueFallback(s, "replacementValue", replacementValue_);

    try {
        s.deserialize("thresholdingStrategyType", thresholdingStrategyType_);
    } catch (SerializationException&) {
        s.removeLastError();
        int thresholdingStrategyType = 0;
        s.deserialize("thresholdingStrategyType", thresholdingStrategyType);
        thresholdingStrategyType_.selectByValue(static_cast<ThresholdingStrategyType>(thresholdingStrategyType));
    }
}

}
