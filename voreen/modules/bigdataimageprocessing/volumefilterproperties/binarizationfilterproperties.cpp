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

#include "binarizationfilterproperties.h"
#include "../volumefiltering/slicereader.h"

#include "../volumefiltering/binarizationfilter.h"

namespace voreen {
BinarizationFilterSettings::BinarizationFilterSettings()
    : threshold_(settingsId<BinarizationFilterSettings>("threshold"), "Binarization Threshold", 0.5f, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
{
}

BinarizationFilterSettings& BinarizationFilterSettings::operator=(const BinarizationFilterSettings& other) {
    copyPropertyValue(other.threshold_, threshold_);

    return *this;
}

std::string BinarizationFilterSettings::getVolumeFilterName() {
    return "Binarization";
}

void BinarizationFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    const auto& mm = input.estimateMinMax();

    threshold_.setMinValue(mm.x);
    threshold_.setMaxValue(mm.y);
}

VolumeFilter* BinarizationFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    return new BinarizationFilter(inputmetadata.getRealworldMapping().realWorldToNormalized(threshold_.get()));
}
void BinarizationFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&threshold_);
}

void BinarizationFilterSettings::serialize(Serializer& s) const {
    s.serialize("threshold", threshold_);
}
void BinarizationFilterSettings::deserialize(Deserializer& s) {
    deserializeTemplatePropertyWithValueFallback(s, "threshold", threshold_);
}

}
