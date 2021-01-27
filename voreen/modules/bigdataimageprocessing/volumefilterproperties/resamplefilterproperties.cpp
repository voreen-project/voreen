/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "resamplefilterproperties.h"
#include "../volumefiltering/resamplefilter.h"
#include "../volumefiltering/slicereader.h"

namespace voreen {

#define RESAMPLE_DIM_DEFAULT_VALUE tgt::ivec3(100000)

ResampleFilterSettings::ResampleFilterSettings()
    : dimensions_(settingsId<ResampleFilterSettings>("dimensions"), "Target Dimensions", RESAMPLE_DIM_DEFAULT_VALUE, tgt::ivec3(2), tgt::ivec3(100000))
{
}
ResampleFilterSettings& ResampleFilterSettings::operator=(const ResampleFilterSettings& other) {
    copyPropertyValue(other.dimensions_, dimensions_);

    return *this;
}

std::string ResampleFilterSettings::getVolumeFilterName() {
    return "Resample Filter";
}

void ResampleFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    if(dimensions_.get() == RESAMPLE_DIM_DEFAULT_VALUE) {
        dimensions_.set(input.getDimensions());
    }
}

VolumeFilter* ResampleFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    return new ResampleFilter(
        dimensions_.get(),
        inputmetadata.getNumChannels()
    );
}
void ResampleFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&dimensions_);
}
void ResampleFilterSettings::serialize(Serializer& s) const {
    s.serialize("dimensions", dimensions_);
}
void ResampleFilterSettings::deserialize(Deserializer& s) {
    deserializeTemplatePropertyWithValueFallback(s, "dimensions", dimensions_);
}
}
