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

#include "valuemapfilterproperties.h"
#include "tgt/memory.h"

namespace voreen {

ValueMapFilterSettings::ValueMapFilterSettings()
    : valueMap_("valueMap", "Value Map (alpha value)")
    , fakeValueVol_(nullptr)
    , minmax_(tgt::vec2(0.0f, 1.0f))
{
}
ValueMapFilterSettings& ValueMapFilterSettings::operator=(const ValueMapFilterSettings& other) {
    copyPropertyValue(other.valueMap_, valueMap_);
    updateValueMapRange(other.minmax_);
    return *this;
}

std::string ValueMapFilterSettings::getVolumeFilterName() {
    return "Value Map";
}
void ValueMapFilterSettings::updateValueMapRange(tgt::vec2 minmax) {
    minmax_ = minmax;

    auto data = tgt::make_unique<VolumeAtomic<float>>(tgt::svec3(2,1,1));
    data->voxel()[0] = minmax[0];
    data->voxel()[1] = minmax[1];
    fakeValueVol_.reset(new Volume(data.release(), tgt::svec3(1), tgt::svec3(0)));
    fakeValueVol_->setRealWorldMapping(RealWorldMapping(tgt::vec2(0.0f, 1.0f), ""));
    valueMap_.setVolume(fakeValueVol_.get());
}

void ValueMapFilterSettings::adjustPropertiesToInput(const SliceReaderMetaData& input) {
    auto minmax = input.estimateMinMax();
    updateValueMapRange(minmax);
}

static std::vector<uint8_t> extractLUT(const TransFunc1DKeys& tf) {
    tgt::Texture* texture = tf.getTexture();

    texture->downloadTexture();
    const int mapDim = texture->getDimensions().x;

    std::vector<uint8_t> valueMap;
    valueMap.reserve(mapDim);
    for(size_t i=0; i<mapDim; ++i) {
        valueMap.push_back(texture->texel<tgt::Vector4<uint8_t>>(i).a);
    }
    return valueMap;
}

VolumeFilter* ValueMapFilterSettings::getVolumeFilter(const SliceReaderMetaData& inputmetadata) const {
    // Currently, only 1D remap is supported.
    auto propVal = valueMap_.get();

    std::vector<uint8_t> valueMap;
    RealWorldMapping tfRwm;

    if(propVal) {
        valueMap = extractLUT(*propVal);
        tfRwm = RealWorldMapping(propVal->getDomain(), "");
    } else {
        valueMap = extractLUT(TransFunc1DKeys());
        tfRwm = RealWorldMapping(tgt::vec2(0.0, 1.0), "");
    }

    auto inputToLut = RealWorldMapping::combine(inputmetadata.getRealWorldMapping(), tfRwm.getInverseMapping());
    return new ValueMapFilter1D(std::move(valueMap), inputToLut);
}
void ValueMapFilterSettings::addProperties(std::vector<Property*>& output) {
    output.push_back(&valueMap_);
}

void ValueMapFilterSettings::serialize(Serializer& s) const {
    s.serialize("valueMap", valueMap_);
}
void ValueMapFilterSettings::deserialize(Deserializer& s) {
    s.deserialize("valueMap", valueMap_);

}

}
