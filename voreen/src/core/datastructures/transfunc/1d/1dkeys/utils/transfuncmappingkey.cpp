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

#include "voreen/core/datastructures/transfunc/1d/1dkeys/utils/transfuncmappingkey.h"

namespace voreen {

TransFuncMappingKey::TransFuncMappingKey(float i, const tgt::col4& color)
    : intensity_(i)
    , colorL_(color)
    , colorR_(color)
    , split_(false)
{}

TransFuncMappingKey::TransFuncMappingKey()
    : intensity_(0.f)
    , colorL_(tgt::vec4(0.f))
    , colorR_(tgt::vec4(0.f))
    , split_(false)
{}

TransFuncMappingKey::~TransFuncMappingKey() {
}

bool TransFuncMappingKey::operator==(const TransFuncMappingKey& key) {
    return (intensity_ == key.intensity_) && (split_ == key.split_) &&
           (colorR_    == key.colorR_) && (colorL_ == key.colorL_);
}

bool TransFuncMappingKey::operator!=(const TransFuncMappingKey& key) {
    return !(*this == key);
}

void TransFuncMappingKey::setColorL(const tgt::col4& color) {
    colorL_ = color;
    if (!split_)
        colorR_ = color;
}

void TransFuncMappingKey::setColorL(const tgt::ivec4& color) {
    colorL_ = tgt::col4(color);
    if (!split_)
        colorR_ = tgt::col4(color);
}

const tgt::col4& TransFuncMappingKey::getColorL() const {
    return colorL_;
}

tgt::col4& TransFuncMappingKey::getColorL() {
    return colorL_;
}

void TransFuncMappingKey::setColorR(const tgt::col4& color) {
    colorR_ = color;
    if (!split_)
        colorL_ = color;
}

void TransFuncMappingKey::setColorR(const tgt::ivec4& color) {
    colorR_ = tgt::col4(color);
    if (!split_)
        colorL_ = tgt::col4(color);
}

const tgt::col4& TransFuncMappingKey::getColorR() const {
    return colorR_;
}

tgt::col4& TransFuncMappingKey::getColorR() {
    return colorR_;
}

void TransFuncMappingKey::setAlphaR(float a) {
    colorR_.a = static_cast<uint8_t>(a*255.f);
    if (!split_)
        colorL_.a = static_cast<uint8_t>(a*255.f);
}

void TransFuncMappingKey::setAlphaL(float a) {
    colorL_.a = static_cast<uint8_t>(a*255.f);
    if (!split_)
        colorR_.a = static_cast<uint8_t>(a*255.f);
}

float TransFuncMappingKey::getAlphaR() {
    return colorR_.a / 255.f;
}

float TransFuncMappingKey::getAlphaL() {
    return colorL_.a / 255.f;
}

bool TransFuncMappingKey::isSplit() const {
    return split_;
}

void TransFuncMappingKey::setSplit(bool split, bool useLeft) {
    if (split_ == split)
        return;
    if (!split) {
        //join colors
        if (useLeft)
            colorR_ = colorL_;
        else
            colorL_ = colorR_;
    }
    split_ = split;
}

float TransFuncMappingKey::getIntensity() const {
    return intensity_;
}

void TransFuncMappingKey::setIntensity(float i) {
    intensity_ = i;
}

void TransFuncMappingKey::serialize(Serializer& s) const {
    s.serialize("intensity", intensity_);
    s.serialize("split", split_);

    s.serialize("colorL", colorL_);
    if (split_)
        s.serialize("colorR", colorR_);
}

void TransFuncMappingKey::deserialize(Deserializer& s) {
    s.deserialize("intensity", intensity_);
    s.deserialize("split", split_);

    tgt::col4 color;
    s.deserialize("colorL", color);
    setColorL(color);
    if (split_) {
        s.deserialize("colorR", color);
        setColorR(color);
    }
}

TransFuncMappingKey* TransFuncMappingKey::clone() const {
    TransFuncMappingKey* k = new TransFuncMappingKey();
    k->colorL_ = colorL_;
    k->colorR_ = colorR_;
    k->intensity_ = intensity_;
    k->split_ = split_;
    return k;
}


} // namespace voreen
