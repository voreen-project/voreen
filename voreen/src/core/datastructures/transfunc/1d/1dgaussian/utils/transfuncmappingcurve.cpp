/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/transfunc/1d/1dgaussian/utils/transfuncmappingcurve.h"
#include "tgt/assert.h"

namespace voreen {

TransFuncMappingCurve::TransFuncMappingCurve(float i, float width, float baseValue, const tgt::col4& color)
    : intensity_(i)
    , width_(width)
    , baseValue_(baseValue)
    , colorL_(color)
    , colorR_(color)
    , baseColorL_(color)
    , baseColorR_(color)
    , split_(false)
    , unicolor_(true)
    , active_(true)
{}

TransFuncMappingCurve::TransFuncMappingCurve()
    : intensity_(0.5f)
    , width_(0.01f)
    , baseValue_(0)
    , colorL_(tgt::vec4(0.f))
    , colorR_(tgt::vec4(0.f))
    , baseColorL_(tgt::vec4(0.f))
    , baseColorR_(tgt::vec4(0.f))
    , split_(false)
    , unicolor_(true)
    , active_(true)
{}

TransFuncMappingCurve::~TransFuncMappingCurve() {
}

bool TransFuncMappingCurve::operator==(const TransFuncMappingCurve& curve) {
    return (intensity_ == curve.intensity_) && (split_ == curve.split_) &&
           (colorR_    == curve.colorR_) && (colorL_ == curve.colorL_) &&
           (width_ == curve.width_) && (baseValue_ == curve.baseValue_) &&
           (unicolor_ == curve.unicolor_) && (active_ == curve.active_);
}

bool TransFuncMappingCurve::operator!=(const TransFuncMappingCurve& curve) {
    return !(*this == curve);
}


/// ---------------------------------------------------------  PEAK Color Getter and Setter ///
void TransFuncMappingCurve::setColorL(const tgt::col4& color) {
    colorL_ = color;
    if (unicolor_) {
        // set rgb of all colors
        setAllColors(color);
    }
    if (!split_)
        colorR_ = color;
}

void TransFuncMappingCurve::setColorL(const tgt::ivec4& color) {
    colorL_ = tgt::col4(color);
    if (unicolor_) {
        // set rgb of all colors
        setAllColors(colorL_);
    }
    if (!split_)
        colorL_ = tgt::col4(color);
}

const tgt::col4& TransFuncMappingCurve::getColorL() const {
    return colorL_;
}

tgt::col4& TransFuncMappingCurve::getColorL() {
    return colorL_;
}

void TransFuncMappingCurve::setColorR(const tgt::col4& color) {
    colorR_ = color;
    if (unicolor_) {
        // set rgb of all colors
        setAllColors(color);
    }
    if (!split_)
        colorL_ = color;
}

void TransFuncMappingCurve::setColorR(const tgt::ivec4& color) {
    colorR_ = tgt::col4(color);
    if (unicolor_) {
        // set rgb of all colors
        setAllColors(colorR_);
    }
    if (!split_)
        colorL_ = colorR_;
}

const tgt::col4& TransFuncMappingCurve::getColorR() const {
    return colorR_;
}

tgt::col4& TransFuncMappingCurve::getColorR() {
    return colorR_;
}

void TransFuncMappingCurve::setAlphaR(float a) {
    colorR_.a = static_cast<uint8_t>(a*255.f);
    if (!split_)
        colorL_.a = static_cast<uint8_t>(a*255.f);
}

void TransFuncMappingCurve::setAlphaL(float a) {
    colorL_.a = static_cast<uint8_t>(a*255.f);
    if (!split_)
        colorR_.a = static_cast<uint8_t>(a*255.f);
}

float TransFuncMappingCurve::getAlphaR() {
    return colorR_.a / 255.f;
}

float TransFuncMappingCurve::getAlphaL() {
    return colorL_.a / 255.f;
}

/// -------------------------------------------------  BASE Color Getter and Setter ///
// the alpha values of the base colors are irrelevant and will be ignored in all calculations

void TransFuncMappingCurve::setBaseColorL(const tgt::col4& color) {
    if (unicolor_)
        setAllColors(color);
    else
        baseColorL_ = color;
}

void TransFuncMappingCurve::setBaseColorL(const tgt::ivec4& color) {
    if (unicolor_)
        setAllColors(tgt::col4(color));
    else
        baseColorL_ = tgt::col4(color);
}

const tgt::col4& TransFuncMappingCurve::getBaseColorL() const {
    return baseColorL_;
}

tgt::col4& TransFuncMappingCurve::getBaseColorL() {
    return baseColorL_;
}

void TransFuncMappingCurve::setBaseColorR(const tgt::col4& color) {
    if (unicolor_)
        setAllColors(color);
    else
        baseColorR_ = color;
}

void TransFuncMappingCurve::setBaseColorR(const tgt::ivec4& color) {
    if (unicolor_)
        setAllColors(tgt::col4(color));
    else
        baseColorR_ = tgt::col4(color);
}

const tgt::col4& TransFuncMappingCurve::getBaseColorR() const {
    return baseColorR_;
}

tgt::col4& TransFuncMappingCurve::getBaseColorR() {
    return baseColorR_;
}


/// --------------------------------------------------------- other Getters and Setters ///

void TransFuncMappingCurve::setAllColors(const tgt::col4& color, bool setAlpha) {
    if (!setAlpha) {
        baseColorL_ = tgt::col4(color.r, color.g, color.b, baseColorL_.a);
        baseColorR_ = tgt::col4(color.r, color.g, color.b, baseColorR_.a);
        colorL_ = tgt::col4(color.r, color.g, color.b, colorL_.a);
        colorR_ = tgt::col4(color.r, color.g, color.b, colorR_.a);
    }
    else {
        baseColorL_ = color;
        baseColorR_ = color;
        colorL_ = color;
        colorR_ = color;
    }
}

bool TransFuncMappingCurve::isSplit() const {
    return split_;
}

void TransFuncMappingCurve::setSplit(bool split, bool useLeft) {
    if (split_ == split)
        return;

    if (!split) {
        // join colors
        if (useLeft)
            colorR_ = colorL_;
        else
            colorL_ = colorR_;
    }
    split_ = split;
}

float TransFuncMappingCurve::getIntensity() const {
    return intensity_;
}

void TransFuncMappingCurve::setIntensity(float i) {
    intensity_ = i;
}

float TransFuncMappingCurve::getWidth() const {
    return width_;
}

void TransFuncMappingCurve::setWidth(float w) {
    tgtAssert(w > 0, "curve width not positve")
    width_ = w;
}

float TransFuncMappingCurve::getBaseValue() const {
    return baseValue_;
}

void TransFuncMappingCurve::setBaseValue(float b) {
    baseValue_ = b;
}

float TransFuncMappingCurve::getOpacityAt(float intensity) {
    if (intensity < intensity_)
        return exp(-(intensity - intensity_)*(intensity - intensity_) / (2 * width_)) * (getAlphaL() - baseValue_) + baseValue_;
    else
        return exp(-(intensity - intensity_)*(intensity - intensity_) / (2 * width_)) * (getAlphaR() - baseValue_) + baseValue_;
}

tgt::col4 TransFuncMappingCurve::getColorAt(float intensity) {

    // peak to base color ratio is based on the Gauss function without the factors for base and peak alpha values
    float yRatio = exp(-(intensity - intensity_)*(intensity - intensity_) / (2 * width_));
    tgt::col4 resultColor;

    // determine the rgb-color values
    // left side of the curve
    if (intensity < intensity_) {
        resultColor.r = static_cast<uint8_t>((colorL_.r * yRatio + baseColorL_.r * (1.f - yRatio)));
        resultColor.g = static_cast<uint8_t>((colorL_.g * yRatio + baseColorL_.g * (1.f - yRatio)));
        resultColor.b = static_cast<uint8_t>((colorL_.b * yRatio + baseColorL_.b * (1.f - yRatio)));
    }
    // right side of the curve
    else {
        resultColor.r = static_cast<uint8_t>((colorR_.r * yRatio + baseColorR_.r * (1.f - yRatio)));
        resultColor.g = static_cast<uint8_t>((colorR_.g * yRatio + baseColorR_.g * (1.f - yRatio)));
        resultColor.b = static_cast<uint8_t>((colorR_.b * yRatio + baseColorR_.b * (1.f - yRatio)));
    }

    // apply this curves's function to the alpha value
    resultColor.a = static_cast<uint8_t>(getOpacityAt(intensity) * 255);

    return resultColor;
}

void TransFuncMappingCurve::setUnicolor(bool unicolor, CurvePart relevantKey) {
    if (unicolor_ == unicolor)
        return;

    if (unicolor) {
        // Set all color values to the color of the relevant key,
        // but without changing the alpha values.
        switch (relevantKey) {
        case CurvePart::PEAK_LEFT:
            setAllColors(colorL_, false);
        case CurvePart::PEAK_RIGHT:
            setAllColors(colorR_, false);
        case CurvePart::BASE_LEFT:
            setAllColors(baseColorL_, false);
        case CurvePart::BASE_RIGHT:
            setAllColors(baseColorR_, false);
        default:
            break;
        }
    }
    unicolor_ = unicolor;
}

bool TransFuncMappingCurve::isUnicolor() {
    return unicolor_;
}

void TransFuncMappingCurve::setActive(bool active) {
    active_ = active;
}

bool TransFuncMappingCurve::isActive() {
    return active_;
}


/// ----------------------------------------------------------- Serialize and Clone ///

void TransFuncMappingCurve::serialize(Serializer& s) const {
    s.serialize("intensity", intensity_);
    s.serialize("width", width_);
    s.serialize("baseValue", baseValue_);
    s.serialize("split", split_);
    s.serialize("unicolor", unicolor_);
    s.serialize("active", active_);

    s.serialize("colorL", colorL_);
    if (split_)
        s.serialize("colorR", colorR_);
    if (!unicolor_) {
        s.serialize("baseColorL", baseColorL_);
        s.serialize("baseColorR", baseColorR_);
    }
}

void TransFuncMappingCurve::deserialize(Deserializer& s) {
    s.deserialize("intensity", intensity_);
    s.deserialize("width", width_);
    s.deserialize("baseValue", baseValue_);
    s.deserialize("split", split_);
    s.deserialize("unicolor", unicolor_);
    s.deserialize("active", active_);

    tgt::col4 color;
    s.deserialize("colorL", color);
    setColorL(color);
    if (split_) {
        s.deserialize("colorR", color);
        setColorR(color);
    }
    if (!unicolor_) {
        s.deserialize("baseColorL", baseColorL_);
        s.deserialize("baseColorR", baseColorR_);
    }
}

TransFuncMappingCurve* TransFuncMappingCurve::clone() const {
    TransFuncMappingCurve* k = new TransFuncMappingCurve();
    k->colorL_ = colorL_;
    k->colorR_ = colorR_;
    k->intensity_ = intensity_;
    k->width_ = width_;
    k->baseValue_ = baseValue_;
    k->baseColorL_ = baseColorL_;
    k->baseColorR_ = baseColorR_;
    k->split_ = split_;
    k->unicolor_ = unicolor_;
    k->active_ = active_;
    return k;
}


} // namespace voreen
