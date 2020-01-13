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

#include "voreen/core/properties/color/colorswitchproperty.h"

namespace voreen {

ColorSwitchProperty::ColorSwitchProperty(const std::string& id, const std::string& guiText,
                     tgt::vec4 activeColor, tgt::vec4 inactiveColor, int invalidationLevel, Property::LevelOfDetail lod)
    : ColorProperty(id, guiText, activeColor, invalidationLevel, lod)
    , activeColor_(activeColor)
    , inactiveColor_(inactiveColor)
    , useActiveColor_(true)
{}

ColorSwitchProperty::ColorSwitchProperty()
    : ColorProperty()
{}


Property* ColorSwitchProperty::create() const {
    return new ColorSwitchProperty();
}


void ColorSwitchProperty::serialize(Serializer& s) const
{
    ColorProperty::serialize(s);
    s.serialize("activeColor", activeColor_);
    s.serialize("inactiveColor", inactiveColor_);
    s.serialize("useActiveColor", useActiveColor_);
}

void ColorSwitchProperty::deserialize(Deserializer& s)
{
    ColorProperty::deserialize(s);
    s.deserialize("activeColor", activeColor_);
    s.deserialize("inactiveColor", inactiveColor_);
    s.deserialize("useActiveColor", useActiveColor_);
}

const tgt::vec4& ColorSwitchProperty::get() const
{
    if (useActiveColor_)
        return getActiveColor();
    else
        return getInactiveColor();
}

void ColorSwitchProperty::set(const tgt::vec4& value)
{
    if (useActiveColor_)
        setActiveColor(value);
    else
        setInactiveColor(value);
}

const tgt::vec4& ColorSwitchProperty::getActiveColor() const
{
    return activeColor_;
}

const tgt::vec4& ColorSwitchProperty::getInactiveColor() const
{
    return inactiveColor_;
}

void ColorSwitchProperty::setActiveColor(tgt::vec4 color)
{
    if (color != activeColor_){
        activeColor_ = color;
        if (useActiveColor_)
            value_ = color;
        invalidate();
    }
}

void ColorSwitchProperty::setInactiveColor(tgt::vec4 color)
{
    if (color != inactiveColor_){
        inactiveColor_ = color;
        if (!useActiveColor_)
            value_ = color;
        invalidate();
    }
}

tgt::col4 ColorSwitchProperty::getActiveColorRGBA() const
{
    return tgt::col4(static_cast<uint8_t>(activeColor_.r*255),
                     static_cast<uint8_t>(activeColor_.g*255),
                     static_cast<uint8_t>(activeColor_.b*255),
                     static_cast<uint8_t>(activeColor_.a*255));
}

tgt::col4 ColorSwitchProperty::getInactiveColorRGBA() const
{
    return tgt::col4(static_cast<uint8_t>(inactiveColor_.r*255),
                     static_cast<uint8_t>(inactiveColor_.g*255),
                     static_cast<uint8_t>(inactiveColor_.b*255),
                     static_cast<uint8_t>(inactiveColor_.a*255));
}

bool ColorSwitchProperty::getUseActiveColor() const
{
    return useActiveColor_;
}

void ColorSwitchProperty::setUseActiveColor(bool use)
{
    if (use != useActiveColor_){
        useActiveColor_ = use;
        if (use)
            value_ = activeColor_;
        else
            value_ = inactiveColor_;
        invalidate();
    }
}

}   // namespace
