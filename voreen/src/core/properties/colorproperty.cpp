/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/properties/colorproperty.h"

namespace voreen {

ColorProperty::ColorProperty(const std::string& id, const std::string& guiText,
                     tgt::vec4 value, int invalidationLevel, Property::LevelOfDetail lod,
                     bool useAlphaChannel)
    : FloatVec4Property(id, guiText, value, tgt::vec4(0.f), tgt::vec4(1.f), invalidationLevel, NumericProperty<tgt::vec4>::STATIC, lod)
    , useAlphaChannel_(useAlphaChannel)
{}

ColorProperty::ColorProperty()
    : FloatVec4Property("", "", tgt::vec4(1.f), tgt::vec4(0.f), tgt::vec4(1.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec4>::STATIC)
    , useAlphaChannel_(true)
{}

void ColorProperty::set(const tgt::vec4& value) {
    //set alpha to 1.f, if alpha channel is disabled
    if(!useAlphaChannel_) {
        tgt::vec4 modified(value.r,value.g,value.b,1.f);
        FloatVec4Property::set(modified);
    } else {
        FloatVec4Property::set(value);
    }
}

tgt::col4 ColorProperty::getRGBA() const {
    tgt::vec4 value = get();
    //transform the normalized value into uint_8
    return tgt::col4(static_cast<uint8_t>(value.r*255),
                     static_cast<uint8_t>(value.g*255),
                     static_cast<uint8_t>(value.b*255),
                     static_cast<uint8_t>(value.a*255));
}

}   // namespace
