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

#include "voreen/core/properties/fontproperty.h"

#include "voreen/core/properties/condition.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/filesystem.h"

namespace voreen {

FontProperty::FontProperty(const std::string& id, const std::string& guiText, tgt::Font* value,
                           int invalidationLevel, Property::LevelOfDetail lod)
                               : TemplateProperty<tgt::Font*>(id, guiText, value, invalidationLevel, lod)
{
    if (!value)
        value_ = new tgt::Font(VoreenApplication::app()->getFontPath("VeraMono.ttf"), 12, 0);
}

FontProperty::FontProperty() {
    value_ = 0;
}

FontProperty::~FontProperty() {
    delete value_;
}

Property* FontProperty::create() const {
    return new FontProperty();
}

void FontProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    tgtAssert(value_, "no font object");


    s.serialize("fontSize", value_->getFontSize());
    s.serialize("fontName", tgt::FileSystem::fileName(value_->getFontName()));
    s.serialize("textAlignment", static_cast<int>(value_->getTextAlignment()));
    s.serialize("lineWidth", value_->getLineWidth());
    s.serialize("fontColor", value_->getFontColor());
}

void FontProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    int textAlignment;
    float lineWidth;
    tgt::vec4 fontColor;
    std::string fontName;
    int fontSize;


    s.optionalDeserialize("fontSize", fontSize, 10);
    s.optionalDeserialize("fontName", fontName, std::string("VeraMono.ttf"));
    s.optionalDeserialize("textAlignment", textAlignment, static_cast<int>(tgt::Font::BottomLeft));
    s.optionalDeserialize("lineWidth", lineWidth, 0.0f);
    s.optionalDeserialize("fontColor", fontColor, tgt::vec4(1.0f));

    if (!fontName.empty()) {
        delete value_;
        set(new tgt::Font(VoreenApplication::app()->getFontPath(fontName), fontSize, lineWidth,
            static_cast<tgt::Font::TextAlignment>(textAlignment)));
        get()->setFontColor(fontColor);

    }
}

} // namespace
