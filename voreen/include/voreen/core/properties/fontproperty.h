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

#ifndef VRN_FONTPROPERTY_H
#define VRN_FONTPROPERTY_H

#include "voreen/core/properties/templateproperty.h"
#include "tgt/font.h"

namespace voreen {

/**
 * Property storing a font.
 * Note: The property is owner of the font object!
 */
class FontProperty : public TemplateProperty<tgt::Font*> {
    using T = tgt::Font*;
public:
    FontProperty(const std::string& id, const std::string& guiText, tgt::Font* value = 0,
                 int invalidationLevel = Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    FontProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FontProperty"; }
    virtual std::string getTypeDescription() const { return "Font"; }

    /**
     * Frees the font object represented by this property.
     */
    virtual ~FontProperty();

    /**
     * Sets a new font object.
     * Note: Frees the former font object represented by this property.
     */
    virtual void set(const T& font) override;

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

};

} // namespace

#endif // VRN_FONTPROPERTY_H
