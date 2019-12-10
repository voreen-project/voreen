/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_COLORSWITCHPROPERTY_H
#define VRN_COLORSWITCHPROPERTY_H

#include "voreen/core/properties/colorproperty.h"

namespace voreen {

    /**
     * Used to represent colors.
     * The color is represented in the rbga space between 0.f - 1.f.
     */
class VRN_CORE_API ColorSwitchProperty : public ColorProperty {
public:
    ColorSwitchProperty(const std::string& id, const std::string& guiText,
                tgt::vec4 activeColor = tgt::vec4(1.f), tgt::vec4 inactiveColor = tgt::vec4(1.f), int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    ColorSwitchProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "ColorSwitchProperty"; }
    virtual std::string getTypeDescription() const { return "FloatVector4"; }

    const tgt::vec4& getActiveColor() const;
    const tgt::vec4& getInactiveColor() const;

    void setActiveColor(tgt::vec4 color);
    void setInactiveColor(tgt::vec4 color);

    tgt::col4 getActiveColorRGBA() const;
    tgt::col4 getInactiveColorRGBA() const;

    bool getUseActiveColor() const;
    void setUseActiveColor(bool use);

    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);

    virtual const tgt::vec4& get() const;

    virtual void set(const tgt::vec4& value);

private:
    tgt::vec4 activeColor_;
    tgt::vec4 inactiveColor_;
    bool useActiveColor_;


};

}   // namespace

#endif
