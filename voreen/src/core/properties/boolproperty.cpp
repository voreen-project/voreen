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

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/condition.h"

namespace voreen {

BoolProperty::BoolProperty(const std::string& id, const std::string& guiText, bool value,
                   int invalidationLevel,  Property::LevelOfDetail lod)
    : TemplateProperty<bool>(id, guiText, value, invalidationLevel, lod)
{}

BoolProperty::BoolProperty() :
    TemplateProperty<bool>()
{}

Property* BoolProperty::create() const {
    return new BoolProperty();
}

void BoolProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("value", value_);
}

void BoolProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    bool value;
    s.deserialize("value", value);
    try {
        set(value);
    }
    catch (Condition::ValidationFailed& e) {
        s.addError(e);
    }
}


}   // namespace
