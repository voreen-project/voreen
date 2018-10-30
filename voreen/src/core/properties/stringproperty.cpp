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

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/condition.h"

#include <sstream>

namespace voreen {

StringProperty::StringProperty(const std::string& id, const std::string& guiText,
                       const std::string& value, int invalidationLevel, Property::LevelOfDetail lod)
    : TemplateProperty<std::string>(id, guiText, value, invalidationLevel, lod)
    , instantUpdate_(true)
{
}

StringProperty::StringProperty()
{}

void StringProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("value", value_);
}

void StringProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    std::string value;
    s.deserialize("value", value);
    try {
        set(value);
    }
    catch (Condition::ValidationFailed& e) {
        s.addError(e);
    }
}

Property* StringProperty::create() const {
    return new StringProperty();
}

void StringProperty::setReadOnly(bool readOnly) {
    setReadOnlyFlag(readOnly);
    updateWidgets();
}

bool StringProperty::isReadOnly() const {
    return isReadOnlyFlagSet();
}

void StringProperty::setInstantUpdate(bool instantUpdate)
{
    if (instantUpdate != instantUpdate_){
        instantUpdate_ = instantUpdate;
        invalidate();
    }
}

bool StringProperty::getInstantUpdate() const
{
    return instantUpdate_;
}

}   // namespace
