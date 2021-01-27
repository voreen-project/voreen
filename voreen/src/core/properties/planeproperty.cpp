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

#include "voreen/core/properties/planeproperty.h"

namespace voreen {

PlaneProperty::PlaneProperty(const std::string& id, const std::string& guiText, const tgt::Plane<float>& value,
                             int invalidationLevel, Property::LevelOfDetail lod)
                             : TemplateProperty<tgt::Plane<float> >(id, guiText, value, invalidationLevel,lod)
{}

void PlaneProperty::serialize(Serializer& s) const {
    TemplateProperty<tgt::Plane<float> >::serialize(s);

    s.serialize("value", value_.toVec4());
}

void PlaneProperty::deserialize(Deserializer& s) {
    TemplateProperty<tgt::Plane<float> >::deserialize(s);
    try {
        tgt::vec4 tmpPlane;
        s.deserialize("value", tmpPlane);
        value_ = tgt::Plane<float>(tmpPlane.x,tmpPlane.y,tmpPlane.z,tmpPlane.a);
    }
    catch(SerializationException&) {
        value_ = tgt::Plane<float>(1.f,0.f,0.f,1.f);
        s.removeLastError();
    }
}

} // namespace voreen
