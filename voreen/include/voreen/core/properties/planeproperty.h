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

#ifndef VRN_PLANEPROPERTY_H
#define VRN_PLANEPROPERTY_H

#include "voreen/core/properties/templateproperty.h"

#include "tgt/plane.h"

namespace voreen {

class VRN_CORE_API PlaneProperty : public TemplateProperty<tgt::Plane<float> > {

public:

    PlaneProperty() {}

    PlaneProperty(const std::string& id, const std::string& guiText, const tgt::Plane<float>& value = tgt::Plane<float>(1.f,0.f,0.f,1.f),
                  int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);

    virtual Property* create() const {
        return new PlaneProperty();
    }

    virtual std::string getClassName() const       { return "PlaneProperty"; }
    virtual std::string getTypeDescription() const { return "Property storing a tgt::Plane"; }

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

};

} // namespace voreen

#endif // VRN_PLANEPROPERTY_H
