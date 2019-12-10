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

#ifndef VRN_LINKEVALUATOR_H
#define VRN_LINKEVALUATOR_H

/*
 * ATTENTION: Do not include voreen/core/properties/property.h due to circular
 *            header-file inclusion. Keep the Property forward declaration instead.
 */
#include "voreen/core/voreenobject.h"
#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

class Property;

/**
 * An object of this type is assigned to each PropertyLink, which uses
 * it for the actual link execution.
 */
class VRN_CORE_API LinkEvaluatorBase : public VoreenSerializableObject {
public:
    virtual ~LinkEvaluatorBase() {}

    virtual std::string getGuiName() const = 0;

    ///Called by PropertyLink for executing the link.
    virtual void eval(Property* src, Property* dst) = 0;

    virtual void propertiesChanged(Property* src, Property* dst);

    //Returns true if the LinkEvaluator can link the two properties.
    virtual bool arePropertiesLinkable(const Property* src, const Property* dst) const = 0;

    /// @see Serializable::serialize
    virtual void serialize(Serializer&) const {}

    /// @see Serializable::deserialize
    virtual void deserialize(Deserializer&) {}
};

} // namespace

#endif // VRN_LINKEVALUATOR_H
