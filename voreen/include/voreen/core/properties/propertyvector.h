/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_PROPERTYVECTOR_H
#define VRN_PROPERTYVECTOR_H

#include "tgt/vector.h"
#include "tgt/logmanager.h"
#include "voreen/core/properties/propertyowner.h"
#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/properties/condition.h"

#include <vector>
#include <string>

namespace voreen {

class VRN_CORE_API PropertyVector : public Property, PropertyOwner {
public:

    PropertyVector(const std::string& id, const std::string& guiText, std::vector<Property*> properties);
    PropertyVector();
    virtual ~PropertyVector();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "PropertyVector"; }
    virtual std::string getTypeDescription() const { return "PropertyVector"; }

    virtual void reset();

    virtual std::string getGuiName() const;
    ///equals the guiname
    virtual std::string getID() const;
    ///not implemented
    virtual void setGuiName(const std::string& guiname);

    void addProperty(Property* prop);

    const std::vector<Property*>& getProperties() const;

    template<typename T>
    T getProperty(int id) const {
        T prop = 0;
        if (id >= 0 && id < size()) {
            prop = dynamic_cast<T>(properties_[id]);
            if (!prop) {
                LWARNING("Property is of invalid type");
            }
        }
        else {
            tgtAssert(false, "Invalid property id.");
            LERROR("Invalid property id: " << id);
        }
        return prop;
    }

    int size() const;

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

    virtual void setOwner(PropertyOwner* owner);

protected:

    /**
     * Initializes the component properties.
     *
     * @see Property::initialize
     */
    void initialize();

    /**
     * Deinitializes the component properties.
     *
     * @see Property::deinitialize
     */
    void deinitialize();

    std::vector<Property*> properties_;

    static const std::string loggerCat_;

};

}   // namespace

#endif
