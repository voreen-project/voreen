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

#ifndef VRN_PROPERTYOWNER_H
#define VRN_PROPERTYOWNER_H

#include "voreen/core/voreenobject.h"
#include "voreen/core/utils/observer.h"
#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

class Property;
class PropertyOwner;
//class PropertyWidget;

class VRN_CORE_API PropertyOwnerObserver : public Observer {
public:
    virtual void preparePropertyRemoval(Property* property);
    /// used in gui layer to determine new added properties (see propertyownerwidget)
    virtual void propertiesChanged(const PropertyOwner* owner);
    /// called, if a property value of the observed owner has been changed (property can be NULL)
    virtual void propertyValueHasBeenModified(Property* property = 0);
};

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API Observable<PropertyOwnerObserver>;
#endif
class VRN_CORE_API PropertyOwner : public VoreenSerializableObject, public Observable<PropertyOwnerObserver>  {

    friend class Property;

public:
    PropertyOwner();
    PropertyOwner(const std::string& id, const std::string& guiName);

    virtual ~PropertyOwner();

    /**
     * Returns the properties assigned to this owner.
     */
    const std::vector<Property*>& getProperties() const;

    /**
     * Returns the property with the given id. If no such property was found,
     * a null-pointer is returned.
     */
    Property* getProperty(const std::string& id) const;

    /**
     * Returns the property with the given GUI name. If no such property was found,
     * a null-pointer is returned.
     */
    Property* getPropertyByName(const std::string& guiName) const;

    /**
     * Returns all properties matching the specified type T,
     * including subtypes.
     */
    template<class T>
    std::vector<T*> getPropertiesByType() const;

    /**
     * The new InvalidationLevel is max(inv, currentInvalidationLevel).
     * Use setValid() to mark the owner as valid.
     *
     * @note the default invalidation level 1 corresponds to Processor::INVALID_RESULT
     */
    virtual void invalidate(int inv = 1);

    /**
     * Returns the property owner's current invalidation level.
     */
    int getInvalidationLevel() const;

    /**
     * Marks the property owner as valid by setting its invalidation level to 0.
     */
    virtual void setValid();

    /**
     * The standard implementation returns true, if getInvalidationLevel() == 0.
     * Override it for custom behavior.
     */
    virtual bool isValid() const;

    /**
     * Assigns a GUI name to a property group.
     *
     * @note Property groups are defined by assigning a group id to
     *  the properties to be grouped.
     *
     * @note this function must be called, if all properties are already assigned to the group.
     *  Adding a new property to the group will not set the guiname automatically.
     *
     * @see Property::setGroupID
     */
    void setPropertyGroupGuiName(const std::string& groupID, const std::string& name);

    /**
     * Returns the GUI name of the property group with the passed ID.
     * If the specified group does not exist or has no name assigned,
     * an empty string is returned.
     */
    std::string getPropertyGroupGuiName(const std::string& groupID) const;

    /**
     * Modifies the visibility of the specified property group's widget.
     * If no group with the passed groupID exists or the group does not have
     * a widget assigned, the call is ignored.
     *
     * @see setPropertyGroupWidget, setPropertyGroupGuiName
     */
    void setPropertyGroupVisible(const std::string& groupID, bool visible);

    /**
     * Resets all properties to their default values
     */
    virtual void resetAllProperties();

    /// @see Serializable::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Serializable::deserialize
    virtual void deserialize(Deserializer& s);

protected:
    /**
     * Registers a property at the owner.
     *
     * @note The property owner does not take ownership of the passed property,
     *  i.e., deleting the property is up to the caller.
     */
    virtual void addProperty(Property* prop);

    /// \overload
    virtual void addProperty(Property& prop);

    /**
     * Unregisters a property at the owner.
     *
     * @note The unregistered property is not deleted.
     */
    virtual void removeProperty(Property* prop);

    /// \overload
    virtual void removeProperty(Property& prop);

    virtual void notifyPropertiesChanged() const;
    virtual void notifyPropertyValueHasBeenModified(Property* prop) const;


    int invalidationLevel_;

    /// category used in logging
    static const std::string loggerCat_;

private:
    /// Stores the owner's properties.
    std::vector<Property*> properties_;
};

//---------------------------------------------------------------------------
// template definitions

template<class T>
std::vector<T*> PropertyOwner::getPropertiesByType() const {
    std::vector<T*> result;
    for (size_t i=0; i<properties_.size(); ++i) {
        if (dynamic_cast<T*>(properties_[i]))
            result.push_back(dynamic_cast<T*>(properties_[i]));
    }
    return result;
}

bool VRN_CORE_API serializeSettings(const PropertyOwner* po, const std::string& filename);
bool VRN_CORE_API deserializeSettings(PropertyOwner* po, const std::string& filename);

} // namespace voreen

#endif // VRN_PROPERTYOWNER_H
