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

#ifndef VRN_PROPERTY_H
#define VRN_PROPERTY_H

#include "voreen/core/voreenobject.h"
#include "voreen/core/properties/link/propertylink.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/datastructures/callback/callback.h"
#include "voreen/core/datastructures/callback/callbackmanager.h"

#include "tgt/exception.h"
#include "tgt/vector.h"

#include <string>
#include <vector>
#include <set>
#include <sstream>

namespace voreen {

class PropetyOwner;
class PropertyWidget;
class PropertyWidgetFactory;
class ProcessorNetwork;

/**
 * Abstract base class for properties that can be assigned to processors
 * and other PropertyOwners.
 *
 * @see TemplateProperty
 */
class VRN_CORE_API Property : public VoreenSerializableObject {

    friend class PropertyLink;
    friend class PropertyVector;
    friend class Processor;
    friend class ProcessorNetwork;

public:
    /**
     * Every property has a level of detail setting which determines, whether the property
     * should be visible in a network editing mode.
     */
    enum LevelOfDetail {
        LOD_DEBUG = 0,      ///< used in network mode (propertylistwidget)
        LOD_ADVANCED = 1,   ///< used in network mode
        LOD_DEFAULT = 2,    ///< used in network mode

        LOD_DEVELOPMENT = LOD_ADVANCED, ///< used in general settings
        LOD_APPLICATION = LOD_DEFAULT   ///< used in general settings
    };

    /**
     * The property view flags specifies how a property should be represented
     * in the user interface. Multiple view flags can be combined.
     */
    enum ViewFlags {
        VF_VISIBLE = 1,
        VF_READ_ONLY = 1 << 1
    };

    /**
     * Constructor - sets standard-values
     *
     * @param property's identifier: Must be not empty
     *        and must be unique among all properties of a PropertyOwner
     * @param guiName textual representation of the property in the GUI
     * @param invalidationLevel the level the owner is invalidated with upon change
     */
    Property(const std::string& id, const std::string& guiName,
             int invalidationLevel = Processor::INVALID_RESULT, LevelOfDetail lod = LOD_DEFAULT);

    /**
     * Default constructor. Used for serialization. Do not call directly!
     */
    Property();

    virtual ~Property();

    /**
    * Override this method for performing initializations
    * of the property. It is usually called by the owning Processor's
    * initialize() function.
    *
    * @note All OpenGL initializations must be done here,
    *       instead of the constructor! Time-consuming operations
    *       should also happen here.
    *
    * @throw tgt::Exception if the initialization failed
    */
    virtual void initialize();

    /**
    * Override this method for performing deinitializations
    * of the property.
    *
    * @note All OpenGL deinitializations must be done here,
    *       instead of the destructor!
    *
    * @throw tgt::Exception if the deinitialization failed
    */
    virtual void deinitialize();

    /// returns true if the property has been initialized and false if it has not been initialized yet or has already been deinitialized
    virtual bool isInitialized() const;

    /**
     * Resets the property to its default value.
     *
     * This method is expected to be re-implemented by each concrete subclass.
     */
    virtual void reset() = 0;

    /**
     * Sets the GuiName and updates the attached widgets
     * @see VoreenObject
     */
    virtual void setGuiName(const std::string& guiName);

    /**
     * Returns a textual description of the property type,
     * usually corresponding to the type of the stored value.
     */
    virtual std::string getTypeDescription() const;

    /**
     * Returns the InvalidationLevel of this property.
     *
     * @return The owner is invalidated with this InvalidationLevel upon change.
     */
    int getInvalidationLevel() const;

    void setInvalidationLevel(int invalidationLevel);

    /**
     * Returns the identifier of the property preceeded by
     * its owner's name, e.g. "Background.color1".
     *
     * If the property is not assigned to an owner,
     * only its id is returned.
     *
     * @see PropertyOwner::getName
     */
    std::string getFullyQualifiedID() const;

    /**
     * Returns the property's level of detail.
     */
    LevelOfDetail getLevelOfDetail() const;

    /**
     * Returns the GUI name of the property preceeded by
     * its owner's name, e.g. "Background.First Color".
     *
     * If the property is not assigned to an owner,
     * only its GUI name is returned.
     *
     * @see PropertyOwner::getName
     * @see getGuiName
     */
    std::string getFullyQualifiedGuiName() const;

    /**
     * Specifies the visibility of the property in the user interface.
     * Internally passes the visibility state to the assigned widgets.
     */
    void setVisibleFlag(bool state);

    /**
     * Returns whether the property (i.e., its widgets) is visible in the GUI.
     */
    bool isVisibleFlagSet() const;

    /**
     * Sets all widgets of this property to read only or not.
     * @note read only properties do not invalidate their owner
     */
    void setReadOnlyFlag(bool state);

    /**
     * Indicates whether the widgets of this property are read only or not.
     * @note read only properties do not invalidate their owner
     */
    bool isReadOnlyFlagSet() const;

    /**
     * Returns all view flags.
     */
    ViewFlags getViewFlags() const;

    /**
     * Sets all view flags.
     * @note read only properties do not invalidate their owner
     */
    void setViewFlags(ViewFlags flags);

    /**
     * Sets the processor this property is assigned to.
     */
    virtual void setOwner(PropertyOwner* owner);

    /**
     * Returns the processor this property is assigned to.
     */
    PropertyOwner* getOwner() const;

    /**
     * Notifies the property that its stored value has changed.
     */
    virtual void invalidate();

    /**
     * Registers a widget at this property. This widget is considered the property
     * in the user interface and is affected by setVisible and setWidgetsEnabled.
     */
    virtual void addWidget(PropertyWidget* widget);

    /**
     * Unregisters the widget from the property without deleting it.
     */
    void removeWidget(PropertyWidget* widget);

    /**
     * Calls updateFromProperty() on the assigned widgets,
     * which causes them to update their state from the property.
     */
    virtual void updateWidgets();

    /**
     * Unregisters all assigned property widgets.
     */
    void disconnectWidgets();

    /**
     * Returns all property widgets assigned to this property.
     */
    const std::set<PropertyWidget*> getPropertyWidgets() const;

    /**
     * Assigns the property to a property group.
     *
     * The group membership of a property has no effect on the property itself,
     * but may be considered in a GUI representation. The default group ID
     * is the empty string, indicating that the property does not belong to a group.
     */
    void setGroupID(const std::string& id);

    /**
     * Sets the guiname of the group the property is in.
     *
     * The group membership of a property has no effect on the property itself,
     * but may be considered in a GUI representation. The default group name
     * is the empty string. For representation the first not empty string of all group members is used.
     */
    void setGroupName(const std::string& name);

    /**
     * Returns the id of the property group to which the property is assigned to.
     * If the property does not belong to a group, an empty string is returned.
     */
    std::string getGroupID() const;

    /**
     * Returns the name of the property group to which the property is assigned to.
     * If the property does not belong to a group, an empty string is returned.
     */
    std::string getGroupName() const;

    /**
     * Returns the property links currently registered
     * at the property.
     */
    const std::vector<PropertyLink*>& getLinks() const;

    /**
     * Returns the property link from this property to the destination property,
     * or null if no link is established between them.
     */
    PropertyLink* getLink(const Property* dest) const;

    /**
     * Returns whether a link from this property to the destination
     * property is established.
     *
     * @param dest the destination property whose link state is to be checked
     * @param transitive if set to true, not only direct links between source
     *      and destination are considered, but also indirect ones over multiple
     *      intermediate hops.
     *
     * @return true, if the dest property is reachable over links or if
     *      this == dest
     */
    bool isLinkedWith(const Property* dest, bool transitive = false) const;

    /**
     * Returns true, if a link evaluator exists between this and the destination
     * property.
     *
     * @param dst the destination property whose link state is to be checked
     */
    bool isLinkableWith(const Property* dst) const;

    /**
     * Returns the gui name and id of all link evaluators available between this and
     * the destination property.
     *
     * @param dst the destination property whose link state is to be checked
     */
    std::vector<std::pair<std::string, std::string> > getCompatibleEvaluators(const voreen::Property* dst) const;

    /**
     * Returns the meta data container of this processor.
     * External objects, such as GUI widgets, can use it
     * to store and retrieve persistent meta data without
     * having to bother with the serialization themselves.
     *
     * @see MetaDataContainer
     */
    MetaDataContainer& getMetaDataContainer() const;

    /// @see Serializable::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Serializable::deserialize
    virtual void deserialize(Deserializer& s);

    ///Serialize the value of the property without the meta data/LOD/guiName
    virtual void serializeValue(Serializer& s) const;

    ///Deserialize the value of the property without the meta data/LOD/guiName
    virtual void deserializeValue(Deserializer& s);

    std::string getDescription() const;

    /// Sets the description
    void setDescription(std::string desc);

    /**
     * Register a callback
     * @param callback Is called if the data of the porperty changes
     */
    void onChange(const Callback& action);

    /**
     * Blocks or unblocks the callback execution.
     */
    void blockCallbacks(bool block);

protected:

    /**
     * Invalidates the owner with the InvalidationLevel set in the constructor
     * @note read only properties do not invalidate their owner
     */
    void invalidateOwner();

    /**
     * Invalidates the owner with a given InvalidationLevel.
     *
     * @param invalidationLevel Use this InvalidationLevel to invalidate
     */
    void invalidateOwner(int invalidationLevel);

    PropertyOwner* owner_;
    int invalidationLevel_;
    std::string groupId_;
    std::string groupName_;

    /// view related enums
    LevelOfDetail levelOfDetail_;
    ViewFlags viewFlags_;

    std::set<PropertyWidget*> widgets_;
    std::vector<PropertyLink*> links_;
    CallbackManager onChangeCallbacks_;
private:
    ///Executes the links currently registered at the property
    void executeLinks();

    /**
     * Adds the passed property link to the property.
     * Is called by the owning ProcessorNetwork and by
     * the PropertyLink object.
     */
    void registerLink(PropertyLink* link);

    /**
     * Removes the passed link from the property, but does not delete it.
     * Is called by the PropertyLink's destructor.
     */
    void removeLink(PropertyLink* link);

    /// flag for marking the initialization status
    bool isInitialized_;

    /// Used for cycle prevention during interaction mode propagation
    bool interactionModeVisited_;

    /// Used for (de-)serializeValue methods
    mutable bool serializeValue_;

    /// Used for cycle prevention during check whether two props are linked
    mutable bool linkCheckVisited_;

    /**
     * Contains the associated meta data.
     *
     * We want to return a non-const reference to it from a const member function
     * and since the MetaDataContainer does not affect the processor itself,
     * mutable appears justifiable.
     */
    mutable MetaDataContainer metaDataContainer_;

    /**
     * Stores the gui name that has been passed to the constructor.
     * The gui name is only serialized, if it differs from the initial one.
     */
    std::string initialGuiName_;

    /// Description for display in GUI etc.
    std::string description_;

    //

};

    //define enum operators
    inline Property::ViewFlags operator~(Property::ViewFlags a)
    {return static_cast<Property::ViewFlags>(~static_cast<int>(a));}
    inline Property::ViewFlags operator|(Property::ViewFlags a, Property::ViewFlags b)
    {return static_cast<Property::ViewFlags>(static_cast<int>(a) | static_cast<int>(b));}
    inline Property::ViewFlags operator&(Property::ViewFlags a, Property::ViewFlags b)
    {return static_cast<Property::ViewFlags>(static_cast<int>(a) & static_cast<int>(b));}

} // namespace voreen

#endif // VRN_PROPERTY_H
