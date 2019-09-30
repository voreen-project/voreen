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

#ifndef VRN_INTERACTIVELISTPROPERTY_H
#define VRN_INTERACTIVELISTPROPERTY_H

#include "voreen/core/properties/property.h"

#include <vector>
#include <string>
#include <functional>

namespace voreen {

/**
 * This property provides a list of items which can be configured interactively
 * using a second list storing available items.
 * Items can be transferred from one list to the other and vice-versa.
 */
class VRN_CORE_API InteractiveListProperty : public Property {
public:

    class VRN_CORE_API Instance : public Serializable {
    public:
        Instance();
        Instance(int itemId, int instanceId);

        virtual void serialize(Serializer& s) const;
        virtual void deserialize(Deserializer& s);

        int getItemId() const;
        int getInstanceId() const;

        bool isActive() const;
        void setActive(bool active);

        const std::string& getName() const;
        void setName(const std::string& name);

    private:

        int itemId_;            ///< id of associated item
        int instanceId_;        ///< unique instance id

        bool active_;           ///< active state
        std::string name_;      ///< instance name
    };

    typedef std::function<std::string(const Instance&)> NameGenerator;

    InteractiveListProperty(const std::string& id, const std::string& guiText, bool allowDuplication = false,
                        int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    InteractiveListProperty();
    virtual ~InteractiveListProperty();

    virtual Property* create() const               { return new InteractiveListProperty(); }
    virtual std::string getClassName() const       { return "InteractiveListProperty"; }
    virtual std::string getTypeDescription() const { return "List"; }
    /** @see Property::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Property::deserialize */
    virtual void deserialize(Deserializer& s);

    /**
     * Removes all output items.
     */
    virtual void reset();

    /**
     * Removes all contained items.
     */
    void clear();

    /**
     * Returns the number of items.
     *
     * @return number of itmes
     */
    size_t getNumItems() const;

    /**
     * Sets (overwrites!) all items of the input list.
     * This will also clear the output list.
     *
     * @param items Items to be set
     */
    void setItems(const std::vector<std::string>& items);

    /**
     * Adds a single item to the input list.
     * The instance list will stay untouched.
     * The item must not have been added beforehand.
     *
     * @param item Item to be added
     */
    void addItem(const std::string& item);

    /**
     * Removes a single item. This will reset input and instance list and the selection!
     * @note if the item is not available, the call will be silently ignored.
     *
     * @param item Item to be removed
     */
    void removeItem(const std::string& item);

    /**
     * Returns all contained items in the order of insertion.
     *
     * @return all contained items
     */
    const std::vector<std::string>& getItems() const;

    /**
     * Returns indices of the available items.
     * The order is implied by the returned list of getItems().
     * In case duplicates are allowed, this list always contains
     * all indices in ascending order.
     *
     * @return indices of available items
     */
    std::vector<int> getInputIndices() const;

    /**
     * Returns all instances of the item list.
     * @return all instances of the item list
     */
    const std::vector<Instance>& getInstances() const;

    /**
     * Returns all instances of the item list.
     * @note In order to apply the changes, the property has to be invalidated first!
     * @return all instances of the item list
     */
    std::vector<Instance>& getInstances();

    /**
     * Returns all instances in the output list of a specific item.
     * @param item item of which all instances should be returned
     * @return all instances of the specified item
     */
    std::vector<Instance> getInstances(const std::string& item) const;

    /**
     * Moves a single item from input to the desired position of the instance list.
     * @note if duplicates are allowed, the item will be copied
     * @note invalid items will be silently ignored
     * @note selection will be updated accordingly
     * @param item item to be moved
     * @param pos position the item should be moved (end, if not specified)
     */
    void addInstance(const std::string& item, int pos = -1);

    /**
     * Removes a single instance.
     * @note invalid instances will be ignored silently
     * @note selection will be update accordingly
     * @param instanceId instance to be removed
     */
    void removeInstance(int instanceId);

    /**
     * Moves a single instance to a new position.
     * @note selection will be update accordingly
     * @param instanceId instance to be moved
     * @param pos position to be moved to
     */
    void moveInstance(int instanceId, int pos);

    /**
     * Moves a single instance to a new position by swapping
     * positions with a target instance.
     * @note selection will be update accordingly
     * @param instanceId instance to be swapped
     * @param pos position to be moved to
     */
    void swapInstances(int instanceId, int pos);

    /**
     * Determines, if the specified item has at least one instance.
     * @param item item to be checked
     * @return true, if the specified item is part of the output, false otherwise
     */
    bool hasInstance(const std::string& item) const;

    /**
     * Returns whether items can have multiple instances.
     */
    bool isDuplicationAllowed() const;

    /**
     * Sets if duplication is allowed or not.
     * This will reset input and instance lists to their defaults.
     */
    void setDuplicationAllowed(bool enabled);

    /**
     * Determines the selected index of the instance list.
     * @note the index refers to the instance list, NOT to the items' indices
     * @note not to be confused with instanceId
     * @return the selected index of the of the instance list
     */
    int getSelectedInstance() const;

    /**
     * Sets the selected index in the instance list.
     * @param index index to be selected
     */
    void setSelectedInstance(int index);

    /**
     * Sets the name generator for instance name generation.
     */
    void setNameGenerator(const NameGenerator& nameGenerator);

    /**
     * Returns the name generator for instance name generation.
     */
    const NameGenerator& getNameGenerator() const;

private:

    int getIndexOfItem(const std::string& item) const;
    int getIndexOfInstance(int instanceId) const;
    std::vector<int> getInstanceIds(const std::string& name) const;
    Instance createInstance(int itemId) const;

    //------------------
    //  Member
    //------------------
    NameGenerator nameGenerator_;
    std::vector<std::string> items_;
    std::vector<Instance> instances_;

    bool allowDuplication_;
    int selectedInstance_;
};

} // namespace voreen

#endif //VRN_INTERACTIVELISTPROPERTY_H
