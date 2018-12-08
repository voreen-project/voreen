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

namespace voreen {

/**
 * This property provides a single list of items which can be configured interactively
 * using a second list storing available items.
 * Items can be transferred from one list to the other and even duplicated in the output list.
 */
class VRN_CORE_API InteractiveListProperty : public Property {
public:

    struct Instance {
        int itemId_;
        int instanceId_;
        std::string name_;
    };

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
     * @return number of times
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
     * The output list will stay untouched.
     * The item must not be added beforehand.
     *
     * @param item Item to be added
     */
    void addItem(const std::string& item);

    /**
     * Removes a single item. This will reset input and output list!
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
     */
    const std::vector<int>& getInputIndices() const;

    /**
     * Returns all instances in the output list.
     * @note in case duplicates are allowed, the items will have the format "Input#1"
     * @return
     */
    const std::vector<Instance>& getInstances() const;

    /**
     * Returns all instances in the output list of a specific item.
     * @note in case duplicates are allowed, the items will have the format "Input#1"
     * @return
     */
    std::vector<Instance> getInstances(const std::string& item) const;

    /**
     * Moves a single item from input to the desired position of output.
     * @note if duplicates are allowed, the item will be copied
     * @note invalid items will be silently ignored
     * @param item item to be moved
     * @param pos position the item should be moved (end, if not specified)
     */
    void addInstance(const std::string& item, int pos = -1);

    /**
     * Removes a single instance.
     * @note invalid instances will be ignored silently
     * @param item instance to be removed
     */
    void removeInstance(const std::string& instanceName);
    void removeInstance(int instanceId);

    /**
     * Moves a single instance to a new position.
     * @param instanceName instance to be moved
     * @param pos position to be moved to
     */
    void moveInstance(const std::string& instanceName, int pos);
    void moveInstance(int instanceId, int pos);

    /**
     * Determines, if the specified item is currently part of output.
     * @param item item to be checked
     * @return true, if the specified item is part of the output, false otherwise
     */
    bool hasInstance(const std::string& item) const;

    /**
     * Returns whether items can be chosen multiple times.
     */
    bool isDuplicationAllowed() const;

    /**
     * Sets if duplication is allowed or not.
     * This will reset input and output lists to their defaults.
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

private:

    int getIndexOfItem(const std::string& item) const;
    std::vector<int> getInstanceIds(const std::string& name) const;
    Instance createInstance(int itemId) const;

    //------------------
    //  Member
    //------------------
    std::vector<std::string> items_;
    std::vector<int> inputItemIds_;
    std::vector<Instance> instances_;

    bool allowDuplication_;
    int selectedInstance_;
};

} // namespace voreen

#endif //VRN_INTERACTIVELISTPROPERTY_H
