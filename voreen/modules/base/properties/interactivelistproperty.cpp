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

#include "interactivelistproperty.h"

namespace voreen {

InteractiveListProperty::Instance::Instance()
    : itemId_(-1)
    , instanceId_(-1)
    , active_(false)
{
}

InteractiveListProperty::Instance::Instance(int itemId, int instanceId)
    : itemId_(itemId)
    , instanceId_(instanceId)
    , active_(true)
{
}

void InteractiveListProperty::Instance::serialize(Serializer& s) const {
    s.serialize("itemId", itemId_);
    s.serialize("instanceId", instanceId_);
    s.serialize("active", active_);
    s.serialize("name", name_);
}

void InteractiveListProperty::Instance::deserialize(Deserializer& s) {
    s.deserialize("itemId", itemId_);
    s.deserialize("instanceId", instanceId_);
    s.deserialize("active", active_);
    s.deserialize("name", name_);
}

int InteractiveListProperty::Instance::getItemId() const {
    return itemId_;
}
int InteractiveListProperty::Instance::getInstanceId() const {
    return instanceId_;
}

bool InteractiveListProperty::Instance::isActive() const {
    return active_;
}
void InteractiveListProperty::Instance::setActive(bool active) {
    active_ = active;
}

const std::string& InteractiveListProperty::Instance::getName() const {
    return name_;
}
void InteractiveListProperty::Instance::setName(const std::string& name) {
    name_ = name;
}


InteractiveListProperty::InteractiveListProperty(const std::string& id, const std::string& guiText, bool allowDuplication,
                            int invalidationLevel, Property::LevelOfDetail lod, bool serializeItems)
    : Property(id, guiText, invalidationLevel, lod)
    , allowDuplication_(allowDuplication)
    , selectedInstance_(-1)
    , serializeItems_(serializeItems)
{
    // Setup default name generator.
    nameGenerator_ =
        [this] (const Instance& instance) {
            std::string name = items_[instance.getItemId()];

            if (allowDuplication_) {
                name += " (" + std::to_string(instance.getInstanceId()) + ")";
            }

            return name;
        };
}

InteractiveListProperty::InteractiveListProperty()
{
}

InteractiveListProperty::~InteractiveListProperty()
{
}

void InteractiveListProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    s.serialize("items", items_); // Serialize to check for changes in item list on next deserialization.
    s.serialize("instances", instances_);

}
void InteractiveListProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    std::vector<std::string> oldItems;
    s.deserialize("items", oldItems);

    if(serializeItems_) {
        items_ = oldItems;
    }

    // Reordering items invalidates instances.
    // Therefore, we need to remap old ids to their current equivalent.

    // Figure out mapping between old and new item ids.
    std::map<int, int> itemIdMappingTable;
    for(size_t oldItemIdx=0; oldItemIdx<oldItems.size(); oldItemIdx++) {
        for(size_t newItemIdx=0; newItemIdx<items_.size(); newItemIdx++) {
            if(oldItems[oldItemIdx] == items_[newItemIdx]) {
                itemIdMappingTable[oldItemIdx] = newItemIdx;
                break;
            }
        }
    }

    std::vector<Instance> oldInstances;
    try {
        s.deserialize("instances", oldInstances);
    }
    catch(SerializationException&) {
        s.removeLastError();
        LERROR("You need to reconfigure InteractiveList: " << getFullyQualifiedGuiName());
    }

    // Only add instances, whose item still exists.
    for(const Instance& oldInstance : oldInstances) {
        auto it = itemIdMappingTable.find(oldInstance.getItemId());
        if(it != itemIdMappingTable.end() && (allowDuplication_ || !hasInstance(items_[it->second]))) {
            // Instance id remains.
            Instance instance(it->second, oldInstance.getInstanceId());
            instance.setActive(oldInstance.isActive());
            instance.setName(oldInstance.getName());
            instances_.push_back(instance);
        }
    }
}

size_t InteractiveListProperty::getNumItems() const {
    return items_.size();
}
void InteractiveListProperty::reset() {
    setItems(items_);
}

void InteractiveListProperty::clear() {
    items_.clear();
    instances_.clear();
    selectedInstance_ = -1;
    invalidate();
}

void InteractiveListProperty::setItems(const std::vector<std::string>& items) {
    instances_.clear();
    items_ = items;

    selectedInstance_ = -1;

    invalidate();
}

void InteractiveListProperty::addItem(const std::string& item) {
    tgtAssert(std::find(items_.begin(), items_.end(), item) == items_.end(), "Item already added");
    items_.push_back(item);

    invalidate();
}

void InteractiveListProperty::removeItem(const std::string& item) {
    auto iter = std::find(items_.begin(), items_.end(), item);
    if(iter == items_.end())
        return;

    // Reset input and output list, since indices will change due to removal.
    reset();

    items_.erase(iter);

    invalidate();
}

const std::vector<std::string>& InteractiveListProperty::getItems() const {
    return items_;
}

std::vector<int> InteractiveListProperty::getInputIndices() const {
    if(allowDuplication_) {
        std::vector<int> inputIndices(items_.size());
        std::iota(inputIndices.begin(), inputIndices.end(), 0);
        return inputIndices;
    }
    else {
        std::vector<int> inputIndices;
        for (size_t i=0; i<items_.size(); i++) {
            if (!hasInstance(items_[i])) {
                inputIndices.push_back(i);
            }
        }
        return inputIndices;
    }
}

const std::vector<InteractiveListProperty::Instance>& InteractiveListProperty::getInstances() const {
    return instances_;
}

std::vector<InteractiveListProperty::Instance>& InteractiveListProperty::getInstances() {
    return instances_;
}

std::vector<InteractiveListProperty::Instance> InteractiveListProperty::getInstances(const std::string& item) const {
    std::vector<Instance> instances;
    for(const Instance& instance : instances_) {
        if(items_[instance.getItemId()] == item) {
            instances.push_back(instance);
        }
    }
    return instances;
}

void InteractiveListProperty::addInstance(const std::string& item, int pos) {

    if(!allowDuplication_ && hasInstance(item))
        return;

    int index = getIndexOfItem(item);
    tgtAssert(index >= 0, "Item not contained");

    Instance instance = createInstance(index);
    if(pos < 0) {
        instances_.push_back(instance);
    }
    else {
        pos = std::min(pos, static_cast<int>(instances_.size()));
        instances_.insert(instances_.begin() + pos, instance);

        if(selectedInstance_ > -1 && pos <= selectedInstance_ ) {
            selectedInstance_++;
        }
    }

    invalidate();
}

void InteractiveListProperty::removeInstance(int instanceId) {

    int idx = getIndexOfInstance(instanceId);
    if(idx == -1)
        return;

    auto instance = instances_.begin() + idx;
    instances_.erase(instance);

    if(selectedInstance_ > -1
       && (selectedInstance_ == static_cast<int>(instances_.size())
           || idx > selectedInstance_)) {
        selectedInstance_--;
    }

    invalidate();
}

void InteractiveListProperty::moveInstance(int instanceId, int pos) {
    tgtAssert(pos >= 0 && pos <= static_cast<int>(instances_.size()), "Position out of range");

    int idx = getIndexOfInstance(instanceId);
    tgtAssert(idx != -1, "Instance not available");
    if(idx == -1 || pos == idx) {
        return;
    }

    // First erase the instance at its old position.
    Instance instance = instances_[idx];
    instances_.erase(instances_.begin() + idx);

    // Insert the instance at its new position.
    if(idx < pos) {
        pos--;
    }
    instances_.insert(instances_.begin() + pos, instance);

    // Update selection.
    if(selectedInstance_ == idx) {
        selectedInstance_ = pos;
    }
    else if(selectedInstance_ == pos) {
        selectedInstance_ = idx;
    }

    invalidate();
}

void InteractiveListProperty::swapInstances(int instanceId, int pos) {
    tgtAssert(pos >= 0 && pos < static_cast<int>(instances_.size()), "Position out of range");

    int idx = getIndexOfInstance(instanceId);
    tgtAssert(idx != -1, "Instance not available");
    if(idx == -1 || pos == idx)
        return;

    // Swap positions.
    std::iter_swap(instances_.begin() + idx, instances_.begin() + pos);

    // Update selection.
    if(selectedInstance_ == idx) {
        selectedInstance_ = pos;
    }
    else if(selectedInstance_ == pos) {
        selectedInstance_ = idx;
    }

    invalidate();
}

bool InteractiveListProperty::hasInstance(const std::string& item) const {
    return !getInstanceIds(item).empty();
}

bool InteractiveListProperty::isDuplicationAllowed() const {
    return allowDuplication_;
}

void InteractiveListProperty::setDuplicationAllowed(bool enabled) {
    if(enabled != allowDuplication_) {
        allowDuplication_ = enabled;
        reset();
    }
}

int InteractiveListProperty::getSelectedInstance() const {
    return selectedInstance_;
}

void InteractiveListProperty::setSelectedInstance(int index) {
    if(selectedInstance_ != index) {
        tgtAssert(selectedInstance_ >= -1 && selectedInstance_ < static_cast<int>(instances_.size()), "Invalid instance index");
        selectedInstance_ = index;
        invalidate();
    }
}

int InteractiveListProperty::getIndexOfItem(const std::string& item) const {
    for(int i = 0; items_.begin() + i != items_.end(); i++) {
        if(*(items_.begin() + i) == item) {
            return i;
        }
    }
    return -1;
}

int InteractiveListProperty::getIndexOfInstance(int instanceId) const {
    for(size_t i = 0; i < instances_.size(); i++) {
        if(instances_[i].getInstanceId() == instanceId) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

std::vector<int> InteractiveListProperty::getInstanceIds(const std::string& name) const {
    std::vector<int> ids;
    for(size_t i = 0; i < instances_.size(); i++) {
        if(items_[instances_[i].getItemId()] == name) {
            ids.push_back(static_cast<int>(i));
        }
    }
    return ids;
}

InteractiveListProperty::Instance InteractiveListProperty::createInstance(int itemId) const {

    //TODO: useful IDs? / Reuse removed instance IDs?
    int instanceId = 0;
    for(const Instance& other : instances_) {
        instanceId = std::max(other.getInstanceId(), instanceId);
    }
    instanceId++;

    Instance instance(itemId, instanceId);
    instance.setName(nameGenerator_(instance));

    return instance;
}

void InteractiveListProperty::setNameGenerator(const NameGenerator& nameGenerator) {
    nameGenerator_ = nameGenerator;
}

const InteractiveListProperty::NameGenerator& InteractiveListProperty::getNameGenerator() const {
    return nameGenerator_;
}

void InteractiveListProperty::setInstanceLabel(const std::string& label) {
    instanceLabel_ = label;
}

const std::string& InteractiveListProperty::getInstanceLabel() const {
    return instanceLabel_;
}

void InteractiveListProperty::setItemLabel(const std::string& label) {
    itemLabel_ = label;
}

const std::string& InteractiveListProperty::getItemLabel() const {
    return itemLabel_;
}

} // namespace voreen
