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

#include "interactivelistproperty.h"

namespace voreen {

InteractiveListProperty::InteractiveListProperty(const std::string& id, const std::string& guiText, bool allowDuplication,
                            int invalidationLevel, Property::LevelOfDetail lod)
    : Property(id, guiText, invalidationLevel, lod)
    , allowDuplication_(allowDuplication)
    , selectedInstance_(-1)
{
}

InteractiveListProperty::InteractiveListProperty()
{
}

InteractiveListProperty::~InteractiveListProperty()
{
}

void InteractiveListProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    s.serialize("items", items_);
    s.serialize("inputItemIds", inputItemIds_);
    s.serializeBinaryBlob("instances", instances_);

}
void InteractiveListProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);
    s.deserialize("items", items_);
    s.deserialize("inputItemIds", inputItemIds_);
    s.deserializeBinaryBlob("instances", instances_);
}

size_t InteractiveListProperty::getNumItems() const {
    return items_.size();
}
void InteractiveListProperty::reset() {
    setItems(items_);
}

void InteractiveListProperty::clear() {
    items_.clear();
    inputItemIds_.clear();
    instances_.clear();
    selectedInstance_ = -1;
    invalidate();
}

void InteractiveListProperty::setItems(const std::vector<std::string>& items) {
    inputItemIds_.clear();
    instances_.clear();
    items_ = items;

    for(size_t i=0; i<getNumItems(); i++) {
        inputItemIds_.push_back(static_cast<int>(i));
    }

    selectedInstance_ = -1;

    invalidate();
}

void InteractiveListProperty::addItem(const std::string& item) {
    tgtAssert(std::find(items_.begin(), items_.end(), item) == items_.end(), "Item already added");
    items_.push_back(item);
    inputItemIds_.push_back(static_cast<int>(items_.size()-1));

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

const std::vector<int>& InteractiveListProperty::getInputIndices() const {
    return inputItemIds_;
}

const std::vector<InteractiveListProperty::Instance>& InteractiveListProperty::getInstances() const {
    return instances_;
}

std::vector<InteractiveListProperty::Instance> InteractiveListProperty::getInstances(const std::string& item) const {
    std::vector<Instance> instances;
    for(const Instance& instance : instances_) {
        if(items_[instance.itemId_] == item) {
            instances.push_back(instance);
        }
    }
    return instances;
}

void InteractiveListProperty::addInstance(const std::string& item, int pos) {
    int index = getIndexOfItem(item);
    tgtAssert(index >= 0, "Item not contained");

    auto iter = std::find(inputItemIds_.begin(), inputItemIds_.end(), index);
    if(iter == inputItemIds_.end())
        return;

    if(!allowDuplication_) {
        inputItemIds_.erase(iter);
    }

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
    int itemId = instance->itemId_;
    instances_.erase(instance);

    if(!allowDuplication_) {
        size_t pos = 0;
        while(pos < inputItemIds_.size() && inputItemIds_[pos] < itemId) pos++;
        inputItemIds_.insert(inputItemIds_.begin()+pos, itemId);
    }

    if(selectedInstance_ > -1
       && (selectedInstance_ == static_cast<int>(instances_.size())
           || idx > selectedInstance_)) {
        selectedInstance_--;
    }

    invalidate();
}

void InteractiveListProperty::moveInstance(int instanceId, int pos) {
    tgtAssert(pos >= 0 && pos < static_cast<int>(instances_.size()), "Position out of range");

    int idx = getIndexOfInstance(instanceId);
    tgtAssert(idx != -1, "Instance not available");
    if(idx == -1)
        return;

    // Swap positions.
    std::iter_swap(instances_.begin() + idx, instances_.begin() + pos);

    // Update selection.
    if(selectedInstance_ == idx) {
        selectedInstance_ = pos;
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
        if(instances_[i].instanceId_ == instanceId) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

std::vector<int> InteractiveListProperty::getInstanceIds(const std::string& name) const {
    std::vector<int> ids;
    for(size_t i = 0; i < instances_.size(); i++) {
        if(items_[instances_[i].itemId_] == name) {
            ids.push_back(static_cast<int>(i));
        }
    }
    return ids;
}

InteractiveListProperty::Instance InteractiveListProperty::createInstance(int itemId) const {
    Instance instance;
    instance.itemId_ = itemId;
    instance.instanceId_ = 0;

    //TODO: useful IDs? / Reuse removed instance IDs?
    for(const Instance& other : instances_) {
        instance.instanceId_ = std::max(other.instanceId_, instance.instanceId_);
    }
    instance.instanceId_++;

    return instance;
}

std::string InteractiveListProperty::getInstanceName(const voreen::InteractiveListProperty::Instance& instance) const {
    std::string name = items_[instance.itemId_];

    if(allowDuplication_) {
        name += " #" + std::to_string(instance.instanceId_);
    }

    return name;
}

} // namespace voreen
