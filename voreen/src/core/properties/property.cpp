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

#include "voreen/core/properties/property.h"
#include "voreen/core/properties/propertywidget.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"

#include "tgt/glcontextmanager.h"

namespace voreen {

Property::Property(const std::string& id, const std::string& guiText, int invalidationLevel, LevelOfDetail lod)
    : VoreenSerializableObject(id, guiText)
    , owner_(0)
    , invalidationLevel_(invalidationLevel)
    , groupId_("")
    , groupName_("")
    /// view related enums
    , levelOfDetail_(lod)
    , viewFlags_(VF_VISIBLE)

    , interactionModeVisited_(false)
    , serializeValue_(false)
    , linkCheckVisited_(false)
    , initialGuiName_(guiText)
    , isInitialized_(false)
{
//    tgtAssert(!id.empty(), "Property's id must not be empty");
}

Property::Property()
    : VoreenSerializableObject()
    , owner_(0)
    , invalidationLevel_(Processor::INVALID_RESULT)
    , groupId_("")
    , groupName_("")
    /// view related enums
    , levelOfDetail_(LOD_DEFAULT)
    , viewFlags_(VF_VISIBLE)

    , interactionModeVisited_(false)
    , serializeValue_(false)
    , linkCheckVisited_(false)
    , initialGuiName_("")
    , isInitialized_(false)
{}

Property::~Property() {
    tgtAssert(!isInitialized_, "Property has not been deinitialized before destruction");
    disconnectWidgets();
    std::set<PropertyWidget*>::iterator it = widgets_.begin();
    for ( ; it != widgets_.end(); ++it) {
        delete (*it);
    }
}

int Property::getInvalidationLevel() const {
    return invalidationLevel_;
}

void Property::setInvalidationLevel(int invalidationLevel) {
    invalidationLevel_ = invalidationLevel;
}

void Property::initialize() {
    tgtAssert(!isInitialized_, "Property is already initialized");
    isInitialized_ = true;
    // currently nothing else to do
}

void Property::deinitialize() {
    tgtAssert(isInitialized_, "Property has not been initialized");
    isInitialized_ = false;
    // currently nothing else to do
}

bool Property::isInitialized() const {
    return isInitialized_;
}

void Property::setReadOnlyFlag(bool state) {
    viewFlags_ = viewFlags_ & ~VF_READ_ONLY;
    viewFlags_ = (state ? viewFlags_ | VF_READ_ONLY : viewFlags_);
    for (std::set<PropertyWidget*>::iterator it = widgets_.begin(); it != widgets_.end(); ++it)
        (*it)->updateViewFlags(viewFlags_);
}

bool Property::isReadOnlyFlagSet() const {
    return viewFlags_ & VF_READ_ONLY;
}

void Property::setVisibleFlag(bool state) {
    viewFlags_ = viewFlags_ & ~VF_VISIBLE;
    viewFlags_ = (state ? viewFlags_ | VF_VISIBLE : viewFlags_);

    // adjust visibility of assigned gui widgets
    std::set<PropertyWidget*>::iterator it = widgets_.begin();
    for ( ; it != widgets_.end(); ++it)
        (*it)->updateViewFlags(viewFlags_);
}

bool Property::isVisibleFlagSet() const {
    return viewFlags_ & VF_VISIBLE;
}

void Property::setViewFlags(ViewFlags flags) {
    viewFlags_ = viewFlags_ & ~VF_VISIBLE & ~VF_READ_ONLY;
    viewFlags_ = viewFlags_ | flags;

    // adjust visibility of assigned gui widgets
    std::set<PropertyWidget*>::iterator it = widgets_.begin();
    for ( ; it != widgets_.end(); ++it)
        (*it)->updateViewFlags(viewFlags_);}

Property::ViewFlags Property::getViewFlags() const {
    return viewFlags_;
}

Property::LevelOfDetail Property::getLevelOfDetail() const {
    return levelOfDetail_;
}

void Property::setGroupID(const std::string& gid) {
    groupId_ = gid;
}

std::string Property::getGroupID() const {
    return groupId_;
}

void Property::setGroupName(const std::string& name) {
    groupName_ = name;
}

std::string Property::getGroupName() const {
    return groupName_;
}

void Property::setOwner(PropertyOwner* processor) {
    owner_ = processor;
}

PropertyOwner* Property::getOwner() const {
    return owner_;
}

void Property::setGuiName(const std::string& guiName) {
    VoreenObject::setGuiName(guiName);
    for(std::set<PropertyWidget*>::iterator it = widgets_.begin(); it != widgets_.end(); it++)
        (*it)->setPropertyGuiName(guiName);
}

void Property::addWidget(PropertyWidget* widget) {
    if (widget) {
        widget->updateViewFlags(viewFlags_);
        widgets_.insert(widget);
    }
}

void Property::removeWidget(PropertyWidget* widget) {
    if (widget)
        widgets_.erase(widget);
}

void Property::disconnectWidgets() {
    std::set<PropertyWidget*>::iterator it = widgets_.begin();
    for ( ; it != widgets_.end(); ++it)
        (*it)->disconnect();
}

void Property::updateWidgets() {
    std::set<PropertyWidget*>::iterator it = widgets_.begin();
    for ( ; it != widgets_.end(); ++it)
        (*it)->updateFromProperty();
}

std::string Property::getFullyQualifiedID() const {
    if (getOwner())
        return getOwner()->getID() + "." + getID();
    else
        return getID();
}

std::string Property::getFullyQualifiedGuiName() const {
    if (getOwner())
        return getOwner()->getGuiName() + "." + getGuiName();
    else
        return getGuiName();
}

void Property::serialize(Serializer& s) const {
    if(serializeValue_)
        return;

    s.serialize("name", id_);

    if (guiName_ != initialGuiName_)
        s.serialize("guiName", guiName_);

    for (std::set<PropertyWidget*>::const_iterator it = widgets_.begin(); it != widgets_.end(); ++it) {
        (*it)->updateMetaData();
        // FIXME What exactly is this supposed to do? The return value is not used... FL
        (*it)->getWidgetMetaData();
    }

    metaDataContainer_.serialize(s);
}

void Property::deserialize(Deserializer& s) {
    if (serializeValue_)
        return;

    // deserialize gui name, if available
    try {
        std::string temp;
        s.deserialize("guiName", temp);
        guiName_ = temp;
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
    }

    metaDataContainer_.deserialize(s);
}

void Property::serializeValue(Serializer& s) {
    serializeValue_ = true;
    serialize(s);
    serializeValue_ = false;
}

void Property::deserializeValue(Deserializer& s) {
    serializeValue_ = true;
    deserialize(s);
    serializeValue_ = false;
}

MetaDataContainer& Property::getMetaDataContainer() const {
    return metaDataContainer_;
}

void Property::registerLink(PropertyLink* link) {
    links_.push_back(link);
}

void Property::removeLink(PropertyLink* link) {
    for (std::vector<PropertyLink*>::iterator it = links_.begin(); it != links_.end(); ++it) {
        if (*it == link) {
            links_.erase(it);
            break;
        }
    }
}

const std::set<PropertyWidget*> Property::getPropertyWidgets() const {
    return widgets_;
}

void Property::invalidateOwner() {
    invalidateOwner(invalidationLevel_);
}

void Property::invalidateOwner(int invalidationLevel) {
    //read only properties do not invalidate their owner
    //if(!isReadOnlyFlagSet()) { has unexpected beahviour, if blocked
        if (PropertyOwner* po = getOwner()) {
            po->notifyPropertyValueHasBeenModified(this);
            po->invalidate(invalidationLevel);
        }
    //}
}

void Property::invalidate() {
    // first execute local stuff
    onChangeCallbacks_.execute();
    // notify widgets of updated values
    updateWidgets();
    // invalidate the owner
    invalidateOwner();
    // tell others about the change
    executeLinks();
}



const std::vector<PropertyLink*>& Property::getLinks() const {
    return links_;
}

PropertyLink* Property::getLink(const Property* dest) const {
    for (size_t i=0; i<links_.size(); ++i) {
        if (links_[i]->getDestinationProperty() == dest)
            return links_[i];
    }
    return 0;
}

bool Property::isLinkedWith(const Property* dest, bool transitive) const {
    // true if dest is the object itself
    if (this == dest)
        return true;

    // check for direct links
    bool linkedDirectly = (getLink(dest) != 0);
    if (linkedDirectly)
        return true;

    // recursive search for indirect links
    if (transitive) {
        if (linkCheckVisited_)
            return false;
        else {
            linkCheckVisited_ = true;
            for (size_t i=0; i<links_.size(); ++i) {
                if (links_[i]->getDestinationProperty()->isLinkedWith(dest, true)) {
                    linkCheckVisited_ = false;
                    return true;
                }
            }
            linkCheckVisited_ = false;
        }
    }

    return false;
}

bool Property::isLinkableWith(const voreen::Property* dst) const{
    return (LinkEvaluatorHelper::arePropertiesLinkable(this, dst));
}

std::vector<std::pair<std::string, std::string> > Property::getCompatibleEvaluators(const voreen::Property* dst) const{
    return LinkEvaluatorHelper::getCompatibleLinkEvaluators(this, dst);
}

std::string Property::getTypeDescription() const {
    return "<unknown>";
}

std::string Property::getDescription() const {
    return description_;
}

void Property::setDescription(std::string desc) {
    description_ = desc;
}

void Property::onChange(const Callback& action) {
    onChangeCallbacks_.registerCallback(action);
}

void Property::blockCallbacks(bool block) {
    onChangeCallbacks_.blockCallbacks(block);
}

void Property::executeLinks() {
    if (links_.empty())
        return;

    // check, if this property is the initiator
    bool isInitiator = PropertyLink::visitedProperties_.empty();
    if (isInitiator) {
        PropertyLink::visitedProperties_.push_back(this);
    }

    // pass change data object to links
    for (std::vector<PropertyLink*>::iterator it = links_.begin(); it != links_.end(); it++) {
        try {
            (*it)->executeLink();
        }
        catch (const VoreenException& e) {
            LERRORC("voreen.Property", "executeLinks(): " << e.what());
        }
    }

    //clean-up vector
    if (isInitiator) {
        PropertyLink::visitedProperties_.clear();
    }
}


} // namespace voreen
