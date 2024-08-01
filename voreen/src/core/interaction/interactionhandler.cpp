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

#include "voreen/core/interaction/interactionhandler.h"

using tgt::Event;

namespace voreen {

InteractionHandler::InteractionHandler(const std::string& id, const std::string& guiName,
        tgt::MouseEvent::MouseButtons mouseButtons,
        tgt::MouseEvent::MouseAction mouseAction,
        tgt::Event::Modifier modifier,
        bool shareEvents, bool enabled)
    : tgt::EventListener()
    , PropertyOwner(id, guiName)
    , owner_(0)
{
    tgtAssert(!id.empty(), "InteractionHandler id must not be empty");

    EventPropertyBase* eventProp = new EventProperty<InteractionHandler>(id, guiName, this, &InteractionHandler::onEvent,
        mouseButtons, mouseAction, tgt::KeyEvent::K_UNKNOWN, modifier, shareEvents, enabled);
    eventProp->setReceivesKeyEvents(false);
    addEventProperty(eventProp);
}

InteractionHandler::InteractionHandler(const std::string& id, const std::string& guiName,
        tgt::KeyEvent::KeyCode keycode, tgt::Event::Modifier modifier,
        bool shareEvents, bool enabled) :
    tgt::EventListener(),
    owner_(0)
{
    id_ = id;
    guiName_ = guiName;
    tgtAssert(!id.empty(), "InteractionHandler id must not be empty");

    EventPropertyBase* eventProp = new EventProperty<InteractionHandler>(id, guiName, this, &InteractionHandler::onEvent,
        tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::ACTION_NONE,
        keycode, modifier, shareEvents, enabled);
    eventProp->setReceivesMouseEvents(false);
    addEventProperty(eventProp);
}

InteractionHandler::InteractionHandler(const std::string& id, const std::string& guiName,
        tgt::MouseEvent::MouseButtons mouseButtons,
        tgt::MouseEvent::MouseAction mouseAction,
        tgt::KeyEvent::KeyCode keycode,
        tgt::Event::Modifier modifier,
        bool shareEvents, bool enabled) :
    tgt::EventListener(),
    owner_(0)
{
    id_ = id;
    guiName_ = guiName;
    tgtAssert(!id.empty(), "InteractionHandler id must not be empty");

    addEventProperty(new EventProperty<InteractionHandler>(id, guiName, this, &InteractionHandler::onEvent,
       mouseButtons, mouseAction, keycode, modifier, shareEvents, enabled));
}

InteractionHandler::InteractionHandler(const std::string& id, const std::string& guiName)
    : PropertyOwner(id, guiName)
    , tgt::EventListener()
    , owner_(0)
{
    tgtAssert(!id.empty(), "InteractionHandler id must not be empty");
}

InteractionHandler::~InteractionHandler() {
    for (size_t i=0; i<eventProperties_.size(); i++)
        delete eventProperties_[i];
    eventProperties_.clear();
}

Processor* InteractionHandler::getOwner() const {
    return owner_;
}

void InteractionHandler::setOwner(Processor* p) {
    if (owner_)
        LWARNINGC("voreen.InteractionHandler", "Setting new owner for InteractionHandler!");
    owner_ = p;
}

const std::vector<EventPropertyBase*>& InteractionHandler::getEventProperties() const {
    return eventProperties_;
}

void InteractionHandler::setEnabled(bool enabled) {
    if (eventProperties_.empty()) {
        LWARNINGC("InteractionHandler", "No event property assigned");
    }
    for (size_t i=0; i<eventProperties_.size(); i++)
        eventProperties_[i]->setEnabled(enabled);
}

void InteractionHandler::setVisible(bool state) {
    if (eventProperties_.empty()) {
        LWARNINGC("InteractionHandler", "No event property assigned");
    }
    for (size_t i=0; i<eventProperties_.size(); i++)
        eventProperties_[i]->setVisibleFlag(state);
}

void InteractionHandler::setSharing(bool sharing) {
    if (eventProperties_.empty()) {
        LWARNINGC("InteractionHandler", "No event property assigned");
    }
    for (size_t i=0; i<eventProperties_.size(); i++)
        eventProperties_[i]->setSharing(sharing);
}

bool InteractionHandler::accepts(tgt::Event* e) const {
    for (size_t i=0; i<eventProperties_.size(); i++)
        if (eventProperties_[i]->accepts(e))
            return true;

    return false;
}

void InteractionHandler::addEventProperty(EventPropertyBase* eventProperty) {
    eventProperties_.push_back(eventProperty);
    addProperty(eventProperty);
}

} // namespace
