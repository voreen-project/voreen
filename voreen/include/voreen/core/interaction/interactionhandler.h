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

#ifndef VRN_INTERACTIONHANDLER_H
#define VRN_INTERACTIONHANDLER_H

#include "tgt/event/eventlistener.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/eventproperty.h"

namespace voreen {

class Processor;

/**
 * \brief Abstract base class for interaction handlers that can be attached to processors
 * and react to user events. Implement onEvent() in a subclass to receive the events.
 *
 * Each InteractionHandler instance is associated with one or multiple event properties
 * that are used for event filtering: Processor base class passes incoming events
 * to its interaction handlers' event properties. For concrete handlers that only
 * need to react to a single mouse or keyboard event, it usually sufficient to call
 * the appropriate constructor of the InteractionHandler base class, which creates
 * the needed EventProperty. More complex handlers (e.g. CameraInteractionHandler),
 * however, may need to create and add event properties on their own.
 *
 * Since InteractionHandler derives from PropertyOwner, a derived class
 * may also use "normal" (non-event) properties.
 */
class VRN_CORE_API InteractionHandler : public PropertyOwner, public tgt::EventListener {

    friend class Processor;
    friend class EventPropertyBase;

public:

    /**
     * Dummy constructor. Since no event property is created, the handler is inactive.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *     a of processor. Must not be empty.
     * @param guiName the string that is to be displayed in the GUI
     */
    InteractionHandler(const std::string& id, const std::string& guiName);

    /**
     * Constructs an mouse interaction handler by internally creating an MouseEventProperty
     * with the passed parameters.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *     a of processor. Must not be empty.
     * @param guiName the string that is to be displayed in the GUI
     * @param mouseButtons Mouse buttons the handler is supposed to react to
     * @param mouseAction Mouse action the handler is supposed to react to
     * @param modifier Modifier that is required to be active for having the interaction handler react to the event
     * @param shareEvents Determines whether the handler shares events, i.e., does not accept events it reacts to
     * @param enabled If set to false, the handler is completely transparent, not reacting to any event
     */
    InteractionHandler(const std::string& id, const std::string& guiName,
        tgt::MouseEvent::MouseButtons mouseButtons, tgt::MouseEvent::MouseAction mouseAction,
        tgt::Event::Modifier modifier, bool shareEvents, bool enabled);

    /**
     * Constructs a keyboard interaction handler by internally creating an KeyboardEventProperty
     * with the passed parameters.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *     a of processor. Must not be empty.
     * @param guiName the string that is to be displayed in the GUI
     * @param keyCode The key code the handler is supposed to react to
     * @param modifier Modifier that is required to be active for having the interaction handler accept the event
     * @param shareEvents Determines whether the handler shares events, i.e., does not accept events it reacts to
     * @param enabled If set to false, the handler is completely transparent, not reacting to any event
     */
    InteractionHandler(const std::string& id, const std::string& guiName,
        tgt::KeyEvent::KeyCode keyCode, tgt::Event::Modifier modifier,
        bool shareEvents, bool enabled);

    /**
     * Constructs an interaction handler that reacts to both mouse and key events.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *     a of processor. Must not be empty.
     * @param guiName the string that is to be displayed in the GUI
     * @param mouseButtons Mouse buttons the handler is supposed to react to
     * @param mouseAction Mouse action the handler is supposed to react to
     * @param keyCode The key code the handler is supposed to react to
     * @param modifier Modifier that is required to be active for having the interaction handler accept the event
     * @param shareEvents Determines whether the handler shares events, i.e., does not accept events it reacts to
     * @param enabled If set to false, the handler is completely transparent, not reacting to any event
     */
     InteractionHandler(const std::string& id, const std::string& guiName,
        tgt::MouseEvent::MouseButtons mouseButtons, tgt::MouseEvent::MouseAction mouseAction,
        tgt::KeyEvent::KeyCode keyCode, tgt::Event::Modifier modifier, bool shareEvents, bool enabled);

    /**
     * Deletes the associated event properties.
     */
    virtual ~InteractionHandler();

    /**
     * Returns the processor this handler is assigned to.
     */
    Processor* getOwner() const;

    /**
     * Returns whether the handler accepts events. Internally calls
     * isEnabled() of its event property.
     */
    bool isEnabled() const;

    /**
     * Returns whether the handler shares events, i.e., does not accept() events it reacts to.
     * Internally calls sharesEvents() on its event property.
     */
    bool isSharing() const;

    /**
     * Returns wether the handler would react to the passed event.
     * Internally calls accepts() on its event property.
     */
    bool accepts(tgt::Event* e) const;

    /**
     * Determines whether the handler accepts events.
     * Internally calls setEnabled() on its event property.
     */
    void setEnabled(bool enabled);

    /**
     * Determines whether its properties are visible.
     * Internally calls setVisible() on its event properties
     * and PropertyOwner::setPropertiesVisible.
     */
    void setVisible(bool state);

    /**
     * Determines whether the handler shares events, i.e., does not accept() events it reacts to.
     * Internally calls setEnabled() on its event property.
     */
    void setSharing(bool sharing);

    /**
     * Returns the associated event properties.
     */
    const std::vector<EventPropertyBase*>& getEventProperties() const;

protected:
    /**
     * To be implemented by the concrete InteractionHandler subclass.
     * Is called via the associated event property, if the event matches.
     */
    virtual void onEvent(tgt::Event* e) = 0;

    /**
     * Assigns the handler's owning processor.
     */
    void setOwner(Processor* p);

    /**
     * Adds an event property, used for passing events to this handler.
     *
     * @note The base class takes ownership of the passed object.
     *      Therefore, the caller must delete it on its own.
     */
    void addEventProperty(EventPropertyBase* eventProperty);

    Processor* owner_;              ///< The processor this interaction handler is associated with

    std::vector<EventPropertyBase*> eventProperties_;  ///< Event filter properties, used by owning processor for passing events
};

} // namespace

#endif // VRN_INTERACTIONHANDLER_H
