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

#ifndef VRN_EVENTPROPERTY_H
#define VRN_EVENTPROPERTY_H

#include "voreen/core/properties/property.h"
#include "tgt/event/event.h"
#include "tgt/event/keyevent.h"
#include "tgt/event/mouseevent.h"
#include "tgt/event/touchevent.h"

namespace voreen {

/**
 * Abstract, non-template base class of EventProperty.
 *
 * @see EventProperty
 */
class VRN_CORE_API EventPropertyBase : public Property {

    friend class Processor;

public:

    /**
     * Constructor.
     *
     * @param id property identifier. Must be unique among the
     *  properties of a PropertyOwner and must not be empty.
     * @param guiName string to be displayed in a user interface
     * @param mouseEvents determines whether the property receives mouse events
     * @param keyEvents determines whether the property receives key events
     * @param mouseButtons mouse buttons the property reacts to
     * @param mouseAction mouse actions the property reacts to
     * @param keyCode code of the key the property reacts to
     * @param modifier modifier that has to pressed for accepting the event
     * @param shareEvents if set to true, the property does ignore the event,
     *  even if it matches and is passed to its owner.
     * @param enabled if set to false, the property does not react to any event
     */
    EventPropertyBase(const std::string& id, const std::string& guiName,
        bool mouseEvents, bool keyEvents,
        tgt::MouseEvent::MouseButtons mouseButtons,
        int mouseAction,
        tgt::KeyEvent::KeyCode keyCode,
        tgt::Event::Modifier modifier,
        bool shareEvents, bool enabled);

    virtual std::string getClassName() const       { return "EventProperty"; }
    virtual std::string getTypeDescription() const { return "EventProperty"; }
    virtual void reset(){}
    /**
     * Returns true, if the event property accepts the event \p e.
     *
     * A property might react to both mouse and key events, or to only one event type:
     *  - a mouse event is accepted, if
     *     -# the event property is enabled, @see isEnabled
     *     -# the property is configured to receive mouse events, @see receivesMouseEvents
     *     -# the event's mouse button, mouse action and modifier matches the property's configuration
     *  - a keyboard event is accepted, if
     *     -# the event property is enabled, @see isEnabled
     *     -# the property is configured to receive key events, @see receivesKeyEvents
     *     -# the event's key code and modifier matches the property's configuration
     */
    virtual bool accepts(tgt::Event* e) const;

    /**
     * Determines whether the property listens to events.
     *
     * @param enabled if false, the property ignores all incoming events.
     */
    void setEnabled(bool enabled);

    /**
     * Returns whether the property listens to events.
     */
    bool isEnabled() const;

    /**
     * Determines whether the event property listens to mouse events.
     */
    void setReceivesMouseEvents(bool mouseEvents);

    /**
     * Returns whether the event property listens to mouse events.
     */
    bool receivesMouseEvents() const;

    /**
     * Determines whether the event property listens to key events.
     */
    void setReceivesKeyEvents(bool keyEvents);

    /**
     * Returns whether the event property listens to key events.
     */
    bool receivesKeyEvents() const;

    /**
     * Determines whether the event property listens to touch events.
     */
    void setReceivesTouchEvents(bool touchEvents);

    /**
     * Returns whether the event property listens to touch events.
     */
    bool receivesTouchEvents() const;

    /**
     * Determines the mouse buttons the property reacts to.
     *
     * @see accepts
     * @see enabled
     * @see setReceivesMouseEvents
     */
    void setMouseButtons(tgt::MouseEvent::MouseButtons buttons);

    /**
     * Returns the mouse buttons the property reacts to.
     *
     * @see setMouseButtons
     */
    tgt::MouseEvent::MouseButtons getMouseButtons() const;

    /**
     * Determines the mouse actions (press, release, move, ...) the property reacts to.
     *
     * @see accepts
     * @see enabled
     * @see setReceivesMouseEvents
     */
    void setMouseAction(int action);

    /**
     * Determines the mouse actions the property reacts to.
     *
     * @see setMouseAction
     */
    int getMouseAction() const;

    /**
     * Determines the code of the key the property reacts to.
     *
     * @see accepts
     * @see enabled
     * @see setReceivesKeyEvents
     */
    void setKeyCode(tgt::KeyEvent::KeyCode keyCode);

    /**
     * Returns the code of the key the property reacts to.
     *
     * @see setKeyCode
     */
    tgt::KeyEvent::KeyCode getKeyCode() const;

    /**
     * Determines the modifier that has to be pressed
     * for allowing the property to receive a mouse or key event.
     *
     * @see accepts
     */
    void setModifier(tgt::Event::Modifier value);

    /**
     * Returns the modifier of an event to be accepted.
     *
     * @see setModifier
     */
    tgt::Event::Modifier getModifier() const;

    /**
     * Determines whether the property allows to further
     * propagate events that are accepted by the property.
     *
     * If sharing is enabled, the property does call ignore()
     * on the event after it has been passed to the owner.
     * Otherwise, the decision to accept an event
     * is completely up to the owner, i.e. the owner has to
     * explicitly call accept() on the event, in order to
     * prevent further propagation.
     */
    void setSharing(bool shareEvents);

    /**
     * Returns whether the property shares events.
     *
     * @see setSharing
     */
    bool isSharing() const;

    /// @see Serializable::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Serializable::deserialize
    virtual void deserialize(Deserializer& s);

protected:

    /**
     * Used by Processor::onEvent of the owning processor
     * in order to relay received events.
     */
    virtual void execute(tgt::Event* e) = 0;

    /// Implemented by generic EventProperty.
    virtual void notifyChangeListener() = 0;

    bool receivesMouseEvents_;  ///< Determines whether mouse event receiving is enabled for this processor.
    bool receivesKeyEvents_;    ///< Determines whether key event receiving is enabled for this processor.
    bool receivesTouchEvents_;  ///< Determines whether touch event receiving is enabled for this processor.
    bool enabled_;              ///< If false, the property does not receive any events.

    tgt::MouseEvent::MouseButtons mouseButtons_; ///< Buttons the property reacts to
    int mouseAction_;                            ///< Actions the property reacts to
    tgt::KeyEvent::KeyCode keyCode_;             ///< Code of the key the property reacts to
    tgt::Event::Modifier modifier_;              ///< Modifier an event needs to have sets for being acceptable.
    bool shareEvents_;          ///< Determines whether the property allows further propagation of an accepted event.
};

//-----------------------------------------------------------------------------------------------------------------------

/**
 * EventProperties are the usual way of accessing user-generated events in Voreen.
 *
 * Each EventProperty is assigned to an owner, usually a Processor, that has to specify the type
 * and configuration of events it intends to receive as well as the member function these events
 * will be passed to. The property may be configured to relay mouse events, key events
 * or both event types. EventProperties allow a Processor author to get access to user-generated
 * events in a convenient way, by performing three basic steps:
 *  - add a function for handling the desired events to the Processor
 *  - declare and create an EventProperty member object
 *  - announce the EventProperty object by passing it to Processor::addEventProperty
 *
 * @tparam T the type of the object whose member function is passed as event receiver
 *
 * @note An alternative way of receiving events is to overwrite Processor::onEvent.
 *       However, the usage of EventProperties is strongly recommended due to several reasons:
 *       - The standardized way of event access provided (or enforced) by the EventProperty class
 *         alleviates the development of user interfaces for the configuration of the event handling
 *         within a processor network, e.g. enabling/disabling certain event receivers or remapping input keys.
 *       - The set of event properties owned by a network's processors contains a significant amount
 *         of information about the event handling in the network, such as which processor reacts to which events.
 *       - EventProperties encourage to group the local event handling of a Processor into logical units,
 *         e.g. there may be one EventProperty whose associated member function handles a trackball navigation,
 *         while another EventProperty/member function pair is used for picking operations.
 *       For an example user interface exploiting the EventProperty concept for providing documentation
 *       as well as customization of the network event handling, see InputMappingDialog.
 *
 * @note Accepting an relayed event is up to the event property's owner. @see setSharing
 *
 * @see Processor::addEventProperty
 *
 * Implementation remarks: Any event generated by a canvas as reaction to a user interaction is propagated
 * through the network, starting at the canvas. A processor propagates each received event to
 * its EventProperties, before further propagating it to its predecessors in the network (Processor::onEvent).
 *
 */
template<class T>
class EventProperty : public EventPropertyBase {

public:
    /**
     * Constructor creating an EventProperty that only listens to \e mouse events.
     *
     * @param id property identifier. Must be unique among the
     *  properties of a PropertyOwner and must not be empty.
     * @param guiName string to be displayed in a user interface
     * @param target the object whose member functions receives the matching events. Must not be null.
     * @param fptMouseEvent function pointer of the function that receives the matching events. Must not be null.
     * @param mouseButtons mouse buttons the property reacts to
     * @param mouseAction mouse actions the property reacts to
     * @param modifier modifier that has to pressed for accepting the event
     * @param shareEvents if set to true, the property does ignore the event,
     *  even if it matches and is passed to its owner. This causes the event
     *  to be further propagated, even if the event property has accepted it. @see setEnabled
     * @param enabled if set to false, the property does not react to any event
     */
    EventProperty(const std::string& id, const std::string& guiName,
        T* target, void (T::*fptMouseEvent)(tgt::MouseEvent*),
        tgt::MouseEvent::MouseButtons mouseButtons, int mouseAction = tgt::MouseEvent::ACTION_ALL,
        tgt::Event::Modifier modifier = tgt::Event::MODIFIER_NONE,
        bool shareEvents = false, bool enabled = true);

    /**
     * Constructor creating an EventProperty that only listens to \e key events.
     *
     * @param id property identifier. Must be unique among the
     *  properties of a PropertyOwner and must not be empty.
     * @param guiName string to be displayed in a user interface
     * @param target the object whose member functions receives the matching events. Must not be null.
     * @param fptKeyEvent function pointer of the function that receives the matching events. Must not be null.
     * @param keyCode code of the key the property reacts to
     * @param modifier modifier that has to pressed for accepting the event
     * @param shareEvents if set to true, the property does ignore the event,
     *  even if it matches and is passed to its owner. This causes the event
     *  to be further propagated, even if the event property has accepted it. @see setEnabled
     * @param enabled if set to false, the property does not react to any event
     */
    EventProperty(const std::string& id, const std::string& guiName,
        T* target, void (T::*fptKeyEvent)(tgt::KeyEvent*),
        tgt::KeyEvent::KeyCode keyCode,
        tgt::Event::Modifier modifier = tgt::Event::MODIFIER_NONE,
        bool shareEvents = false, bool enabled = true);

    /**
     * Constructor creating an EventProperty that listens to both \e mouse and \e key events.
     *
     * @param id property identifier. Must be unique among the
     *  properties of a PropertyOwner and must not be empty.
     * @param guiName string to be displayed in a user interface
     * @param target the object whose member functions receives the matching events. Must not be null.
     * @param fptEvent function pointer of the function that receives the matching events. Must not be null.
     * @param mouseButtons mouse buttons the property reacts to
     * @param mouseAction mouse actions the property reacts to
     * @param keyCode code of the key the property reacts to
     * @param modifier modifier that has to pressed for accepting the event
     * @param shareEvents if set to true, the property does ignore the event,
     *  even if it matches and is passed to its owner. This causes the event
     *  to be further propagated, even if the event property has accepted it. @see setEnabled
     * @param enabled if set to false, the property does not react to any event
     */
    EventProperty(const std::string& id, const std::string& guiName,
        T* target, void (T::*fptEvent)(tgt::Event*),
        tgt::MouseEvent::MouseButtons mouseButtons, int mouseAction,
        tgt::KeyEvent::KeyCode keyCode,
        tgt::Event::Modifier modifier = tgt::Event::MODIFIER_NONE,
        bool shareEvents = false, bool enabled = true);

    /**
     * Constructor creating an EventProperty that only listens to \e touch events.
     *
     * @param id property identifier. Must be unique among the
     *  properties of a PropertyOwner and must not be empty.
     * @param guiName string to be displayed in a user interface
     * @param target the object whose member functions receives the matching events. Must not be null.
     * @param fptTouchEvent function pointer of the function that receives the matching events. Must not be null.
     * @param shareEvents if set to true, the property does ignore the event,
     *  even if it matches and is passed to its owner. This causes the event
     *  to be further propagated, even if the event property has accepted it. @see setEnabled
     * @param enabled if set to false, the property does not react to any event
     */
    EventProperty(const std::string& id, const std::string& guiName,
        T* target, void (T::*fptTouchEvent)(tgt::TouchEvent*),
        bool shareEvents = false, bool enabled = true);

    /**
     * Default constructor. Needed for serialization. Do not use directly!
     */
    EventProperty();

    virtual Property* create() const;

    /**
     * Registers a change listener function that is called
     * after the property state (keys/buttons) has changed.
     *
     * @note The passed function pointer must point to
     *      a function of the target type that has been
     *      passed to the constructor. Only one listener
     *      can be registered. If this function is called
     *      multiple times, a previous listener is overwritten.
     */
    void onChange(void (T::*fptOnChange)());

protected:
    /**
     * Called by Processor::onEvent of the owning processor for
     * relaying incoming events.
     *
     * Internally, queries accepts(e) and passes the event
     * to the specified function, if true is returned.
     */
    virtual void execute(tgt::Event* e);

    /**
     * Implementation of EventPropertyBase::notifyChangeListener,
     * just calls notifyChangeListenerImpl.
     */
    virtual void notifyChangeListener();

private:
    /// Actual generic notification implementation.
    void notifyChangeListenerImpl();

    T* target_;                                   ///< The object matching events are relayed to
    void (T::*fptMouseEvent_)(tgt::MouseEvent*);  ///< Target function in case of an mouse-only property
    void (T::*fptKeyEvent_)(tgt::KeyEvent*);      ///< Target function in case of an key-only property
    void (T::*fptTouchEvent_)(tgt::TouchEvent*);  ///< Target function in case of an touch-only property
    void (T::*fptEvent_)(tgt::Event*);            ///< Target function in case of an property receiving all types

    void (T::*fptOnChange_)();                    ///< Called in case the property state (keys/buttons) changes
};

//-----------------------------------------------------------------------------------------------------------------------
// template definitions

template<class T>
EventProperty<T>::EventProperty(const std::string& id, const std::string& guiName,
                                T* target, void (T::*fptMouseEvent)(tgt::MouseEvent*),
                                tgt::MouseEvent::MouseButtons mouseButtons, int mouseAction,
                                tgt::Event::Modifier modifier,
                                bool shareEvents, bool enabled)
    : EventPropertyBase(id, guiName,
        true, false,
        mouseButtons, mouseAction,
        tgt::KeyEvent::K_UNKNOWN,
        modifier, shareEvents, enabled)
    , target_(target)
    , fptMouseEvent_(fptMouseEvent)
    , fptKeyEvent_(0)
    , fptTouchEvent_(0)
    , fptEvent_(0)
    , fptOnChange_(0)
{
    tgtAssert(target_ && fptMouseEvent_, "Passed target or function pointer invalid");
}

template<class T>
EventProperty<T>::EventProperty(const std::string& id, const std::string& guiName,
                                T* target, void (T::*fptKeyEvent)(tgt::KeyEvent*),
                                tgt::KeyEvent::KeyCode keyCode,
                                tgt::Event::Modifier modifier,
                                bool shareEvents, bool enabled)
    : EventPropertyBase(id, guiName,
        false, true,
        tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::ACTION_NONE,
        keyCode,
        modifier, shareEvents, enabled)
    , target_(target)
    , fptMouseEvent_(0)
    , fptKeyEvent_(fptKeyEvent)
    , fptTouchEvent_(0)
    , fptEvent_(0)
    , fptOnChange_(0)
{
    tgtAssert(target_ && fptKeyEvent_, "Passed target or function pointer invalid");
}

template<class T>
EventProperty<T>::EventProperty(const std::string& id, const std::string& guiName,
                                T* target, void (T::*fptTouchEvent)(tgt::TouchEvent*),
                                bool shareEvents, bool enabled)
    : EventPropertyBase(id, guiName,
        false, false,
        tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::ACTION_NONE,
        tgt::KeyEvent::K_UNKNOWN, tgt::Event::MODIFIER_NONE,
        shareEvents, enabled)
    , target_(target)
    , fptMouseEvent_(0)
    , fptKeyEvent_(0)
    , fptTouchEvent_(fptTouchEvent)
    , fptEvent_(0)
    , fptOnChange_(0)
{
    // TODO pass this via constructor
    receivesTouchEvents_ = true;
    tgtAssert(target_ && fptTouchEvent_, "Passed target or function pointer invalid");
}

template<class T>
EventProperty<T>::EventProperty(const std::string& id, const std::string& guiName,
                                T* target, void (T::*fptEvent)(tgt::Event*),
                                tgt::MouseEvent::MouseButtons mouseButtons, int mouseAction,
                                tgt::KeyEvent::KeyCode keyCode,
                                tgt::Event::Modifier modifier,
                                bool shareEvents, bool enabled)
    : EventPropertyBase(id, guiName,
        true, true,
        mouseButtons, mouseAction,
        keyCode,
        modifier, shareEvents, enabled)
    , target_(target)
    , fptMouseEvent_(0)
    , fptKeyEvent_(0)
    , fptTouchEvent_(0)
    , fptEvent_(fptEvent)
    , fptOnChange_(0)
{
    tgtAssert(target_ && fptEvent_, "Passed target or function pointer invalid");
}

template<class T>
voreen::EventProperty<T>::EventProperty()
    : EventPropertyBase("", "", false, false,
      tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::ACTION_NONE,
      tgt::KeyEvent::K_UNKNOWN, tgt::Event::MODIFIER_NONE, false, false)
    , target_(0)
    , fptMouseEvent_(0)
    , fptKeyEvent_(0)
    , fptTouchEvent_(0)
    , fptEvent_(0)
    , fptOnChange_(0)
{}

template<class T>
Property* voreen::EventProperty<T>::create() const {
    return new EventProperty<T>();
}

template<class T>
void voreen::EventProperty<T>::onChange( void (T::*fptOnChange)() ) {
    fptOnChange_ = fptOnChange;
}

template<class T>
void EventProperty<T>::execute(tgt::Event* e) {

    if (accepts(e)) {

        // mouse event
        if (tgt::MouseEvent* mouseEvent = dynamic_cast<tgt::MouseEvent*>(e)) {
            tgtAssert(fptMouseEvent_ || fptEvent_, "No valid function pointer");
            if (fptMouseEvent_)
                (target_->*fptMouseEvent_)(mouseEvent);
            else if (fptEvent_)
                (target_->*fptEvent_)(e);
        }

        // key event
        if (tgt::KeyEvent* keyEvent = dynamic_cast<tgt::KeyEvent*>(e)) {
            tgtAssert(fptKeyEvent_ || fptEvent_, "No valid function pointer");
            if (fptKeyEvent_)
                (target_->*fptKeyEvent_)(keyEvent);
            else if (fptEvent_)
                (target_->*fptEvent_)(e);
        }

        // touch event
        if (tgt::TouchEvent* touchEvent = dynamic_cast<tgt::TouchEvent*>(e)) {
            tgtAssert(fptTouchEvent_ || fptEvent_, "No valid function pointer");
            if (fptTouchEvent_)
                (target_->*fptTouchEvent_)(touchEvent);
            else if (fptEvent_)
                (target_->*fptEvent_)(e);
        }

        // in principal accepting the event is up to the owner,
        // but if sharing is enabled the event is ignored by force
        if (isSharing())
            e->ignore();
    }
}

template<class T>
void EventProperty<T>::notifyChangeListener() {
    notifyChangeListenerImpl();
}

template<class T>
void EventProperty<T>::notifyChangeListenerImpl() {
    updateWidgets();
    tgtAssert(target_, "No target");
    if (fptOnChange_)
        (target_->*fptOnChange_)();
}

} // namespace

#endif // VRN_EVENTPROPERTY_H
