/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/event/eventlistener.h"

#include "tgt/event/event.h"
#include "tgt/event/keyevent.h"
#include "tgt/event/mouseevent.h"
#include "tgt/event/timeevent.h"
#include "tgt/event/touchevent.h"
#include "tgt/event/touchpoint.h"
#include <typeinfo>

namespace tgt {

EventListener::EventListener() {
    eventTypes_ = Event::ALLEVENTTYPES;
}

EventListener::~EventListener() {
}

void EventListener::onEvent(Event* e) {
    /*
        since we already checked which type we have we can safely use
        the faster static_cast (rl)
    */
    if (typeid(*e) == typeid(MouseEvent)) {
        MouseEvent* me = static_cast<MouseEvent*>(e);
        if (me->action() == MouseEvent::PRESSED)
            mousePressEvent(me);
        else if (me->action() == MouseEvent::RELEASED)
            mouseReleaseEvent(me);
        else if (me->action() == MouseEvent::MOTION)
            mouseMoveEvent(me);
        else if (me->action() == MouseEvent::DOUBLECLICK)
            mouseDoubleClickEvent(me);
        else if (me->action() == MouseEvent::WHEEL)
            wheelEvent(me);
    }
    else if (typeid(*e) == typeid(KeyEvent)) {
        keyEvent(static_cast<KeyEvent*>(e));
    }
    else if (typeid(*e) == typeid(TimeEvent)) {
        timerEvent(static_cast<TimeEvent*>(e));
    }
    else if (typeid(*e) == typeid(TouchEvent)) {
        touchEvent(static_cast<TouchEvent*>(e));
    }
}

void EventListener::setEventTypes(int eventTypes){
    eventTypes_ = eventTypes;
}

int EventListener::getEventTypes() const{
    return eventTypes_;
}

void EventListener::mousePressEvent(MouseEvent* e) {
    e->ignore();
}

void EventListener::mouseReleaseEvent(MouseEvent* e) {
    e->ignore();
}

void EventListener::mouseMoveEvent(MouseEvent* e) {
    e->ignore();
}

void EventListener::mouseDoubleClickEvent(MouseEvent* e) {
    e->ignore();
}

void EventListener::wheelEvent(MouseEvent* e) {
    e->ignore();
}

void EventListener::timerEvent(TimeEvent* e) {
    e->ignore();
}

void EventListener::keyEvent(KeyEvent* e) {
    e->ignore();
}

void EventListener::touchEvent(TouchEvent* e) {
    e->ignore();
}


} // namespace tgt
