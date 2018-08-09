/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
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

#include "tgt/event/eventhandler.h"

#include "tgt/assert.h"

#include <algorithm>
#include <vector>

using std::vector;

namespace tgt {

EventHandler::EventHandler()
  : listeners_(std::deque<EventListener*>()) {
}

void EventHandler::addListenerToBack(EventListener* e) {
    tgtAssert(e != 0, "Trying to add null pointer as event listener");
    listeners_.push_back(e);
}

void EventHandler::addListenerToFront(EventListener* e) {
    tgtAssert(e != 0, "Trying to add null pointer as event listener");
    listeners_.push_front(e);
}

void EventHandler::removeListener(EventListener* e) {
    std::deque<EventListener*>::iterator pos;
    pos = std::find(listeners_.begin(), listeners_.end(), e);

    if (pos != listeners_.end())
        listeners_.erase(pos);
}

void EventHandler::clear() {
    listeners_.clear();
}

void EventHandler::broadcast(Event* e) {
    for (size_t i = 0 ; i < listeners_.size() ; ++i) {
        // check if current listener listens to the eventType of e
        if(listeners_[i]->getEventTypes() & e->getEventType() ){
            listeners_[i]->onEvent(e);
            if (e->isAccepted())
                break;
        }
    }
    delete e;
}

size_t EventHandler::getListenerNumber() const {
    return listeners_.size();
}

}

