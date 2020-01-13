/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#include "tgt/event/touchevent.h"

namespace tgt {


TouchEvent::TouchEvent(Event::Modifier mod, TouchPoint::State touchPointStates, DeviceType deviceType, const std::deque<TouchPoint> &touchPoints, const tgt::ivec2& viewport)
    : Event()
    , deviceType_(deviceType)
    , touchPointStates_(touchPointStates)
    , touchPoints_(touchPoints)
    , viewport_(viewport)
    , mod_(mod)
{
}

TouchEvent::~TouchEvent() {

}

TouchEvent::DeviceType TouchEvent::deviceType() const
{
    return deviceType_;
}

TouchPoint::State TouchEvent::touchPointStates() const
{
    return touchPointStates_;
}

const std::deque<TouchPoint> &TouchEvent::touchPoints() const
{
    return touchPoints_;
}

int TouchEvent::getEventType() {
    return TOUCHEVENT;
}

Event::Modifier TouchEvent::modifiers() const {
    return mod_;
}

}
