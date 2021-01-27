/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
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

#ifndef TGT_TOUCHEVENT_H
#define TGT_TOUCHEVENT_H

#include "tgt/event/event.h"
#include "tgt/vector.h"
#include "tgt/event/touchpoint.h"
#include "tgt/types.h"
#include <deque>
//#include "tgt/qt/qtcanvas.h"
//#include "tgt/qt/qtapplication.h"

namespace tgt {

class TGT_API TouchEvent : public Event {
public:
    enum DeviceType{
        TouchScreen,
        TouchPad
    };

    TouchEvent(Event::Modifier mod, TouchPoint::State touchPointStates, DeviceType deviceType, const std::deque<TouchPoint> &touchPoints, const tgt::ivec2& viewport = tgt::ivec2(-1));
    ~TouchEvent();
    DeviceType deviceType() const;
    TouchPoint::State touchPointStates() const;

    const std::deque<TouchPoint> &touchPoints() const;

    Event::Modifier modifiers() const;

    virtual int getEventType();
    const tgt::ivec2& getViewport() {
        return viewport_;
    }

private:
    DeviceType deviceType_;
    TouchPoint::State touchPointStates_;
    const std::deque<TouchPoint> &touchPoints_;

protected:
    ivec2 coord_;
    ivec2 viewport_;
    Event::Modifier mod_;

};

}

#endif //TGT_MOUSEEVENT_H
