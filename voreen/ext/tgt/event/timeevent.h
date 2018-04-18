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

#ifndef TGT_TIMEEVENT_H
#define TGT_TIMEEVENT_H

#include "tgt/event/event.h"
#include "tgt/types.h"

namespace tgt {

class Timer;

class TGT_API TimeEvent : public Event {
public:
    /// A default constructor. Note that there is no tgt::Timer set this way.
    TimeEvent();
    TimeEvent(Timer* timer);

    /// Getter for the tgt::Timer that throw this Event
    /// @return tgt::Timer
    Timer* getTimer();

    /// Sets the tgt::Timer that is used to determine who sent this event.
    void setTimer(Timer* timer);

    virtual int getEventType();

protected:
    Timer* timer_;
};

}

#endif //TGT_TIMEEVENT_H
