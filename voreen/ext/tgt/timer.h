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

#ifndef TGT_GUITIMER_H_
#define TGT_GUITIMER_H_

#include "tgt/types.h"

namespace tgt {

class EventHandler;

/// This is the base class for all tgt timer classes.
/// Derived classes use its functionality to to initiate gui-toolkit
/// dependend timers which gui-toolkit events are then translated
/// into tgt::TimeEvent objects and passed to an tgt::EventHandler.
class TGT_API Timer {
public:
    /// A constructor.
    /// @param eh       EventHandler to pass tgt::TimeEvent objects to.
    Timer( EventHandler* eh );

    /// A destructor
    virtual ~Timer();

    /// This abstract method is use to start this Timer. It will pass a tgt::TimeEvent to the
    /// actually set tgt::EventHandler.
    /// @param msec     Event will be created after given period of time in msec.
    /// @param limit    How many events to launch. Default is 0 (= unlimited).
    virtual void start( const int msec, const int limit = 0 ) = 0;

    /// Abstract method. Used to stop the timer.
    virtual void stop() = 0;


    // * getter / setter *

    /// Adjust the timer period
    virtual int getTickTime() { return tickTime_; };
    virtual void setTickTime( const int msec ) = 0;

    /// Redjust the counter limit.
    virtual void setLimit( const int limit ) { limit_ = limit; };
    virtual int getLimit() { return limit_; };

    /// Get the number of time events since calling tgt::Timer::start().
    virtual int getCount() { return count_; };
    virtual void setCount( const int count ) { count_ = count; };

    virtual bool isStopped() { return stopped_; }

    /// Set the tgt::EventHandler used to send tgt::TimeEvent objects to.
    void setEventHandler( EventHandler* eh ) { eventHandler_ = eh; }
    EventHandler* getEventHandler() { return eventHandler_; }

protected:
    /// How many events to fire till we stop.
    int limit_;
    /// Counts events fired since last start.
    int count_;
    /// How much time between the events to fire in msec.
    int tickTime_;
    /// Holds the tgt::EventHandler in which the tgt::TimeEvent
    /// objects are broadcasted.
    EventHandler* eventHandler_;
    bool stopped_;
};

}

#endif //TGT_GUITIMER_H
