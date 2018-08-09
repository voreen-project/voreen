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

#ifndef TGT_STOPWATCH_H
#define TGT_STOPWATCH_H

#include <ctime>
#ifdef WIN32
    #include <windows.h>
#else
    #include <sys/time.h>
#endif

#include "tgt/types.h"

namespace tgt {

/**
 *   This class can be used to time events within a program.  For example,
 *   a framecounter can be easily implemented by letting the stopwatch run until
 *   it reaches one second and, while it's running, incrementing a counting
 *   variable each frame which then yields the fps.  A stopwatch can be stopped and
 *   started multiple times; the runtime is not affected by pauses in which the
 *   stopwatch was stopped.
 */
class TGT_API Stopwatch {
protected:
    /// is the stopwatch currently running or not?
    bool isStopped_;
    /// this is the point of time the stopwatch was started
    clock_t startTime_;
    /// this is the amount of time the stopwatch has been running until now
    clock_t runTime_;

public:
    /// The constructor.  It sets all the uints to 0 and stops the stopwatch.
    Stopwatch();

    /// This function should return the number of elapsed milliseconds since some reference point in the past.
    /// Under windows this is system start time. Be aware: it counts something like 47 days.
    ///
    /// DO NOT CONFUSE WITH getRuntime() !!
    static uint64_t getTicks();

    /// This is the function that returns the current runtime of the stopwatch in milliseconds.
    clock_t getRuntime();

    /// Stop the stopwatch.
    void stop();

    /// Start the stopwatch.
    void start();

    /// Returns the time the stopwatch was started.
    inline clock_t getStartTime() const { return startTime_; }

    /// Reset the stopwatch.
    void reset();
};

}
#endif  // TGT_STOPWATCH_H
