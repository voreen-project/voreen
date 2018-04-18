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

#ifndef TGT_GLUTTIMER_H
#define TGT_GLUTTIMER_H

#include "tgt/timer.h"

namespace tgt {

class EventHandler;

class GLUTTimer : public Timer {
public:
    // FIXME: strange comment ...
    /// used to hold the canvases registered in our app.
    /// must be global for the glut wrapper.
    /// see GLUTApplication for more details
    static GLUTTimer* GLUTTimers_[];
    static int GLUTTimerCounter_;

    GLUTTimer( EventHandler* eh );
    virtual ~GLUTTimer();
    virtual void start( const int msec, const int limit = 0 );
    virtual void stop();
    virtual void setTickTime( const int msec );
    virtual void timerEvent();
    static void timerEventFunc( int id );

protected:
    int id_;
};

} // namespace tgt

#endif //TGT_GLUTTIMER_H
