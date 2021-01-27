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

#ifndef TGT_QTTIMER_H_
#define TGT_QTTIMER_H_

#include "tgt/event/eventhandler.h"
#include "tgt/timer.h"
#include <QObject>

namespace tgt {

class QtTimer : public Timer, public QObject {
public:
    QtTimer( EventHandler* eh );
    virtual ~QtTimer();
    virtual void start( const int msec, const int limit = 0 );
    virtual void stop();
    virtual void setTickTime( const int msec );
    virtual void timerEvent( QTimerEvent* e );

private:
    int id_;
};

}

#endif //TGT_QTTIMER_H
