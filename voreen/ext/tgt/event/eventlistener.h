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

#ifndef TGT_EVENTLISTENER_H
#define TGT_EVENTLISTENER_H

#include "tgt/types.h"

namespace tgt {

class Event;
class KeyEvent;
class MouseEvent;
class TimeEvent;
class TouchEvent;

class TGT_API EventListener {
public:
    EventListener();
    virtual ~EventListener();

    virtual void mousePressEvent(MouseEvent* e);
    virtual void mouseReleaseEvent(MouseEvent* e);
    virtual void mouseMoveEvent(MouseEvent* e);
    virtual void mouseDoubleClickEvent(MouseEvent* e);
    virtual void wheelEvent(MouseEvent* e);
    virtual void timerEvent(TimeEvent* e);
    virtual void keyEvent(KeyEvent* e);
    virtual void touchEvent(TouchEvent* e);

    virtual void onEvent(Event* e);

    void setEventTypes(int eventTypes);
    int getEventTypes() const;

private:
    /*
     * bitmask in which all eventTypes this EventListener is listening to are stored
     * could be used to deactivate an EventListener by setting eventTypes_ = 0
     */
    int eventTypes_;
};

}

#endif //TGT_EVENTLISTENER_H
