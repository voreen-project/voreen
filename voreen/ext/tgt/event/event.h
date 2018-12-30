/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
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

#ifndef TGT_EVENT_H
#define TGT_EVENT_H

#include "tgt/types.h"



namespace tgt {
/**
 *  The base class for all events. Provides the accepted/ignore mechanism.
 */
class TGT_API Event {
public:
    enum Modifier {
        MODIFIER_NONE   = 0x0000,
        SHIFT           = 0x0001,
        CTRL            = 0x0002,
        ALT             = 0x0004,
        META            = 0x0008,
        NUM             = 0x0010,
        CAPS            = 0x0020,
        MODE            = 0x0040,
        RESERVED        = 0x8000
    };

    /*
     * enum to select to which eventType has to be broadcast, see EventHandler::broadcast
     */
    enum eventType {
        MOUSEPRESSEVENT       = 0x0001,
        MOUSERELEASEEVENT     = 0x0002,
        MOUSEMOVEEVENT        = 0x0004,
        MOUSEDOUBLECLICKEVENT = 0x0008,
        MOUSEENTEREVENT       = 0x0010,
        MOUSEEXITEVENT        = 0x0020,
        WHEELEVENT            = 0x0040,
        TIMEREVENT            = 0x0080,
        KEYEVENT              = 0x0100,
        ROTATIONEVENT         = 0x0200,
        ZOOMEVENT             = 0x0400,
        SIZINGEVENT           = 0x0800,
        TRANSLATIONEVENT      = 0x1000,
        TWOFINGEREVENT        = 0x2000,
        TOUCHEVENT            = 0x4000,

        ALLEVENTTYPES = MOUSEPRESSEVENT | MOUSERELEASEEVENT | MOUSEMOVEEVENT | MOUSEDOUBLECLICKEVENT | MOUSEENTEREVENT | MOUSEEXITEVENT | WHEELEVENT | TIMEREVENT | KEYEVENT | ROTATIONEVENT | ZOOMEVENT | SIZINGEVENT | TRANSLATIONEVENT | TWOFINGEREVENT | TOUCHEVENT
    };


    /**
     * The Constructor. The use of the \a accepted_ variable is modeled after qt-events; \a accepted_ being on true
     * as the default-value just opens the possibility to pass on the event using ignore().
     */
    Event() {
        accepted_ =  true;
    }

    virtual ~Event() {}

    bool isAccepted() {
        return accepted_;
    }
    void accept() {
        accepted_ = true;
    }

    void ignore() {
        accepted_ = false;
    }

    virtual int getEventType() = 0;

    bool accepted_;

};

}

#endif //TGT_EVENT_H
