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

#include "tgt/event/keyevent.h"

namespace tgt {

KeyEvent::KeyEvent(KeyCode keyCode, int mod, bool autoRepeat, bool pressed)
    : Event()
    , keyCode_(keyCode)
    , mod_(mod)
    , autoRepeat_(autoRepeat)
    , pressed_(pressed)
{}

int KeyEvent::getEventType() {
    return KEYEVENT;
}

KeyEvent::KeyCode KeyEvent::keyCode() {
    return keyCode_;
}

int KeyEvent::modifiers() {
    return mod_;
}

bool KeyEvent::autoRepeat() {
    return autoRepeat_;
}

bool KeyEvent::pressed() {
    return pressed_;
}

}
