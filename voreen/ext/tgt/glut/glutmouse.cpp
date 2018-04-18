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

#include "tgt/glut/glutmouse.h"

namespace tgt {

GLUTMouse::GLUTMouse(float rate, bool show)
  : Mouse(rate, show)
{
    if (show)
        glutSetCursor( GLUT_CURSOR_INHERIT );
    else
        glutSetCursor( GLUT_CURSOR_NONE );
}

ivec2 GLUTMouse::getPosition() const {
    return position_;
}

void GLUTMouse::setPosition(const ivec2& pos) {
    position_ = pos;
    glutWarpPointer( pos.x, pos.y );
}

void GLUTMouse::hideCursor() {
    glutSetCursor( GLUT_CURSOR_INHERIT );
    showCursor_ = false;
}

void GLUTMouse::showCursor() {
    glutSetCursor( GLUT_CURSOR_NONE );
    showCursor_ = true;
}

}
