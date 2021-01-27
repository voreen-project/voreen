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

#ifndef TGT_GLUTMOUSE_H
#define TGT_GLUTMOUSE_H

#include "tgt/mouse.h"
#include <GL/glut.h>


/**
    This class implements the abstract class Mouse by using GLUT-Functions.
*/
namespace tgt {

class GLUTMouse : public Mouse
{
protected:
    ivec2 position_;

public:

    /// Constructor / Destructor
    GLUTMouse(float rate = 5.f, bool show = true);
    ~GLUTMouse() {};

    /// These methods implement the abstract functions of super-class Mouse
    virtual ivec2 getPosition() const;
    virtual void setPosition(const ivec2& pos);
    virtual void hideCursor();
    virtual void showCursor();
};

}

#endif //tgt_GLUTMOUSE_H
