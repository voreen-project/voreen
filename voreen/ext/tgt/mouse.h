/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#ifndef TGT_MOUSE_H
#define TGT_MOUSE_H

#include "tgt/vector.h"

namespace tgt {

/**
    This abstract class capsules the basic funcionality of a mouse.
*/
class Mouse
{
protected:

    /// the rate determines how sensible the mouse is
    float rate_;
    /// is the mouse-cursor currently hidden?
    bool showCursor_;

    /// these booleans capture the state of the mouse
    bool leftButtonPressed_;
    bool middleButtonPressed_;
    bool rightButtonPressed_;
    bool hasMoved_;

public:

    /// the Constructor
    Mouse(float rate, bool show) {
        rate_ = rate;
        showCursor_ = show;

        leftButtonPressed_   = false;
        middleButtonPressed_ = false;
        rightButtonPressed_  = false;
        hasMoved_            = false;
    }

    virtual ~Mouse() {};

    /// classes derived from Mouse have at least to implement the following
    /// functions (which are pretty much self-explanatory)
    virtual ivec2 getPosition() const = 0;
    virtual void setPosition(const ivec2& v) = 0;
    virtual void hideCursor() = 0;
    virtual void showCursor() = 0;

    /// switches Cursor on or off
    void toggleCursor() {
        showCursor_ ? hideCursor() : showCursor();
    }

    /// getter / setter
    float getRate() const {return rate_;}
    void  setRate(float rate) {rate_ = rate;}

    void setLeft  (bool enable) { leftButtonPressed_    = enable; }
    void setMiddle(bool enable) { middleButtonPressed_  = enable; }
    void setRight (bool enable) { rightButtonPressed_   = enable; }
    void setMoved (bool enable) { hasMoved_             = enable; }

    bool getLeft()   { return leftButtonPressed_;   }
    bool getMiddle() { return middleButtonPressed_; }
    bool getRight()  { return rightButtonPressed_;  }
    bool getMoved()  { return hasMoved_;            }
};

} // namespace tgt

#endif //TGT_MOUSE_H
