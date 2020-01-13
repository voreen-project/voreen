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

#include "tgt/event/touchpoint.h"

namespace tgt {

TouchPoint::TouchPoint() {

}

TouchPoint::~TouchPoint() {

}

int    TouchPoint::id () const {
    return id_; //< TODO
}

bool TouchPoint::isPrimary() const {
    return isPrimary_; //< TODO
}

tgt::vec2 TouchPoint::lastNormalizedPos() const {
    return lastNormalizedPos_;
}

tgt::vec2 TouchPoint::lastPos() const {
    return lastPos_;
}

tgt::vec2 TouchPoint::lastScenePos() const {
    return lastScenePos_;
}

tgt::vec2 TouchPoint::lastScreenPos() const {
    return lastScreenPos_;
}

tgt::vec2 TouchPoint::normalizedPos() const {
    return normalizedPos_;
}

tgt::vec2 TouchPoint::pos() const {
    return pos_;
}

double TouchPoint::pressure() const {
    return pressure_;
}

/*QRectF    TouchPoint::rect() const {

}*/

tgt::vec2 TouchPoint::scenePos() const {
    return scenePos_;
}

/*QRectF    TouchPoint::sceneRect() const {

}*/

tgt::vec2 TouchPoint::screenPos() const {
    return screenPos_;
}

/*QRectF    TouchPoint::screenRect() const {

}*/

tgt::vec2 TouchPoint::startNormalizedPos() const {
    return startNormalizedPos_;
}

tgt::vec2 TouchPoint::startPos() const {
    return startPos_;
}

tgt::vec2 TouchPoint::startScenePos() const {
    return startScenePos_;
}

tgt::vec2 TouchPoint::startScreenPos() const {
    return startScreenPos_;
}

TouchPoint::State    TouchPoint::state() const {
    return state_;
}

//Internal
void TouchPoint::setId(int id){ id_ = id; }
void TouchPoint::setPrimary(bool primary){ isPrimary_ = primary; }
void TouchPoint::setState(TouchPoint::State state){ state_ = state; }
void TouchPoint::setPos(const vec2  &pos){ pos_ = pos; }
void TouchPoint::setScenePos(const vec2 &scenePos){ scenePos_ = scenePos; }
void TouchPoint::setScreenPos(const vec2 &screenPos){ screenPos_ = screenPos; }
void TouchPoint::setNormalizedPos(const vec2 &normalizedPos){normalizedPos_ = normalizedPos;}
void TouchPoint::setStartPos(const vec2 &startPos){ startPos_ = startPos; }
void TouchPoint::setStartScenePos(const vec2 &startScenePos){ startScenePos_ = startScenePos; }
void TouchPoint::setStartScreenPos(const vec2 &startScreenPos){ startScreenPos_ = startScreenPos; }
void TouchPoint::setStartNormalizedPos(const vec2 &startNormalizedPos){ startNormalizedPos_ = startNormalizedPos; }
void TouchPoint::setLastPos(const vec2 &lastPos){ lastPos_ = lastPos; }
void TouchPoint::setLastScenePos(const vec2 &lastScenePos){ lastScenePos_ = lastScenePos; }
void TouchPoint::setLastScreenPos(const vec2 &lastScreenPos){ lastScreenPos_ = lastScreenPos; }
void TouchPoint::setLastNormalizedPos(const vec2 &lastNormalizedPos){ lastNormalizedPos_ = lastNormalizedPos; }
//void setRect(const QRectF &rect);
//void setSceneRect(const QRectF &sceneRect);
//void setScreenRect(const QRectF &screenRect);
void TouchPoint::setPressure(double pressure){ pressure_ = pressure; }

}

