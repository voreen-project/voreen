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

#ifndef TGT_TOUCHPOINT_H
#define TGT_TOUCHPOINT_H

#include "tgt/event/event.h"
#include "tgt/vector.h"
//#include "tgt/qt/qtcanvas.h"
//#include "tgt/qt/qtapplication.h"

namespace tgt {
class TGT_API TouchPoint {
public:

    enum State {
        TouchPointPressed    = 0x01,
        TouchPointMoved      = 0x02,
        TouchPointStationary = 0x04,
        TouchPointReleased   = 0x08,
        TouchPointDblTouch   = 0x10
    };

    TouchPoint();
    ~TouchPoint();


    int    id () const;
    bool        isPrimary () const;
    vec2   lastNormalizedPos () const;
    vec2    lastPos () const;
    vec2    lastScenePos () const;
    vec2    lastScreenPos () const;
    vec2    normalizedPos () const;
    vec2    pos () const;
    double    pressure () const;
    //QRectF    rect () const;
    vec2    scenePos () const;
    //QRectF    sceneRect () const;
    vec2    screenPos () const;
    //QRectF    screenRect () const;
    vec2    startNormalizedPos () const;
    vec2    startPos () const;
    vec2    startScenePos () const;
    vec2    startScreenPos () const;
    State state() const;

    // internal
    void setId(int id);
    void setPrimary(bool primary);
    void setState(State state);
    void setPos(const vec2  &pos);
    void setScenePos(const vec2 &scenePos);
    void setScreenPos(const vec2 &screenPos);
    void setNormalizedPos(const vec2 &normalizedPos);
    void setStartPos(const vec2 &startPos);
    void setStartScenePos(const vec2 &startScenePos);
    void setStartScreenPos(const vec2 &startScreenPos);
    void setStartNormalizedPos(const vec2 &startNormalizedPos);
    void setLastPos(const vec2 &lastPos);
    void setLastScenePos(const vec2 &lastScenePos);
    void setLastScreenPos(const vec2 &lastScreenPos);
    void setLastNormalizedPos(const vec2 &lastNormalizedPos);
    //void setRect(const QRectF &rect);
    //void setSceneRect(const QRectF &sceneRect);
    //void setScreenRect(const QRectF &screenRect);
    void setPressure(double pressure_);
    TouchPoint &operator=(const TouchPoint &other_);

private:
    int id_;
    bool isPrimary_;
    State state_;
    vec2 pos_;
    vec2 scenePos_;
    vec2 screenPos_;
    vec2 normalizedPos_;
    vec2 startPos_;
    vec2 startScenePos_;
    vec2 startScreenPos_;
    vec2 startNormalizedPos_;
    vec2 lastPos_;
    vec2 lastScenePos_;
    vec2 lastScreenPos_;
    vec2 lastNormalizedPos_;
    double pressure_;
};

}

#endif //TGT_MOUSEEVENT_H
