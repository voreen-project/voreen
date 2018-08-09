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

#ifndef TGT_NAVIGATION_H
#define TGT_NAVIGATION_H

#include "tgt/tgt_gl.h"

#include "tgt/camera.h"
#include "tgt/glcanvas.h"
#include "tgt/types.h"
#include "tgt/vector.h"
#include "tgt/event/mouseevent.h"
#include "tgt/event/eventlistener.h"

/**
    This is the base class for navigation metaphores.
    Derived classes offer solutions to use typical high-level camera movements in object space
    (e.g. moving camera upon a trackball by mouse, make tracking shot following a spline) in a
    simple way.

    In addition, class it implements basic camera operations like moving or rotating the camera
    in object space. Most deriving classes use these operations to do advanced navigations.
*/

namespace tgt {

class TGT_API Navigation : virtual public EventListener {

protected:

    Camera* camera_;

public:

    Navigation(Camera* cam) {
        camera_ = cam;
    }

    virtual ~Navigation() {}

    /// The following functions may be used to rotate the Camera about
    /// an arbitrary axis.
    void rotateView(float angle, float x, float y, float z);
    void rotateView(float angle, const vec3& axis);

    /// This function rotates the view about the up- and strafe-vector
    void rotateViewHV(float anglehorz, float anglevert);

    /// The following functions may be used to rotate the Camera-View about
    /// the Up- and the Strafe-Vector; they exist for reasons of convenience.
    void rotateViewVert(float angle);
    void rotateViewHorz(float angle);

    /// The following functions may be used to rotate the Up-Vector about
    /// the Strafe- and the Look-Vector.  Use this with care since it may
    /// leave the Camera with a "strange" orientation.
    void rollCameraVert(float angle);
    void rollCameraHorz(float angle);

    /// The following functions may be used to move the camera a given
    /// distance in a certain direction.
    void moveCameraForward(float length);
    void moveCameraBackward(float length);
    void moveCameraUp(float length);
    void moveCameraDown(float length);
    void moveCameraRight(float length);
    void moveCameraLeft(float length);
    void moveCamera(float length, float x, float y, float z);
    void moveCamera(float length, const vec3& axis);
    void moveCamera(const vec3& motionvector);

};

} // namespace tgt

#endif //TGT_NAVIGATION_H
