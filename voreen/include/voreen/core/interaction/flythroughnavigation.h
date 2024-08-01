/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#ifndef _VOREEN_FLYTHROUGH_NAVIGATION_H_
#define _VOREEN_FLYTHROUGH_NAVIGATION_H_

#include "tgt/event/eventlistener.h"
#include "tgt/vector.h"
#include "tgt/event/mouseevent.h"
#include "tgt/glcanvas.h"
#include "tgt/navigation/navigation.h"

#include "voreen/core/utils/voreenpainter.h"

namespace voreen {

/**
 * A class that makes it possible to use a trackball-metaphor to rotate a dataset.
 */
class VRN_CORE_API FlythroughNavigation : public tgt::Navigation {

public:

    /**
     * The Constructor.
     *
     * @param canvas the canvas whose camera will be modified
     */
    FlythroughNavigation(tgt::GLCanvas* canvas, tgt::Camera* camera);
    virtual ~FlythroughNavigation() {}

    virtual void keyEvent(tgt::KeyEvent* e);

    /**
     * React to a press-Event.  In this case, remember the current mouse-position, so we can calculate relative movement in mouseMoveEvent.
     *
     * @param e The event to be processed.
     */
    virtual void mousePressEvent(tgt::MouseEvent* e);

    /**
     * React to a release-Event.  In this case, reset \a leftPressed_ if necessary and switch coarseness off.
     *
     * @param e The event to be processed.
     */
    virtual void mouseReleaseEvent(tgt::MouseEvent* e);

    /**
     * React to a move-Event.  In this case, this actually causes the object to rotate.
     *
     * @param e The event to be processed.
     */
    virtual void mouseMoveEvent(tgt::MouseEvent* e);

    /**
     * React to a double-click-Event.  In this case, this actually resets the camera's position and orientation
     *
     * @param e The event to be processed.
     */
    virtual void mouseDoubleClickEvent(tgt::MouseEvent* e);

   /**
     * React to a mouse-wheel-Event.  In this case, this actually causes the object to rotate sideways.
     *
     * @param e The event to be processed.
     */
    virtual void wheelEvent(tgt::MouseEvent* e);

   /**
     * React to a time-Event.  In this case, this does a number of things:
     *      - If necessary auto-spin the trackball
     *      - If the mouse-wheel was used, switch coarseness off after a certain amount of time
     *      - If auto-spinning of the trackball is active, don't spin the trackball if the user has waited too long after moving the mouse
     *
     * @param e The event to be processed.
     */
    virtual void timerEvent(tgt::TimeEvent* e);

protected:

    tgt::MouseEvent::MouseButtons resetButton_; ///< If this button is double-clicked, reset trackball to initial
                                                /// position and orientation

    int wheelCounter_; ///< Counts how many time-ticks have passed since the mouse-wheel was used.
    int moveCounter_;  ///< Counts how many time-ticks have passed since the user has moved the mouse.

    int wheelID_;
    int spinID_;
    int moveID_;
    int overlayTimerID_;

    bool initialCameraValid_;
    tgt::vec3 initialPosition_;
    tgt::vec3 initialFocusPoint_;
    tgt::vec3 initialUpVector_;

    tgt::GLCanvas* canvas_;

    static const std::string loggerCat_;

};

}

#endif //_VOREEN_FLYTHROUGH_NAVIGATION_H_
