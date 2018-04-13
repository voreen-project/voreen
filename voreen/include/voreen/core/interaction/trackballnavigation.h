/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_TRACKBALL_NAVIGATION_H
#define VRN_TRACKBALL_NAVIGATION_H

#include "tgt/event/eventlistener.h"
#include "tgt/vector.h"
#include "tgt/event/mouseevent.h"
#include "tgt/event/touchevent.h"
#include "tgt/event/touchpoint.h"
#include "tgt/glcanvas.h"
#include "voreen/core/interaction/voreentrackball.h"

#include "voreen/core/utils/voreenpainter.h"

namespace tgt {
class Camera;
}

namespace voreen {

class CameraProperty;

/**
 * A class that makes it possible to use a trackball-metaphor to rotate, zoom, shift, roll a dataset.
 */
class VRN_CORE_API TrackballNavigation : public tgt::EventListener {
public:

    enum Mode {
        ROTATE_MODE,
        ZOOM_MODE,
        SHIFT_MODE,
        ROLL_MODE
    };

    /**
     * The Constructor.
     *
     * @param camera the camera property that is to be modified by the navigation
     * @param minDist the minimum allowed orthogonal distance to the center of the trackball
     * @param maxDist the maximum allowed orthogonal distance to the center of the trackball
     */
    TrackballNavigation(CameraProperty* cameraProperty, Mode mode = ROTATE_MODE, float minDist = 0.01f);
    virtual ~TrackballNavigation();

    void setMode(Mode mode);

    Mode getMode() const;

    bool isTracking() {
        return tracking_;
    }

    /**
     * React to a press-Event.
     *
     * @param e The event to be processed.
     */
    virtual void mousePressEvent(tgt::MouseEvent* e);

     /**
     * React to a touch press-Event.
     *
     * @param e The event to be processed.
     */

    virtual void touchPressEvent(tgt::TouchEvent* e);

     /**
     * React to a touch release event.
     *
     * @param e The event to be processed.
     */

    virtual void touchReleaseEvent(tgt::TouchEvent* e);

    /**
     * React to a touch move-Event. Actually this causes rotation or zoom.
     *
     * @param e The event to be processed.
     */
    virtual void touchMoveEvent(tgt::TouchEvent* e);

    /**
     * React to a release event.
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
     * React to a double-click-Event.
     *
     * @param e The event to be processed.
     */
    virtual void mouseDoubleClickEvent(tgt::MouseEvent* e);

   /**
     * React to a mouse-wheel-Event.
     *
     * @param e The event to be processed.
     */
    virtual void wheelEvent(tgt::MouseEvent* e);

   /**
     * React to a time-Event. This does a number of things:
     *
     *  - If necessary auto-spin the trackball.
     *  - If the mouse-wheel was used, switch coarseness off after a certain amount of time.
     *  - If auto-spinning of the trackball is active, don't spin the trackball if the user has
     *    waited too long after moving the mouse.
     *
     * @param e The event to be processed.
     */
    virtual void timerEvent(tgt::TimeEvent* e);

    /**
     * React to key event.
     */
    virtual void keyEvent(tgt::KeyEvent* e);

    /**
     * Returns the internally used trackball.
     */
    VoreenTrackball* getTrackball();

    /**
     *  Let trackball react on key pressures with rotating
     *  @param acuteness The per-keypress angle of rotation will be smaller at greater acuteness. Use acuteness = 0.f to disable key rotation.
     *  @param left, right, up, down keycode which should cause suitable rotation
     *  @param mod which modifiers must be set
     *  @param pressed rotate on key-pressed-event when true, on key-release-event when false
     *
     */
    void setKeyRotate(float acuteness = 10.f,
                      tgt::KeyEvent::KeyCode left  = tgt::KeyEvent::K_LEFT,
                      tgt::KeyEvent::KeyCode right = tgt::KeyEvent::K_RIGHT,
                      tgt::KeyEvent::KeyCode up    = tgt::KeyEvent::K_UP,
                      tgt::KeyEvent::KeyCode down  = tgt::KeyEvent::K_DOWN,
                      int mod = tgt::Event::MODIFIER_NONE,
                      bool pressed = false);


    /**
     * Let trackball react on key presses with moving
     * @param acuteness The per-keypress length of movement will be smaller at greater
     *        acuteness. Use acuteness = 0.f to disable key rotation.
     * @param left, right, up, down keycode which should cause suitable movement
     * @param mod which modifiers must be set
     * @param pressed move on key-pressed-event when true, on key-release-event when false
     */
    void setKeyMove(float acuteness = 100.f,
                    tgt::KeyEvent::KeyCode left  = tgt::KeyEvent::K_LEFT,
                    tgt::KeyEvent::KeyCode right = tgt::KeyEvent::K_RIGHT,
                    tgt::KeyEvent::KeyCode up    = tgt::KeyEvent::K_UP,
                    tgt::KeyEvent::KeyCode down  = tgt::KeyEvent::K_DOWN,
                    int mod = tgt::Event::SHIFT,
                    bool pressed = false);


    /**
     * Defines the mouse zoom direction.
     *
     * @param zoomInDirection specifies in which direction mouse must be moved to zoom in.
     *        The greater the length of this vector, the longer the way the mouse must be moved
     *        to achieve a certain zoom factor.
     */
    void setMouseZoom(tgt::vec2 zoomInDirection);

    /**
     * Defines mouse wheel zoom acuteness.
     *
     * @param acuteness The zoom factor will be smaller at greater acuteness. Use acuteness = 0.f to disable key rotation.
     * @param wheelUpZoomIn zoom in on wheel up and out on wheel down when true, otherwise when false
     */
    void setMouseWheelZoom(float acuteness = 10.f, bool wheelUpZoomIn = true);

    /// @param acuteness The zoom factor will be smaller at greater acuteness. Use acuteness = 0.f to disable key rotation.
    void setKeyZoom(float acuteness = 10.f,
                    tgt::KeyEvent::KeyCode in  = tgt::KeyEvent::K_UP,
                    tgt::KeyEvent::KeyCode out = tgt::KeyEvent::K_DOWN,
                    int mod = tgt::Event::CTRL,
                    bool pressed = false);

    /**
     * Defines the mouse roll acuteness.
     *
     * @param acuteness greater the acuteness means less tiling
     */
    void setMouseRoll(float acuteness);

    /**
     * Defines the mouse wheel roll acuteness and direction.
     *
     * @param acuteness Less tiling at greater acuteness. Use acuteness = 0.f to disable key rotation.
     * @param wheelUpRollLeft roll left on wheel up and right on wheel down when true, otherwise when false
     */
    void setMouseWheelRoll(float acuteness, bool wheelUpRollLeft = true);

    /// @param acuteness Less tiling at greater acuteness. Use acuteness = 0.f to disable key rotation.
    void setKeyRoll(float acuteness = 10.f,
                    tgt::KeyEvent::KeyCode left = tgt::KeyEvent::K_LEFT,
                    tgt::KeyEvent::KeyCode right = tgt::KeyEvent::K_RIGHT,
                    int mod = tgt::Event::ALT,
                    bool pressed = false);

protected:

    void startMouseDrag(tgt::MouseEvent* eve);
    void endMouseDrag(tgt::MouseEvent* eve);

    /// The following functions may be used to rotate the Up-Vector about
    /// the Strafe- and the Look-Vector.  Use this with care since it may
    /// leave the Camera with a "strange" orientation.
    void rollCameraVert(float angle);
    void rollCameraHorz(float angle);

    void initializeEventHandling();

    /// scale screen-coodinates of mouse to intervall [-1, 1]x[-1, 1]
    tgt::vec2 scaleMouse(const tgt::ivec2& coords, const tgt::ivec2& viewport) const;

    float getRotationAngle(float acuteness) const;
    float getMovementLength(float acuteness) const;
    float getZoomFactor(float acuteness, bool zoomIn) const;
    float getRollAngle(float acuteness, bool left) const;

    /// Resets position and orientation of the trackball's camera to the initial parameters the camera had
    /// when passed to the trackball.
    /// Projective parameters (frustum) are not touched.
    virtual void resetTrackball();

    CameraProperty* cameraProperty_;     ///< Camera property that is modified
    VoreenTrackball* trackball_;        ///< The trackball that is modified when navigating

    Mode mode_;     ///< current trackball mode: rotate, zoom, shift, roll

    /// Last mouse coordinates to allow calculation of the relative change.
    /// Ought to be relative coordinates within range [-1, 1]x[-1, 1].
    tgt::vec2 lastMousePosition_;

    // Last distance between two TouchPoints to allow calculation of zoomFactor
    float lastDistance_;

    // Last connection vector between two touch points
    tgt::vec2 lastConnection_;

    //float minDistance_;     ///< minimal allowed orthogonal distance to center of trackball
    //float maxDistance_;     ///< maximal allowed orthogonal distance to center of trackball (now retrieved from camera property)

    int wheelCounter_; ///< Counts how many time-ticks have passed since the mouse-wheel was used.
    int spinCounter_;  ///< Counts how many time-ticks have passed since the trackball was spun.
    int moveCounter_;  ///< Counts how many time-ticks have passed since the user has moved the mouse.

    int wheelID_;
    int spinID_;
    int moveID_;
    int overlayTimerID_;

    bool spinit_; ///< Helper member to control auto-spinning of trackball.
    bool trackballEnabled_; ///< Is the trackball enabled?

    bool tracking_; ///< Are we tracking mouse move events? Only when we received a mousePressEvent before.

    /// used to store settings how to react on certain events
    float keyRotateAcuteness_                   ; ///< acuteness of rotation steps of key presses
    tgt::KeyEvent::KeyCode keyRotateLeft_       ; ///< at which key code rotate left
    tgt::KeyEvent::KeyCode keyRotateRight_      ; ///< at which key code rotate right
    tgt::KeyEvent::KeyCode keyRotateUp_         ; ///< at which key code rotate up
    tgt::KeyEvent::KeyCode keyRotateDown_       ; ///< at which key code rotate down
    int keyRotateMod_                           ; ///< at which modifiers to rotate by keys
    bool keyRotatePressed_                      ;
    float keyMoveAcuteness_                     ;
    tgt::KeyEvent::KeyCode keyMoveLeft_         ;
    tgt::KeyEvent::KeyCode keyMoveRight_        ;
    tgt::KeyEvent::KeyCode keyMoveUp_           ;
    tgt::KeyEvent::KeyCode keyMoveDown_         ;
    int keyMoveMod_                             ;
    bool keyMovePressed_                        ;
    tgt::vec2 mouseZoomInDirection_             ;
    float mouseWheelZoomAcuteness_              ;
    bool mouseWheelUpZoomIn_                    ;
    float keyZoomAcuteness_                     ;
    tgt::KeyEvent::KeyCode keyZoomIn_           ;
    tgt::KeyEvent::KeyCode keyZoomOut_          ;
    int keyZoomMod_                             ;
    bool keyZoomPressed_                        ;
    float mouseRollAcuteness_                   ;
    float mouseWheelRollAcuteness_              ;
    bool mouseWheelUpRollLeft_                  ;
    float keyRollAcuteness_                     ;
    tgt::KeyEvent::KeyCode keyRollLeft_         ;
    tgt::KeyEvent::KeyCode keyRollRight_        ;
    int keyRollMod_                             ;
    bool keyRollPressed_                        ;
    tgt::MouseEvent::MouseButtons resetButton_  ; ///< If this button is double-clicked, reset trackball to initial
                                                  /// position and orientation

    static const std::string loggerCat_;
};

} // namespace

#endif //VRN_TRACKBALL_NAVIGATION_H
