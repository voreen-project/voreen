/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
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

#ifndef TGT_TRACKBALL_H
#define TGT_TRACKBALL_H

#include <cmath>
#include <vector>

#include "tgt/matrix.h"
#include "tgt/quaternion.h"
#include "tgt/vector.h"
#include "tgt/types.h"
#include "tgt/navigation/navigation.h"
#include "tgt/timer.h"
#include "tgt/stopwatch.h"


namespace tgt {

/**
    This class implements a Trackball which can be used to freely rotate an object
    around a given center (or, if you will, to rotate the camera around and focus it
    on that center while it is continuously moving on a sphere).
*/
class Trackball : public Navigation {
    public:

        /// Constructor
        ///
        /// Center of Trackball is set to camera focus of the given canvas.
        /// \param canvas The canvas the Trackball is associated with
        /// \param defaultEventHandling You can specify the mouse buttons and keys used to
        ///                             interact with trackball. If you want default behavior,
        ///                             set this to true and trackball should easyly work. If you
        ///                             want to specify different behavior than default, give false
        ///                             here. You will then have to call some setMouseFoo and
        ///                             setKeyBar methods to make the trackball work.
        /// \param continousSpinTimer To offer continous spin feature, trackball needs a
        ///                           tgt::Timer-object. As we need a toolkit-specific timer,
        ///                           trackball cannot create it itself, the user has to provide
        ///                           it to the constructor.
        Trackball(GLCanvas* canvas, Camera* camera, bool defaultEventHandling = true, Timer* continousSpinTimer = NULL);

        /// Destructor
        virtual ~Trackball();

        /// Resets the trackball to the initial configuration of the canvas' camera.
        void reset();

        /// Rotate the trackball according to Quaternion quat.
        /// @param quat Quaternion represention rotation.
        void rotate(Quaternion<float> quat);
        /// Rotate the trackball by angle phi around axis axis.
        /// @param axis axis in camera coordinates.
        virtual void rotate(vec3 axis, float phi);
        /// Rotate the trackball according to new mouse position.
        /// @param mouse coodinates of mouse, scaled to [-1, 1]x[-1, 1]
        void rotate(vec2 mouse);

        /// Move the trackball along axis.
        /// @param length Distance relative to frustum dimensions at trackball center.
        ///               Trackball center will get out of sight when >> 1.
        /// @param axis Axis in camera coordinates along which to move.
        virtual void move(float length, vec3 axis);
        /// Move the trackball according to new mouse position.
        /// @param mouse coodinates of mouse, scaled to [-1, 1]x[-1, 1]
        void move(vec2 mouse);

        /// Zoom in by factor factor.
        virtual void zoom(float factor);
        /// Zoom in according to new mouse position.
        /// @param mouse coodinates of mouse, scaled to [-1, 1]x[-1, 1]
        void zoom(vec2 mouse);
        /// set an absolute Distance from Focus
        void zoomAbsolute(float focallength);

        /// getter / setter

        /// whether this Trackball can be set to continious rotation by sweeping it
        void setContinuousSpin(bool b) { continuousSpin_ = b; }
        bool getContinuousSpin() const { return continuousSpin_; }

        /// last rotation applied to trackball
        quat getLastOrientationChange() const { return lastOrientationChange_; }

        /// If the center of the trackball should be moved relative to objects.
        /// This has large influence on trackball behavior. If you look at one object located at
        /// world coordinates origin, you might want the center not to move. This is used in voreen.
        /// If you use trackball to look at mutiple objects, located at varying position, you might
        /// want the center to move.
        /// See trackball demo in samples folder.
        void setMoveCenter(const bool& b = true) { moveCenter_ = b; };
        bool getMoveCenter() const { return moveCenter_; };

        /// Set trackballs center. Most people do not need this method. Use with caution!
        void setCenter(const vec3& center) { center_ = center; };
        const vec3& getCenter() const { return center_; };

        /// set the radius of the trackball to a given value
        void setSize(const float& size) { size_ = size; };
        float getSize() const { return size_; };

        /// Returns the orthogonal distance between camera position and center of the trackball
        /// with respect to the camera's look vector.
        float getCenterDistance();

        GLCanvas* getCanvas() const { return canvas_; }
        void setCanvas(GLCanvas* canvas) { canvas_ = canvas; }

        Timer* getContinuousSpinTimer() const { return continuousSpinTimer_; }

// Event Handling /////////////////////////////////////////////////////////////////////////////

        /// let trackball react on mouse movements when button is pressed and modifiers are suitable
        /// to mod with rotating
        ///
        /// use MouseEvent::NO_MOUSE_BUTTON as button to disable mouse rotation
        void setMouseRotate(MouseEvent::MouseButtons button = MouseEvent::MOUSE_BUTTON_LEFT,
                            int mod = Event::MODIFIER_NONE);
        /// let trackball react on key pressures with rotating
        /// @param acuteness The per-keypress angle of rotation will be smaller at greater acuteness. Use acuteness = 0.f to disable key rotation.
        /// @param left, right, up, down keycode which should cause suitable rotation
        /// @param mod which modifiers must be set
        /// @param pressed rotate on key-pressed-event when true, on key-release-event when false
        void setKeyRotate(float acuteness = 10.f,
                          KeyEvent::KeyCode left  = KeyEvent::K_LEFT,
                          KeyEvent::KeyCode right = KeyEvent::K_RIGHT,
                          KeyEvent::KeyCode up    = KeyEvent::K_UP,
                          KeyEvent::KeyCode down  = KeyEvent::K_DOWN,
                          int mod = Event::MODIFIER_NONE,
                          bool pressed = false);
        /// let trackball react on mouse movements when button is pressed and modifiers are suitable
        /// to mod with moving
        ///
        /// use MouseEvent::NO_MOUSE_BUTTON as button to disable mouse movement
        void setMouseMove(MouseEvent::MouseButtons button = MouseEvent::MOUSE_BUTTON_LEFT,
                          int mod = Event::SHIFT);
        /// let trackball react on key presses with moving
        /// @param acuteness The per-keypress length of movement will be smaller at greater
        ///        acuteness. Use acuteness = 0.f to disable key rotation.
        /// @param left, right, up, down keycode which should cause suitable movement
        /// @param mod which modifiers must be set
        /// @param pressed move on key-pressed-event when true, on key-release-event when false
        void setKeyMove(float acuteness = 100.f,
                        KeyEvent::KeyCode left  = KeyEvent::K_LEFT,
                        KeyEvent::KeyCode right = KeyEvent::K_RIGHT,
                        KeyEvent::KeyCode up    = KeyEvent::K_UP,
                        KeyEvent::KeyCode down  = KeyEvent::K_DOWN,
                        int mod = Event::SHIFT,
                        bool pressed = false);

        /// let trackball react on mouse movements when button is pressed and modifiers are suitable
        /// to mod with zooming (changing the trackball's size)
        ///
        /// use MouseEvent::NO_MOUSE_BUTTON as button to disable mouse zooming
        ///
        /// @param zoomInDirection specifies in which direction mouse must be moved to zoom in.
        ///        The greater the length of this vector, the longer the way the mouse must be moved
        ///        to achieve a certain zoom factor.
        void setMouseZoom(MouseEvent::MouseButtons button = MouseEvent::MOUSE_BUTTON_LEFT,
                          vec2 zoomInDirection = vec2( 0.f, 1.f),
                          int mod = Event::CTRL);
        /// let trackball react on mouse wheel events when modifiers are suitable to mod with
        /// zooming (changing the trackball's size)
        /// @param acuteness The zoom factor will be smaller at greater acuteness. Use acuteness = 0.f to disable key rotation.
        /// @param wheelUpZoomIn zoom in on wheel up and out on wheel down when true, otherwise when
        ///        false
        void setMouseWheelZoom(float acuteness = 10.f, bool wheelUpZoomIn = true,
                               int mod = Event::MODIFIER_NONE);
        /// @param acuteness The zoom factor will be smaller at greater acuteness. Use acuteness = 0.f to disable key rotation.
        void setKeyZoom(float acuteness = 10.f,
                        KeyEvent::KeyCode in  = KeyEvent::K_UP,
                        KeyEvent::KeyCode out = KeyEvent::K_DOWN,
                        int mod = Event::CTRL,
                        bool pressed = false);

        /// let trackball react on mouse movements when button is pressed and modifiers are suitable
        /// to mod with rolling (tiling camera to left or right)
        ///
        /// use MouseEvent::NO_MOUSE_BUTTON as button to disable mouse rolling
        ///
        /// @param acuteness greater the acuteness means less tiling
        void setMouseRoll(MouseEvent::MouseButtons button = MouseEvent::MOUSE_BUTTON_LEFT,
                          float acuteness = .5f,
                          int mod = Event::CTRL);
        /// let trackball react on mouse wheel events when modifiers are suitable to mod with
        /// rolling  (tiling camera to left or right)
        /// @param acuteness Less tiling at greater acuteness. Use acuteness = 0.f to disable key rotation.
        /// @param wheelUpRollLeft roll left on wheel up and right on wheel down when true,
        ///                        otherwise when false
        void setMouseWheelRoll(float acuteness = 10.f, bool wheelUpRollLeft = true,
                               int mod = Event::SHIFT);
        /// @param acuteness Less tiling at greater acuteness. Use acuteness = 0.f to disable key rotation.
        void setKeyRoll(float acuteness = 10.f,
                        KeyEvent::KeyCode left = KeyEvent::K_LEFT,
                        KeyEvent::KeyCode right = KeyEvent::K_RIGHT,
                        int mod = Event::ALT,
                        bool pressed = false);

       MouseEvent::MouseButtons getRotateButton() const { return mouseRotateButton_; }
       MouseEvent::MouseButtons getZoomButton() const { return mouseZoomButton_; }
       MouseEvent::MouseButtons getMoveButton() const { return mouseMoveButton_; }
       MouseEvent::MouseButtons getRollButton() const { return mouseRollButton_; }

        /// implement EventListener functions
        virtual void mousePressEvent(MouseEvent* e);
        virtual void mouseReleaseEvent(MouseEvent* e);
        virtual void mouseMoveEvent(MouseEvent* e);
        virtual void wheelEvent(MouseEvent* e);
        virtual void keyEvent(KeyEvent* e);
        virtual void timerEvent(TimeEvent* e);

        /**
         * Set the initial values for the camera. This is required, whenever a
         * camera / trackball setup has been saved (e.g. to disk as XML-file)
         * and needs to be reloaded.
         */
        void reinitializeCamera(const vec3& position, const vec3& focus, const vec3& upVector) {
            cameraPosition_ = position;
            cameraFocus_ = focus;
            cameraUpVector_ = upVector;
            if (camera_ != 0)
                reset();
        }

///////////////////////////////////////////////////////////////////////////////////////////////

    protected:
        //the canvas
        GLCanvas* canvas_;
        /// is this a continuously spinning Trackball?
        bool continuousSpin_;
        /// timer used to do continuous spinning
        Timer* continuousSpinTimer_;
        /// stopwatch to get the duration between mouse events
        Stopwatch* continuousSpinStopwatch_;
        /// msecs between the last two mouse motion events (the time the last orientation change took)
        clock_t continuousSpinLastOrientationChangeMSecs_;

        /// Old mouse coordinates to allow calculation of the relative change.
        /// Ought to be relative coordinates within range [-1, 1]x[-1, 1].
        vec2 mouse_;

        /// this holds the center around which the camera will be rotated
        vec3 center_;
        /// Wheather to change the center in move method.
        /// center_ will alway be (0,0,0) if false.
        bool moveCenter_;
        /// the trackball size effects the trackball sensitivity
        float size_;

        /// store the initial position, focus point and up vector of camera used to be able to
        /// reset camera to those values
        void saveCameraParameters();
        vec3 cameraPosition_;
        vec3 cameraFocus_;
        vec3 cameraUpVector_;

        /// last rotation applied to trackball
        quat lastOrientationChange_;

        /// Stores mouse coordinates to be able to track relative mouse motions.
        /// Should be called when mouse buttons get pressed.
        void startMouseDrag(MouseEvent* e);

        /// Should be called when mouse buttons get released.
        void endMouseDrag(MouseEvent* e);

        /// scale screen-coodinates of mouse to intervall [-1, 1]x[-1, 1]
        vec2 scaleMouse(const ivec2& sc) const;

        /// projection math
        vec3 projectToSphere(const vec2 xy) const;

        /// transform vectors given in camera coordinates to vectors in world coordinates
        vec3 coordTransform(vec3 const axis) const;

        /// prepare trackball for usage of continuous spin feature
        void initializeContinuousSpin();

// Event Handling /////////////////////////////////////////////////////////////////////////////
        /// used to store settings how to react on certain events
        MouseEvent::MouseButtons mouseRotateButton_ ; ///< at which mouse button to rotate
        int mouseRotateMod_                         ; ///< at which modifiers to rotate by mouse
        float keyRotateAcuteness_                   ; ///< acuteness of rotation steps of key presses
        KeyEvent::KeyCode keyRotateLeft_            ; ///< at which key code rotate left
        KeyEvent::KeyCode keyRotateRight_           ; ///< at which key code rotate right
        KeyEvent::KeyCode keyRotateUp_              ; ///< at which key code rotate up
        KeyEvent::KeyCode keyRotateDown_            ; ///< at which key code rotate down
        int keyRotateMod_                           ; ///< at which modifiers to rotate by keys
        bool keyRotatePressed_                      ; // FIXME: cannot remember. not in use anyway...
        MouseEvent::MouseButtons mouseMoveButton_   ;
        int mouseMoveMod_                           ;
        float keyMoveAcuteness_                     ;
        KeyEvent::KeyCode keyMoveLeft_              ;
        KeyEvent::KeyCode keyMoveRight_             ;
        KeyEvent::KeyCode keyMoveUp_                ;
        KeyEvent::KeyCode keyMoveDown_              ;
        int keyMoveMod_                             ;
        bool keyMovePressed_                        ; // FIXME: cannot remember, not in use
        MouseEvent::MouseButtons mouseZoomButton_   ;
        vec2 mouseZoomInDirection_                  ;
        int mouseZoomMod_                           ;
        bool mouseWheelZoom_                        ;
        float mouseWheelZoomAcuteness_              ;
        bool mouseWheelUpZoomIn_                    ;
        int mouseWheelZoomMod_                      ;
        float keyZoomAcuteness_                     ;
        KeyEvent::KeyCode keyZoomIn_                ;
        KeyEvent::KeyCode keyZoomOut_               ;
        int keyZoomMod_                             ;
        bool keyZoomPressed_                        ; // FIXME: cannot remember, not in use
        MouseEvent::MouseButtons mouseRollButton_   ;
        float mouseRollAcuteness_                   ;
        int mouseRollMod_                           ;
        bool mouseWheelRoll_                        ;
        float mouseWheelRollAcuteness_              ;
        bool mouseWheelUpRollLeft_                  ;
        int mouseWheelRollMod_                      ;
        float keyRollAcuteness_                     ;
        KeyEvent::KeyCode keyRollLeft_              ;
        KeyEvent::KeyCode keyRollRight_             ;
        int keyRollMod_                             ;
        bool keyRollPressed_                        ; // FIXME: cannot remember, not in use

        bool tracking_; ///< Are we tracking mouse move events? Only when we received a
                        ///< mousePressEvent before.


        void initializeEventHandling() {
            mouseRotateButton_ = mouseMoveButton_ = mouseZoomButton_ = mouseRollButton_
                = MouseEvent::MOUSE_BUTTON_NONE;
            keyRotateLeft_ = keyRotateRight_ = keyRotateUp_ = keyRotateDown_
                = keyMoveLeft_ = keyMoveRight_ = keyMoveUp_ = keyMoveDown_
                = keyZoomIn_ = keyZoomOut_
                = keyRollLeft_ = keyRollRight_
                = KeyEvent::K_LAST;
            mouseWheelZoom_ = mouseWheelRoll_ = false;
        }

        float getRotationAngle(const float& acuteness) const ;
        float getMovementLength(const float& acuteness) const ;
        float getZoomFactor(const float& acuteness, const bool& zoomIn) const ;
        float getRollAngle(const float& acuteness, const bool& left) const ;


///////////////////////////////////////////////////////////////////////////////////////////////

};

} // namespace tgt

#endif // TGT_TRACKBALL_H
