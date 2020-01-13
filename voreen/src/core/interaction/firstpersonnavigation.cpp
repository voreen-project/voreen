/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/interaction/firstpersonnavigation.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/cameraproperty.h"

#include "tgt/camera.h"
#include "tgt/event/keyevent.h"
#include "tgt/event/mouseevent.h"
#include "tgt/stopwatch.h"
#include "tgt/timer.h"

namespace voreen{

    FirstPersonNavigation::FirstPersonNavigation(CameraProperty* camera)
            : camera_(camera)
            , movingDirection_(tgt::vec3::zero)
            , firstMove_(true)
    {
        initialize();
    }

    FirstPersonNavigation::~FirstPersonNavigation(){
        delete stopwatch_;
    }

    void FirstPersonNavigation::initialize(){
        stopwatch_ = new tgt::Stopwatch;

        loadDefaultKeysetting();
        loadDefaultMousesetting();
        loadDefaultCamerasetting();
    }

    void FirstPersonNavigation::keyEvent(tgt::KeyEvent* keyEve){
        if (keyEve->pressed()){
            if (keyMoveForward_ == keyEve->keyCode()){
                firstMove_ = !isMoving();
                movingDirection_.z = 1;
            }
            else if (keyMoveBackward_ == keyEve->keyCode()){
                firstMove_ = !isMoving();
                movingDirection_.z = -1;
            }
            else if (keyMoveLeft_ == keyEve->keyCode()){
                firstMove_ = !isMoving();
                movingDirection_.x = 1;
            }
            else if (keyMoveRight_ == keyEve->keyCode()){
                firstMove_ = !isMoving();
                movingDirection_.x = -1;
            }
            else if (keyMoveDown_ == keyEve->keyCode()){
                firstMove_ = !isMoving();
                movingDirection_.y = 1;
            }
            else if (keyMoveUp_ == keyEve->keyCode()){
                firstMove_ = !isMoving();
                movingDirection_.y = -1;
            }

            move();
        }
        else {
            if (keyMoveForward_ == keyEve->keyCode() || keyMoveBackward_ == keyEve->keyCode())
                movingDirection_.z = 0;
            else if (keyMoveRight_ == keyEve->keyCode() || keyMoveLeft_ == keyEve->keyCode())
                movingDirection_.x = 0;
            else if (keyMoveUp_ == keyEve->keyCode() || keyMoveDown_ == keyEve->keyCode())
                movingDirection_.y = 0;

            if (!isMoving()){
                curAccuteness_ = minAccuteness_;
                stopwatch_->stop();
                firstMove_ = true;
            }
        }
    }

    void FirstPersonNavigation::mousePressEvent(tgt::MouseEvent* mouseEve){
        lastMousePosition_ = tgt::ivec2(mouseEve->x(), mouseEve->y());
        mouseEve->accept();
    }

    void FirstPersonNavigation::mouseMoveEvent(tgt::MouseEvent* mouseEve){
        tgt::ivec2 newMouse = tgt::ivec2(mouseEve->x(), mouseEve->y());

        if ((mouseEve->button() == buttonRotateUpAndStrafe_)){
            if (lastMousePosition_.x != newMouse.x)
                rollCameraAroundUp(-computeRotationAngle(lastMousePosition_.x, newMouse.x));
            if (lastMousePosition_.y != newMouse.y)
                rollCameraAroundStrafe(-computeRotationAngle(lastMousePosition_.y, newMouse.y));
        }
        else if (mouseEve->button() == buttonRotateLook_)
            if (lastMousePosition_.x != newMouse.x)
                rollCameraAroundLook(computeRotationAngle(lastMousePosition_.x, newMouse.x));

        lastMousePosition_ = newMouse;

        // consider movement for simultaneous use of mouse and keyboard
        move();

        mouseEve->accept();
    }

    void FirstPersonNavigation::timerEvent(tgt::TimeEvent* /* e */){
        move();
    }

    void FirstPersonNavigation::loadDefaultKeysetting(){
        setKeysetting(tgt::KeyEvent::K_W,
                      tgt::KeyEvent::K_S,
                      tgt::KeyEvent::K_A,
                      tgt::KeyEvent::K_D,
                      tgt::KeyEvent::K_R,
                      tgt::KeyEvent::K_F);
    }

    void FirstPersonNavigation::loadDefaultMousesetting(){
        setMousesetting(tgt::MouseEvent::MOUSE_BUTTON_LEFT,
                        tgt::MouseEvent::MOUSE_BUTTON_RIGHT);
    }


    void FirstPersonNavigation::loadDefaultCamerasetting(){
        setCameraSetting(1, 40, 0.0005f, 0.15f, 0.005f, 0.4f);
        curAccuteness_ = minAccuteness_;
    }

    void FirstPersonNavigation::setKeysetting(tgt::KeyEvent::KeyCode keyMoveForward,
                                              tgt::KeyEvent::KeyCode keyMoveBackward,
                                              tgt::KeyEvent::KeyCode keyMoveLeft,
                                              tgt::KeyEvent::KeyCode keyMoveRight,
                                              tgt::KeyEvent::KeyCode keyMoveUp,
                                              tgt::KeyEvent::KeyCode keyMoveDown){
        keyMoveForward_  = keyMoveForward;
        keyMoveBackward_ = keyMoveBackward;
        keyMoveLeft_     = keyMoveLeft;
        keyMoveRight_    = keyMoveRight;
        keyMoveUp_       = keyMoveUp;
        keyMoveDown_     = keyMoveDown;
    }

    void FirstPersonNavigation::setMousesetting(tgt::MouseEvent::MouseButtons buttonRotateUpAndStrafe,
                                                tgt::MouseEvent::MouseButtons buttonRotateLook){
        buttonRotateUpAndStrafe_ = buttonRotateUpAndStrafe;
        buttonRotateLook_        = buttonRotateLook;
    }

    void FirstPersonNavigation::setCameraSetting(int minMouseSteps,
                                                 int maxMouseSteps,
                                                 float minRotationAngle,
                                                 float maxRotationAngle,
                                                 float minAccuteness,
                                                 float maxAccuteness){
        minMouseSteps_    = minMouseSteps;
        maxMouseSteps_    = maxMouseSteps;
        minRotationAngle_ = minRotationAngle;
        maxRotationAngle_ = maxRotationAngle;
        minAccuteness_    = minAccuteness;
        maxAccuteness_    = maxAccuteness;
    }

    float FirstPersonNavigation::computeRotationAngle(int oldMouse, int newMouse){
        int distance = std::abs(newMouse - oldMouse);

        float rotationAngle;

        if (distance == 0)
            rotationAngle = 0;
        else if (distance >= maxMouseSteps_)
            rotationAngle = maxRotationAngle_;
        else if (distance <= minMouseSteps_)
            rotationAngle = minRotationAngle_;
        else{
            float convexFactor = static_cast<float>(distance) / (maxMouseSteps_ - minMouseSteps_);
            rotationAngle = minRotationAngle_ + convexFactor * (maxRotationAngle_ - minRotationAngle_);
        }

        return newMouse > oldMouse ? rotationAngle : -rotationAngle;
    }

    void FirstPersonNavigation::move(){
        move(curAccuteness_, movingDirection_);
    }

    void FirstPersonNavigation::move(float length, tgt::vec3 axis) {
        if (length == 0.f || tgt::length(axis) == 0.f)
            return;

        float frustFactor = camera_->get().getFocalLength() / camera_->get().getFrustum().getNearDist();
        float frustWidth  = ( camera_->get().getFrustum().getRight()
                            - camera_->get().getFrustum().getLeft() );
        float frustHeight = ( camera_->get().getFrustum().getTop()
                              - camera_->get().getFrustum().getBottom() );

        float timestep;

        if (firstMove_)
            timestep = 1;
        else
            timestep = static_cast<float>(stopwatch_->getRuntime())/40;

        axis = normalize(axis) * length * frustFactor * timestep;
        axis.x *= frustWidth;
        axis.y *= frustHeight;
        axis.z *= (frustWidth+frustHeight)/2;

        tgt::mat4 rotation;
        camera_->get().getRotateMatrix().invert(rotation);

        tgt::Camera cam = camera_->get();
        cam.setPosition(cam.getPosition() - rotation*axis);
        cam.setFocus(cam.getFocus() - rotation*axis);
        camera_->set(cam);

        if (curAccuteness_ < maxAccuteness_)
            curAccuteness_ += 0.005f*timestep;

        firstMove_ = false;

        stopwatch_->reset();
        stopwatch_->start();
    }

    void FirstPersonNavigation::rollCameraAroundUp(float angle){
        float currentFocusDist = tgt::length(camera_->get().getPosition() - camera_->get().getFocus());
        tgt::vec3 look=normalize(tgt::quat::rotate(camera_->get().getLook(), angle, camera_->get().getUpVector()));
        tgt::Camera cam = camera_->get();
        cam.setFocus(cam.getPosition() + look*currentFocusDist);
        camera_->set(cam);
    }

    void FirstPersonNavigation::rollCameraAroundLook(float angle){
        tgt::vec3 up = normalize(tgt::quat::rotate(camera_->get().getUpVector(), angle, camera_->get().getLook()) );
        tgt::Camera cam = camera_->get();
        cam.setUpVector(up);
        camera_->set(cam);
    }

    void FirstPersonNavigation::rollCameraAroundStrafe(float angle){
        float currentFocusDist = tgt::length(camera_->get().getPosition() - camera_->get().getFocus());
        tgt::vec3 up = normalize(tgt::quat::rotate(camera_->get().getUpVector(), angle, camera_->get().getStrafe()) );
        tgt::Camera cam = camera_->get();
        cam.setUpVector(up);
        tgt::vec3 look = cross(up, camera_->get().getStrafe());
        cam.setFocus(cam.getPosition() + look*currentFocusDist);
        camera_->set(cam);
    }

    CameraProperty* FirstPersonNavigation::getCamera() const{
        return camera_;
    }

    void FirstPersonNavigation::setCamera(CameraProperty* camera){
        tgtAssert(camera, "No camera");
        camera_ = camera;
    }

    bool FirstPersonNavigation::isMoving(){
        return movingDirection_ != tgt::vec3::zero;
    }
}
