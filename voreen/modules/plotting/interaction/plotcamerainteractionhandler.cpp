/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "plotcamerainteractionhandler.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "tgt/event/mouseevent.h"
#include "tgt/camera.h"
#include "tgt/matrix.h"
#include "voreen/core/properties/eventproperty.h"

namespace voreen {

PlotCameraInteractionHandler::PlotCameraInteractionHandler() :
    InteractionHandler("dummy", "dummy")
{}

PlotCameraInteractionHandler::PlotCameraInteractionHandler(const std::string& id, const std::string& guiName,
                                                   CameraProperty* cameraProp, bool sharing, bool enabled)
    : InteractionHandler(id, guiName)
    , cameraProp_(cameraProp)
{
    tgtAssert(cameraProp, "No camera property");
    cameraProp_ = cameraProp;

    // event property
    rotateEvent_ = new EventProperty<PlotCameraInteractionHandler>(id + ".rotate", guiName + " Rotate", this,
        &PlotCameraInteractionHandler::rotateEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(rotateEvent_);
}

void PlotCameraInteractionHandler::rotateEvent(tgt::MouseEvent* e) {
        tgtAssert(cameraProp_, "No camera property");

        if (e->action() == tgt::MouseEvent::PRESSED) {
            lastMousePosition_ = tgt::vec2(
                 static_cast<float>(e->x())/static_cast<float>(e->viewport().x)
                ,static_cast<float>(e->y())/static_cast<float>(e->viewport().y));
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
        }
        else if (e->action() == tgt::MouseEvent::RELEASED) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
        }
        else if (e->action() == tgt::MouseEvent::MOTION) {
            tgt::vec2 newMouse = tgt::vec2(
                 static_cast<float>(e->x())/static_cast<float>(e->viewport().x)
                ,static_cast<float>(e->y())/static_cast<float>(e->viewport().y));
            rotate(120.f*(lastMousePosition_.x-newMouse.x), 120.f*(newMouse.y-lastMousePosition_.y));
            lastMousePosition_ = newMouse;
            e->accept();
        }

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (e->isAccepted()) {
        cameraProp_->invalidate();
    }
}

void PlotCameraInteractionHandler::rotate(float horAngle, float vertAngle) {
    //convert to radian
    vertAngle = tgt::deg2rad(vertAngle);
    horAngle = tgt::deg2rad(horAngle);

    //first we get the old camera position
    tgt::vec3 pos = cameraProp_->get().getPosition();

    // only rotate vertically when camera z position is in [-0.9, 0.9] to avoid flipping phenomenon
    double z = pos.z;
    if ((z < -0.95 || vertAngle > 0) && (z > 0.95 || vertAngle < 0))
        vertAngle = 0;

    //read backwards! this is the second rotation, vertical from the point of the camera
    pos = tgt::Matrix::createRotation(vertAngle, tgt::cross(pos, tgt::vec3(0,0,1))) * pos;

    //this is the first rotation, horizontal from the point of the camera (and the
    //world coordinate system)
    pos = tgt::Matrix::createRotationZ(horAngle) * pos;

    tgt::Camera cam = cameraProp_->get();
    cam.setPosition(pos);
    cam.setFocus(tgt::vec3(0,0,0));
    cam.setUpVector(tgt::vec3(0,0,1));
    //cameraProp_->setPosition(pos);
    cameraProp_->set(cam);
}

void PlotCameraInteractionHandler::onEvent(tgt::Event* eve) {
    tgtAssert(cameraProp_, "No camera property");
    // invalidate processor and update camera prop widgets, if event has been accepted
    if (eve->isAccepted()) {
        cameraProp_->invalidate();
    }
}

} // namespace
