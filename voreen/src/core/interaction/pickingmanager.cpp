/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/interaction/pickingmanager.h"
#include "voreen/core/interaction/idmanager.h"

#include "tgt/glmath.h"
#include "tgt/material.h"
#include "tgt/camera.h"

using tgt::vec3;

namespace voreen {

PickingManager::PickingManager()
    : idManager_(NULL)
    , initialized_(false)
{}

PickingManager::~PickingManager(){
    delete idManager_;

    for (std::vector<PickableObject*>::iterator it = pickableObjects_.begin(); it != pickableObjects_.end(); ++it)
        delete (*it);
}

void PickingManager::setRenderTarget(RenderTarget* renderTarget){
    if (renderTarget != NULL){
        idManager_ = new IDManager();

        idManager_->setRenderTarget(renderTarget);
        idManager_->initializeTarget();

        initialized_ = true;
    }
}

void PickingManager::prepareDrawing(){
    idManager_->clearRegisteredObjects();
    idManager_->activateTarget();
    idManager_->clearTarget();

    std::vector<PickableObject*> tempPicked;

    for (std::vector<PickableObject*>::iterator it = pickableObjects_.begin(); it != pickableObjects_.end(); ++it)
        if ((*it)->isClicked())
            tempPicked.push_back((*it));
        else
            delete (*it);

    pickableObjects_.clear();
    nonDynamicPickableObjects_.clear();
    pickableObjects_ = tempPicked;
}

void PickingManager::addPickableObject(PickableObject* pickableObject){
    if (pickableObject != NULL){
        pickableObjects_.push_back(pickableObject);

        idManager_->registerObject(pickableObject);

        pickableObject->setRenderTargetDimensions(idManager_->getRenderTarget()->getSize());
    }
}

void PickingManager::addNonDynamicPickableObject(PickableObject* pickableObject){
    if (pickableObject != NULL){
        nonDynamicPickableObjects_.push_back(pickableObject);

        idManager_->registerObject(pickableObject);

        pickableObject->setRenderTargetDimensions(idManager_->getRenderTarget()->getSize());
    }
}

bool PickingManager::setGLColor(PickableObject* pickableObject){
    if (idManager_->isRegistered(pickableObject)){
        idManager_->setGLColor(pickableObject);

        return true;
    }

    return false;
}

bool PickingManager::onMouseEvent(tgt::MouseEvent* mouseEve, tgt::Camera* sceneCamera){
    tgt::ivec2 position(mouseEve->x(), mouseEve->viewport().y - mouseEve->y());

    bool changes = false;

    const void* pickedOne = idManager_->getObjectAtPos(position);

    switch (mouseEve->action()){
        case tgt::MouseEvent::MOTION:
            for (std::vector<PickableObject*>::iterator it = pickableObjects_.begin(); it != pickableObjects_.end(); ++it){
                if ((*it)->isClicked()){
                    (*it)->mouseMove(mouseEve, sceneCamera);
                    changes = true;
                }
            }

            for (std::vector<PickableObject*>::iterator it = nonDynamicPickableObjects_.begin(); it != nonDynamicPickableObjects_.end(); ++it){
                if ((*it)->isClicked()){
                    (*it)->mouseMove(mouseEve, sceneCamera);
                    changes = true;
                }
            }

            break;

        case tgt::MouseEvent::RELEASED:
             for (std::vector<PickableObject*>::iterator it = pickableObjects_.begin(); it != pickableObjects_.end(); ++it){
                if ((*it)->isClicked()){
                    (*it)->mouseRelease(mouseEve, sceneCamera);
                    changes = true;
                }
             }


            for (std::vector<PickableObject*>::iterator it = nonDynamicPickableObjects_.begin(); it != nonDynamicPickableObjects_.end(); ++it){
                if ((*it)->isClicked()){
                    (*it)->mouseRelease(mouseEve, sceneCamera);
                    changes = true;
                }
            }
            break;

        case tgt::MouseEvent::PRESSED:
        {
             pickedOne = idManager_->getObjectAtPos(position);

             if (pickedOne != 0){
                 PickableObject* pickedObject = const_cast<PickableObject*>(reinterpret_cast<const PickableObject*>(pickedOne));

                 pickedObject->mousePress(mouseEve, sceneCamera);
                 changes = mouseEve->isAccepted();
             }
             break;
        }

        default:
            break;
    }

    return changes;
}

bool PickingManager::isInitialized() const{
    return initialized_;
}

// ================================================================================================

PickableObject::PickableObject()
        : isClicked_(false)
{}

PickableObject::~PickableObject()
{}

void PickableObject::mousePress(tgt::MouseEvent* mouseEve, tgt::Camera*  sceneCamera){
    isClicked_ = true;

    pressedButton_ = mouseEve->button();

    oldScreenPosition_ = tgt::ivec2(mouseEve->x(), mouseEve->viewport().y - mouseEve->y());

    newScreenPosition_ = tgt::ivec2(mouseEve->x(), mouseEve->viewport().y - mouseEve->y());

    //GLint deltaX, deltaY;

    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLdouble posX, posY, posZ;

    //deltaX = newScreenPosition_.x - oldScreenPosition_.x;
    //deltaY = newScreenPosition_.y - oldScreenPosition_.y;

    tgt::mat4 projection_tgt = sceneCamera->getProjectionMatrix(mouseEve->viewport());
    tgt::mat4 modelview_tgt = sceneCamera->getViewMatrix();

    for (int i = 0; i < 4; ++i) {
        modelview[i+0]   = modelview_tgt[i].x;
        modelview[i+4]   = modelview_tgt[i].y;
        modelview[i+8]   = modelview_tgt[i].z;
        modelview[i+12]  = modelview_tgt[i].w;
        projection[i+0]  = projection_tgt[i].x;
        projection[i+4]  = projection_tgt[i].y;
        projection[i+8]  = projection_tgt[i].z;
        projection[i+12] = projection_tgt[i].w;
    }

    viewport[0] = 0;
    viewport[1] = 0;
    viewport[2] = static_cast<GLint>(renderTargetDimensions_.x);
    viewport[3] = static_cast<GLint>(renderTargetDimensions_.y);

    GLfloat WindowPosZ;

    glReadPixels(oldScreenPosition_.x, oldScreenPosition_.y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &WindowPosZ);

    LGL_ERROR;
    tgt::gluUnProject(oldScreenPosition_.x, oldScreenPosition_.y, WindowPosZ, modelview, projection, viewport, &posX, &posY, &posZ);
    LGL_ERROR;

    oldWorldPosition_ = tgt::vec3(static_cast<float>(posX), static_cast<float>(posY), static_cast<float>(posZ));
    newWorldPosition_ = tgt::vec3(static_cast<float>(posX), static_cast<float>(posY), static_cast<float>(posZ));

    mouseEve->accept();
}

void PickableObject::mouseRelease(tgt::MouseEvent*  mouseEve, tgt::Camera* /* sceneCamera */){
    isClicked_ = false;
    mouseEve->accept();
}

void PickableObject::mouseMove(tgt::MouseEvent* mouseEve, tgt::Camera* sceneCamera){
    newScreenPosition_ = tgt::ivec2(mouseEve->x(), mouseEve->viewport().y - mouseEve->y());

    GLint deltaX, deltaY;

    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLdouble winX, winY, winZ;
    GLdouble posX, posY, posZ;

    deltaX = newScreenPosition_.x - oldScreenPosition_.x;
    deltaY = newScreenPosition_.y - oldScreenPosition_.y;

    tgt::mat4 projection_tgt = sceneCamera->getProjectionMatrix(mouseEve->viewport());
    tgt::mat4 modelview_tgt = sceneCamera->getViewMatrix();

    for (int i = 0; i < 4; ++i) {
        modelview[i+0]   = modelview_tgt[i].x;
        modelview[i+4]   = modelview_tgt[i].y;
        modelview[i+8]   = modelview_tgt[i].z;
        modelview[i+12]  = modelview_tgt[i].w;
        projection[i+0]  = projection_tgt[i].x;
        projection[i+4]  = projection_tgt[i].y;
        projection[i+8]  = projection_tgt[i].z;
        projection[i+12] = projection_tgt[i].w;
    }

    viewport[0] = 0;
    viewport[1] = 0;
    viewport[2] = static_cast<GLint>(renderTargetDimensions_.x);
    viewport[3] = static_cast<GLint>(renderTargetDimensions_.y);

    posX = oldWorldPosition_.x;
    posY = oldWorldPosition_.y;
    posZ = oldWorldPosition_.z;

    LGL_ERROR;
    tgt::gluProject(posX, posY, posZ, modelview, projection, viewport, &winX, &winY, &winZ);

    winX = winX + deltaX;
    winY = winY + deltaY;

    LGL_ERROR;
    tgt::gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
    LGL_ERROR;

    newWorldPosition_ = tgt::vec3(static_cast<float>(posX), static_cast<float>(posY), static_cast<float>(posZ));
    mouseEve->accept();
}

void PickableObject::setRenderTargetDimensions(tgt::ivec2 renderTargetDimensions){
    renderTargetDimensions_ = renderTargetDimensions;
}

bool PickableObject::isClicked() const{
    return isClicked_;
}

void PickableObject::refreshPositions(){
    oldScreenPosition_ = newScreenPosition_;
    oldWorldPosition_ = newWorldPosition_;
}

}
