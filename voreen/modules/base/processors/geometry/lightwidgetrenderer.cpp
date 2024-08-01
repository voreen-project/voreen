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

#include "lightwidgetrenderer.h"
#include "voreen/core/utils/voreenqualitymode.h"

namespace voreen {

using tgt::vec4;
using tgt::vec3;

LightWidgetRenderer::LightWidgetRenderer()
    : GeometryRendererBase()
    , showLightWidget_("set.showLightWidget", "Show Light Widget", true)
    , isClicked_(false)
    , sphereRadius_("sphereRadius", "Radius", 0.03f, 0.0f, 10000.0f)
    , sphereColor_("sphereColor", "Color", tgt::vec4(0.24f,0.2f,0.07f,1), Processor::INVALID_RESULT, Property::LOD_DEFAULT, false)
    , lightPosition_("lightPosition", "Light Source Position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.f),
                     tgt::vec4(-10000.f), tgt::vec4(10000.f))
    , mesh_()
    , lightSource_()
    , material_()
{

    moveSphereProp_ = new EventProperty<LightWidgetRenderer>(
        "mouseEvent.moveSphere", "Light widget motion",
        this, &LightWidgetRenderer::moveSphere,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::MOTION | tgt::MouseEvent::RELEASED,
        tgt::Event::MODIFIER_NONE);
    addProperty(showLightWidget_);
    addProperty(lightPosition_);
    addProperty(sphereColor_);
    addProperty(sphereRadius_);
        sphereRadius_.setNumDecimals(3);
    addEventProperty(moveSphereProp_);

    mesh_.setSphereGeometry(1, tgt::vec3::zero, tgt::vec4::one, 30);

    //lightPosition_.setViews(Property::View(Property::LIGHT_POSITION | Property::DEFAULT));

    // light parameters
    lightSource_.position = tgt::vec4(0,1,1,0);
    lightSource_.ambientColor = tgt::vec3(1,1,1);
    lightSource_.diffuseColor = tgt::vec3(1,1,1);
    lightSource_.specularColor = tgt::vec3(1,1,1);

    // parameters for yellow plastic
    material_.shininess = 51.0f;
}

LightWidgetRenderer::~LightWidgetRenderer() {
    delete moveSphereProp_;
}

Processor* LightWidgetRenderer::create() const {
    return new LightWidgetRenderer();
}

void LightWidgetRenderer::moveSphere(tgt::MouseEvent* e) {
    LGL_ERROR;
    if (!idManager_)
        return;

    if (e->action() & tgt::MouseEvent::PRESSED) {
        if (idManager_->isHit(tgt::ivec2(e->x(), e->viewport().y - e->y() ), this)) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            e->accept();
            invalidate();
            isClicked_ = true;
            lightPositionAbs_.x = lightPosition_.get().x;
            lightPositionAbs_.y = lightPosition_.get().y;
            lightPositionAbs_.z = lightPosition_.get().z;
            startCoord_.x = e->coord().x;
            startCoord_.y = e->coord().y;
        }
        return;
    }

    if (e->action() & tgt::MouseEvent::MOTION) {
        if (isClicked_) {
            e->accept();

            GLint deltaX, deltaY;

            GLint viewport[4];
            GLdouble winX, winY, winZ;
            GLdouble posX, posY, posZ;

            deltaX = e->coord().x - startCoord_.x;
            deltaY = startCoord_.y - e->coord().y;

            tgt::dmat4 projection_transposed = tgt::transpose(camera_.getProjectionMatrix(idManager_->getRenderTarget()->getSize()));
            tgt::dmat4 modelview_transposed = tgt::transpose(camera_.getViewMatrix());
            viewport[0] = 0;
            viewport[1] = 0;
            viewport[2] = static_cast<GLint>(idManager_->getRenderTarget()->getSize().x);
            viewport[3] = static_cast<GLint>(idManager_->getRenderTarget()->getSize().y);

            posX = lightPositionAbs_.x;
            posY = lightPositionAbs_.y;
            posZ = lightPositionAbs_.z;

            tgt::gluProject(posX, posY, posZ, modelview_transposed.elem, projection_transposed.elem, viewport, &winX, &winY, &winZ);

            winX = winX + deltaX;
            winY = winY + deltaY;

            tgt::gluUnProject(winX, winY, winZ, modelview_transposed.elem, projection_transposed.elem, viewport, &posX, &posY, &posZ);

            lightPosition_.set(vec4(static_cast<float>(posX), static_cast<float>(posY), static_cast<float>(posZ), 1.f));

            invalidate();
        }
        return;
    }

    if (e->action() & tgt::MouseEvent::RELEASED) {
        if (isClicked_) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            e->accept();
            isClicked_ = false;
        }
        return;
    }
}

void LightWidgetRenderer::render() {
    if (showLightWidget_.get()) {
        IMode.setLightSource(lightSource_);
        IMode.color(sphereColor_.get());
        IMode.setMaterial(material_);

        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.translate(lightPosition_.get().xyz());
        MatStack.scale(tgt::vec3(sphereRadius_.get()));

        mesh_.render();
        LGL_ERROR;

        MatStack.popMatrix();
    }
}

void LightWidgetRenderer::renderPicking() {
    if (!idManager_)
        return;
    if (showLightWidget_.get()) {
        LGL_ERROR;

        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.translate(lightPosition_.get().xyz());
        MatStack.scale(tgt::vec3(sphereRadius_.get()));

        idManager_->setGLColor(this);
        GlMeshGeometryUInt16Simple::renderDefault(mesh_);
        IMode.color(tgt::vec4::one);
        LGL_ERROR;

        MatStack.popMatrix();
    }
}

tgt::Bounds LightWidgetRenderer::getBoundingBox() const {
    return mesh_.getBoundingBox();
}

void LightWidgetRenderer::setIDManager(IDManager* idm) {
    if (idManager_ == idm)
        return;

    idManager_ = idm;
    if (idManager_) {
        idm->registerObject(this);
    }
}

}

