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

#include "distancemeasure.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/textureunit.h"

#include <sstream>

using tgt::TextureUnit;

namespace voreen {


DistanceMeasure::DistanceMeasure()
    : ImageProcessor("image/distance")
    , imgInport_(Port::INPORT, "image", "Image Input")
    , fhpInport_(Port::INPORT, "fhp", "First-hit-points Input", false, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F)
    , refInport_(Port::INPORT, "refvol", "Reference Volume", false)
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , mouseEventProp_("mouseEvent.measure", "Distance measure", this, &DistanceMeasure::measure, tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED, tgt::Event::CTRL, false)
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , renderSpheres_("renderSpheres", "Render Spheres", true)
    , font_(VoreenApplication::app()->getFontPath("VeraMono.ttf"), 16)
    , mouseCurPos2D_(0.0f)
    , mouseCurPos3D_(0.0f)
    , mouseStartPos2D_(0.0f)
    , mouseStartPos3D_(0.0f)
    , mouseDown_(false)
    , distance_(0)
    , mesh_()           //
    , lightSource_()    // Are initialized below
    , material_()       //
{
    addPort(imgInport_);
    addPort(fhpInport_);
    addPort(refInport_);
    addPort(outport_);

    addProperty(camera_);
    addProperty(renderSpheres_);

    addEventProperty(&mouseEventProp_);

    // light parameters
    lightSource_.position = tgt::vec4(0,1,1,0);
    lightSource_.ambientColor = tgt::vec3(1,1,1);
    lightSource_.diffuseColor = tgt::vec3(1,1,1);
    lightSource_.specularColor = tgt::vec3(1,1,1);

    // material parameters
    material_.shininess = 20.0f;

    // initialize sphere geometry
    mesh_.setSphereGeometry(1.0f, tgt::vec3::zero, tgt::vec4::one, 40);
}

DistanceMeasure::~DistanceMeasure() {
}

Processor* DistanceMeasure::create() const {
    return new DistanceMeasure();
}

bool DistanceMeasure::isReady() const {
    if (!isInitialized() || !imgInport_.isReady() || !fhpInport_.isReady() || !outport_.isReady())
        return false;

    return true;
}

tgt::ivec2 DistanceMeasure::clampToViewport(tgt::ivec2 mousePos) {
    tgt::ivec2 result = mousePos;
    tgt::ivec2 size = imgInport_.getSize();
    if (result.x < 0) result.x = 0;
    else if (result.x > size.x-1) result.x = size.x-1;
    if (result.y < 0) result.y = 0;
    else if (result.y > size.y-1) result.y = size.y-1;
    return result;
}

void DistanceMeasure::measure(tgt::MouseEvent* e) {
    const VolumeBase* refVolume = refInport_.getData();
    if(!refVolume) {
        LERROR("No reference volume");
        return;
    }

    if (e->action() & tgt::MouseEvent::PRESSED) {
        if (!mouseDown_) {
            mouseStartPos2D_ = clampToViewport(tgt::ivec2(e->coord().x, e->viewport().y-e->coord().y));
            tgt::vec4 fhp;
            fhpInport_.activateTarget();
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            glReadPixels(mouseStartPos2D_.x, mouseStartPos2D_.y, 1, 1, GL_RGBA, GL_FLOAT, &fhp);
            if (length(fhp) > 0.0f) {
                mouseDown_ = true;
                distance_ = 0.0f;
                mouseStartPos3D_ = refVolume->getTextureToWorldMatrix() * fhp;
                e->accept();
            }
            fhpInport_.deactivateTarget();
        }
    }
    if (e->action() & tgt::MouseEvent::MOTION) {
        if (mouseDown_) {
            // check domain
            tgt::ivec2 oldMouseCurPos2D_ = mouseCurPos2D_;
            mouseCurPos2D_ = clampToViewport(tgt::ivec2(e->coord().x, e->viewport().y-e->coord().y));
            tgt::vec4 fhp;
            fhpInport_.activateTarget();
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            glReadPixels(mouseCurPos2D_.x, mouseCurPos2D_.y, 1, 1, GL_RGBA, GL_FLOAT, &fhp);
            if (length(fhp) > 0.0f) {
                mouseCurPos3D_ = refVolume->getTextureToWorldMatrix() * fhp;
                distance_ = length(mouseStartPos3D_-mouseCurPos3D_);
                invalidate();
            } else {
                mouseCurPos2D_ = oldMouseCurPos2D_;
            }
            e->accept();
            fhpInport_.deactivateTarget();
        }
    }
    if (e->action() & tgt::MouseEvent::RELEASED) {
        if (mouseDown_) {
            mouseDown_ = false;
            mouseStartPos2D_ = tgt::ivec2(0, 0);
            mouseStartPos3D_ = tgt::vec4(0.0f);
            invalidate();
            e->accept();
        }
    }
}

void DistanceMeasure::process() {
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();

    const VolumeBase* refVolume = refInport_.getData();
    if(!refVolume) {
        LERROR("No reference volume");
        return;
    }


    outport_.activateTarget();
    outport_.clearTarget();

    TextureUnit colorUnit, depthUnit;
    imgInport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    imgInport_.setTextureParameters(program_, "textureParameters_");

    renderQuad();

    program_->deactivate();
    TextureUnit::setZeroUnit();

    if (mouseDown_) {
        // render text
        glDisable(GL_DEPTH_TEST);

        // construct the string for display
        std::string label = formatSpatialLength(distance_);

        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.loadIdentity();
        MatStack.translate(-1.0f, -1.0f, 0.0f);
        float scaleFactorX = 2.0f / static_cast<float>(imgInport_.getSize().x);
        float scaleFactorY = 2.0f / static_cast<float>(imgInport_.getSize().y);
        MatStack.scale(scaleFactorX, scaleFactorY, 1);
        tgt::ivec2 screensize = imgInport_.getSize();

        font_.setFontColor(tgt::vec4(0.0f, 0.0f, 0.0f, 1.f));
        font_.render(tgt::vec3(static_cast<float>(mouseCurPos2D_.x+11), static_cast<float>(mouseCurPos2D_.y+11), 0.0f), label, screensize);
        font_.setFontColor(tgt::vec4(0.7f, 0.7f, 0.7f, 1.f));
        font_.render(tgt::vec3(static_cast<float>(mouseCurPos2D_.x+10), static_cast<float>(mouseCurPos2D_.y+10), 0.0f), label, screensize);

        glLineWidth(2.f);
        IMode.begin(tgt::ImmediateMode::LINES);
            IMode.color(0.0f, 0.0f, 0.0f, 1.f);
            IMode.vertex(mouseStartPos2D_.x, mouseStartPos2D_.y);
            IMode.vertex(mouseCurPos2D_.x, mouseCurPos2D_.y);
            IMode.color(1.0f, 1.0f, 1.0f, 1.f);
            IMode.vertex(mouseStartPos2D_.x+1, mouseStartPos2D_.y+1);
            IMode.vertex(mouseCurPos2D_.x+1, mouseCurPos2D_.y+1);
        IMode.end();
        glLineWidth(1.f);

        MatStack.popMatrix();
        glEnable(GL_DEPTH_TEST);

        if(renderSpheres_.get()) {
            IMode.setLightSource(lightSource_);
            IMode.color(0.5f, 0.5f, 0.5f, 0.5f);
            IMode.setMaterial(material_);
            float sphereRadius = tgt::length(refVolume->getCubeSize())*0.005;

            // set modelview and projection matrices
            MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
            MatStack.pushMatrix();
            MatStack.loadMatrix(camera_.get().getProjectionMatrix(outport_.getSize()));

            MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
            MatStack.pushMatrix();
            MatStack.loadIdentity();
            MatStack.loadMatrix(camera_.get().getViewMatrix());

            MatStack.pushMatrix();
            MatStack.translate(mouseStartPos3D_.x, mouseStartPos3D_.y, mouseStartPos3D_.z);
            MatStack.scale(tgt::vec3(sphereRadius));
            mesh_.render();
            MatStack.popMatrix();

            MatStack.pushMatrix();
            MatStack.translate(mouseCurPos3D_.x, mouseCurPos3D_.y, mouseCurPos3D_.z);
            MatStack.scale(tgt::vec3(sphereRadius));
            mesh_.render();
            MatStack.popMatrix();

            MatStack.popMatrix();

            MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
            MatStack.popMatrix();
        }

        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
    }

    outport_.deactivateTarget();
    LGL_ERROR;

}

} // namespace voreen
