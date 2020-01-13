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

#include "interactiveregistrationwidget.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "tgt/textureunit.h"
#include "tgt/glmath.h"
#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

using tgt::ivec2;
using tgt::vec3;
using tgt::vec4;
using tgt::mat4;

const std::string InteractiveRegistrationWidget::loggerCat_("voreen.InteractiveRegistrationWidget");

InteractiveRegistrationWidget::InteractiveRegistrationWidget()
    : RenderProcessor()
    , inport_(Port::INPORT, "input", "Image Input")
    , outport_(Port::OUTPORT, "output", "Image Output")
    , pickingPort_(Port::OUTPORT, "picking")
    , textPort_(Port::OUTPORT, "text", "Text Output")
    , transformMatrix_("transformMatrix", "Transformation Matrix", tgt::mat4::identity, tgt::mat4(-2000.0), tgt::mat4(2000.0))
    , point_("point", "Point", vec3(0.0f), vec3(-999999.9f), vec3(999999.9f))
    , plane_("plane", "Plane", vec3(1.0f, 0.0f, 0.0f), vec3(-5.0f), vec3(5.0f), Processor::VALID)
    , planeDist_("planeDist", "Plane Distance", 0.0f, -1000.0f, 1000.0f, Processor::VALID)
    , render_("render", "Render", true)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
    , sphereRadius_("sphereRadius", "Sphere Radius", 5.0f, 1.0f, 100.0f)
    , ringRadius_("ringRadius", "Ring Radius", 50.0f, 1.0f, 300.0f)
    , ringColor_("ringColor", "Ring Color", vec4(0.0f, 1.0f, 0.0f, 0.8f))
    , sphereColor_("sphereColor", "Sphere Color", vec4(0.0f, 0.0f, 1.0f, 0.8f))
    , centerPoint_("centerWidget", "Center Widget", Processor::VALID)
    , copyShader_(0)
    , lastCoord_(0)
    , startDragCoord_(0.0f)
    , curDragCoord_(0.0f)
    , rotAngle_(0.0f)
    , mouseDown_(-1)
    , sphereShader_(0)
{
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(pickingPort_);
    addPort(textPort_);

    addProperty(render_);
    addProperty(transformMatrix_);
    addProperty(point_);
    addProperty(plane_);
    addProperty(planeDist_);
    addProperty(camera_);
    addProperty(sphereRadius_);
    addProperty(ringRadius_);
    addProperty(ringColor_);
    addProperty(sphereColor_);
    addProperty(centerPoint_);

    plane_.onChange(MemberFunctionCallback<InteractiveRegistrationWidget>(this, &InteractiveRegistrationWidget::planeChanged));
    planeDist_.onChange(MemberFunctionCallback<InteractiveRegistrationWidget>(this, &InteractiveRegistrationWidget::planeChanged));
    centerPoint_.onChange(MemberFunctionCallback<InteractiveRegistrationWidget>(this, &InteractiveRegistrationWidget::centerPoint));
}

InteractiveRegistrationWidget::~InteractiveRegistrationWidget() {}

void InteractiveRegistrationWidget::initialize() {
    RenderProcessor::initialize();

    try {
        copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
            RenderProcessor::generateHeader(), false);
        copyShader_->deactivate();
        sphereShader_ = ShdrMgr.load("pointlistrenderer_spheres", "", false);
    } catch(tgt::Exception&) {
        processorState_ = PROCESSOR_STATE_NOT_INITIALIZED;
        throw VoreenException("Failed to load shader: passthrough.vert/copyimage.frag");
    }
}

void InteractiveRegistrationWidget::deinitialize() {
    ShdrMgr.dispose(copyShader_);
    copyShader_ = 0;
    ShdrMgr.dispose(sphereShader_);
    sphereShader_ = 0;
    RenderProcessor::deinitialize();
}

Processor* InteractiveRegistrationWidget::create() const {
    return new InteractiveRegistrationWidget();
}

bool InteractiveRegistrationWidget::isReady() const {
    return (inport_.isReady() && outport_.isReady());
}

void InteractiveRegistrationWidget::process() {
    tgtAssert(outport_.isReady(), "Outport not ready");
    tgtAssert(inport_.isReady(), "Inport not ready");

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // bind input image to tex unit
    tgt::TextureUnit imageUnit, imageUnitDepth;
    inport_.bindTextures(imageUnit.getEnum(), imageUnitDepth.getEnum());

    // 1. copy input image to outport
    copyShader_->activate();
    setGlobalShaderParameters(copyShader_);
    inport_.setTextureParameters(copyShader_, "texParams_");
    copyShader_->setUniform("colorTex_", imageUnit.getUnitNumber());
    copyShader_->setUniform("depthTex_", imageUnitDepth.getUnitNumber());
    renderQuad();
    copyShader_->deactivate();
    LGL_ERROR;

    tgt::vec2 c = camera_.get().project(outport_.getSize(), point_.get()).xy();
    int steps = 30;
    float lw = 5.0f;
    float r = sphereRadius_.get();
    float as = 0.8f; //arrow size
    float cs = 0.4f; //cross size

    GlMeshGeometryUInt32Simple* sphere = 0;
    if (render_.get()) {
        //make new sphere
        sphere = new GlMeshGeometryUInt32Simple();
        sphere->setSphereGeometry(sphereRadius_.get(),tgt::vec3::zero,tgt::vec4::zero,32);

        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadMatrix(tgt::mat4::createOrtho(0, outport_.getSize().x, 0, outport_.getSize().y, -100.0f, 100.0f));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        LGL_ERROR;

        glDepthFunc(GL_ALWAYS);
        glEnable(GL_BLEND);

        // Center of rotation sphere:
        MatStack.translate(c.x, c.y, 0.0f);
        // activate shader
        sphereShader_->activate();

        // set matrix stack uniforms
        sphereShader_->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        sphereShader_->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
        tgt::mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        sphereShader_->setUniform("modelViewMatrixInverseStack_", viewInverse);

        // set lighting parameters
        sphereShader_->setUniform("lightingEnabled_",false);

        sphereShader_->setUniform("color_",sphereColor_.get());

        sphere->render();

        sphereShader_->deactivate();

        // rotation circle:
        glLineWidth(lw);
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.color(ringColor_.get());
        for(int i=0; i<steps; i++) {
            float angle = tgt::PIf * 2.0f * (float)i/(float)steps;
            IMode.vertex(cos(angle)*ringRadius_.get(), sin(angle)*ringRadius_.get());
        }
        IMode.color(tgt::vec4::one);
        IMode.end();

        // rotation angle indicator:
        if(mouseDown_ == 2) {
            tgt::vec2 a = camera_.get().project(outport_.getSize(), startDragCoord_).xy();
            tgt::vec2 b = camera_.get().project(outport_.getSize(), curDragCoord_).xy();
            a -= c;
            b -= c;

            glLineWidth(lw*0.5f);
            IMode.begin(tgt::ImmediateMode::LINES);
                IMode.color(ringColor_.get());
                IMode.vertex(0.0f, 0.0f);
                IMode.vertex(a);
                IMode.vertex(0.0f, 0.0f);
                IMode.vertex(b);
                IMode.color(tgt::vec4::one);
            IMode.end();
        }
        glLineWidth(1.0f);

        // translation cross-arrow-thing:
        MatStack.translate(sphereRadius_.get()*3.0f, sphereRadius_.get()*3.0f, 0.0f);

        IMode.begin(tgt::ImmediateMode::TRIANGLES);
            IMode.color(sphereColor_.get());
            // Arrows:
            IMode.vertex(r*2.0f, 0.0f);
            IMode.vertex(r, -r*as);
            IMode.vertex(r, r*as);

            IMode.vertex(r*-2.0f, 0.0f);
            IMode.vertex(-r, -r*as);
            IMode.vertex(-r, r*as);

            IMode.vertex(0.0f, r*-2.0f);
            IMode.vertex(-r*as, -r);
            IMode.vertex(r*as, -r);

            IMode.vertex(0.0f, r*2.0f);
            IMode.vertex(-r*as, r);
            IMode.vertex(r*as, r);

            // Cross:
            IMode.vertex(r, cs*r);
            IMode.vertex(r, -cs*r);
            IMode.vertex(-r, -cs*r);

            IMode.vertex(-r, cs*r);
            IMode.vertex(-r, -cs*r);
            IMode.vertex(r, cs*r);

            IMode.vertex(cs*r, r);
            IMode.vertex(-cs*r, r);
            IMode.vertex(-cs*r, -r);

            IMode.vertex(cs*r, -r);
            IMode.vertex(-cs*r, -r);
            IMode.vertex(cs*r, r);
            IMode.color(tgt::vec4::one);
        IMode.end();

        glDisable(GL_BLEND);
        glDepthFunc(GL_LESS);

        // delete is called in the picking mode
        // delete sphere;

        LGL_ERROR;

        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
    }

    outport_.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();
    LGL_ERROR;

    if (render_.get()) {
        // sphere is still present
        //GlMeshGeometryUInt32Color* sphere = new GlMeshGeometryUInt32Color();
        //sphere->setSphereGeometry(sphereRadius_.get(),tgt::vec3::zero,tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f),32);


        pickingPort_.activateTarget("picking");
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDepthFunc(GL_ALWAYS);
        glEnable(GL_BLEND);

        // Center of rotation sphere:
        MatStack.translate(c.x, c.y, 0.0f);
        // activate shader
        sphereShader_->activate();

        // set matrix stack uniforms
        sphereShader_->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        sphereShader_->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
        tgt::mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        sphereShader_->setUniform("modelViewMatrixInverseStack_", viewInverse);

        // set lighting parameters
        sphereShader_->setUniform("lightingEnabled_",false);

        sphereShader_->setUniform("color_",tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f));
        sphere->render();
        sphereShader_->deactivate();

        // rotation circle:
        glLineWidth(lw);
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.color(0.0f, 1.0f, 0.0f, 1.0f);
        for(int i=0; i<steps; i++) {
            float angle = tgt::PIf * 2.0f * (float)i/(float)steps;
            IMode.vertex(cos(angle)*ringRadius_.get(), sin(angle)*ringRadius_.get());
        }
        IMode.color(tgt::vec4::one);
        IMode.end();
        glLineWidth(1.0f);

        // translation cross-arrow-thing:
        MatStack.translate(sphereRadius_.get()*3.0f, sphereRadius_.get()*3.0f, 0.0f);

        IMode.begin(tgt::ImmediateMode::TRIANGLES);
            IMode.color(0.0f, 0.0f, 1.0f, 1.0f);
            // Arrows:
            IMode.vertex(r*2.0f, 0.0f);
            IMode.vertex(r, -r*as);
            IMode.vertex(r, r*as);

            IMode.vertex(r*-2.0f, 0.0f);
            IMode.vertex(-r, -r*as);
            IMode.vertex(-r, r*as);

            IMode.vertex(0.0f, r*-2.0f);
            IMode.vertex(-r*as, -r);
            IMode.vertex(r*as, -r);

            IMode.vertex(0.0f, r*2.0f);
            IMode.vertex(-r*as, r);
            IMode.vertex(r*as, r);

            // Cross:
            IMode.vertex(r, cs*r);
            IMode.vertex(r, -cs*r);
            IMode.vertex(-r, -cs*r);

            IMode.vertex(-r, cs*r);
            IMode.vertex(-r, -cs*r);
            IMode.vertex(r, cs*r);

            IMode.vertex(cs*r, r);
            IMode.vertex(-cs*r, r);
            IMode.vertex(-cs*r, -r);

            IMode.vertex(cs*r, -r);
            IMode.vertex(-cs*r, -r);
            IMode.vertex(cs*r, r);
            IMode.color(tgt::vec4::one);
        IMode.end();

        glDisable(GL_BLEND);
        glDepthFunc(GL_LESS);

        //now the sphere is deleted
        delete sphere;

        pickingPort_.deactivateTarget();
        LGL_ERROR;
    }

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;
}

void InteractiveRegistrationWidget::translate(tgt::vec3 v) {
    mat4 m = transformMatrix_.get();
    mat4 m2 = mat4::createTranslation(v);
    transformMatrix_.set(m2 * m);
}

void InteractiveRegistrationWidget::rotate(tgt::vec3 center, tgt::vec3 v, float angle) {
    mat4 m = transformMatrix_.get();

    mat4 t1 = mat4::createTranslation(-center);
    mat4 t2 = mat4::createTranslation(center);
    mat4 rot = mat4::createRotation(angle, v);

    mat4 m2 = t2 * rot * t1;
    transformMatrix_.set(m2 * m);
}

void InteractiveRegistrationWidget::planeChanged() {
    textPort_.setData("");
    vec3 p = point_.get();

    vec3 n = normalize(plane_.get());
    tgt::plane pl(-n, planeDist_.get());
    float dist = pl.distance(p);
    p -= dist * n;

    if(p != point_.get())
        point_.set(p);
}

void InteractiveRegistrationWidget::centerPoint() {
    tgt::line3 l = camera_.get().getViewRay(outport_.getSize(), outport_.getSize() / 2);
    tgt::plane p(normalize(plane_.get()), planeDist_.get());

    float t;
    if(p.intersect(l, t))
        point_.set(l.getFromParam(t));
}

void InteractiveRegistrationWidget::onEvent(tgt::Event* e) {
    e->ignore();

    if(render_.get()) {
        tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);

        if(me) {
            if (me->action() & tgt::MouseEvent::PRESSED) {
                vec4 c = pickingPort_.getRenderTarget()->getColorAtPos(ivec2(me->x(), me->viewport().y - me->y()));

                if(c.a > 0.5f) {
                    if(c.r > 0.5f) {
                        mouseDown_ = 0;
                        e->accept();
                    }
                    else if(c.g > 0.5f) {
                        mouseDown_ = 2;
                        e->accept();
                        rotAngle_ = 0.0f;
                    }
                    else if(c.b > 0.5f) {
                        mouseDown_ = 1;
                        e->accept();
                    }
                    textPort_.setData("");

                    tgt::line3 l = camera_.get().getViewRay(me->viewport(), me->coord());
                    tgt::plane p(normalize(plane_.get()), planeDist_.get());

                    float t;
                    if(p.intersect(l, t)) {
                        // Calculate intersection points:
                        startDragCoord_ = l.getFromParam(t);
                        curDragCoord_ = startDragCoord_;
                    }
                }
            }
            else if (me->action() & tgt::MouseEvent::MOTION) {
                if(mouseDown_ >= 0) {
                    tgt::line3 l = camera_.get().getViewRay(me->viewport(), me->coord());
                    tgt::line3 prevL = camera_.get().getViewRay(me->viewport(), lastCoord_);
                    tgt::plane p(normalize(-plane_.get()), planeDist_.get());

                    float t, prevT;
                    if(p.intersect(l, t) && p.intersect(prevL, prevT)) {
                        // Calculate intersection points:
                        tgt::vec3 i = l.getFromParam(t);
                        tgt::vec3 prevI = prevL.getFromParam(prevT);
                        curDragCoord_ = i;

                        if(mouseDown_ == 0) {
                            point_.set(point_.get() + i - prevI);
                        }
                        else if(mouseDown_ == 1) {
                            translate(i - prevI);
                            point_.set(point_.get() + i - prevI);
                        }
                        else if(mouseDown_ == 2) {
                            vec3 a = normalize(i - point_.get());
                            vec3 b = normalize(prevI - point_.get());
                            float angle = acos(dot(a, b));

                            // Check in which direction we need to rotate:
                            vec3 test = cross(a, b);
                            vec3 toCam = camera_.get().getPosition() - point_.get();
                            if(dot(test, toCam) < 0.0f)
                                angle *= -1.0f;

                            rotAngle_ += angle;
                            rotate(point_.get(), plane_.get(), angle);

                            std::stringstream ss;
                            ss << (rotAngle_ / (2.0f * tgt::PIf)) * 360.0f << "°";

                            textPort_.setData(ss.str());
                        }
                    }
                    e->accept();
                }
            }
            else if (me->action() & tgt::MouseEvent::RELEASED) {
                if(mouseDown_ >= 0) {
                    mouseDown_ = -1;
                    e->accept();
                }
            }

            lastCoord_ = me->coord();
        }
    }

    if(!e->accepted_)
        RenderProcessor::onEvent(e);
}
}   // namespace
