/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "slicepointrenderer3d.h"

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "tgt/vector.h"
#include "tgt/immediatemode/immediatemode.h"

namespace {
    tgt::ivec3 INVISIBLE = -tgt::ivec3::one;
}

namespace voreen {

SlicePointRenderer3D::SlicePointRenderer3D()
    : GeometryRendererBase()
    //port
    , inport_(Port::INPORT, "volume", "Volume Input")
    //properties
    , enable_("enable", "Enable", true)
    , pointRadius_("pointRadius","Point Radius",5,2,100)
    , renderLines_("renderLines","Render Lines",true,Processor::INVALID_RESULT,Property::LOD_DEFAULT)
    // points
    , renderPoint0_("renderPoint0", "Render Point1", false)
    , pointColor0_("pointColor0", "Point1 Color", tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f))
    , pointPos0_("pointPos0", "Point1 Position", INVISIBLE, tgt::ivec3(-1), tgt::ivec3(INT_MAX-1),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , renderPoint1_("renderPoint1", "Render Point2", false)
    , pointColor1_("pointColor1", "Point2 Color", tgt::vec4(0.0f, 1.0f, 0.0f, 1.0f))
    , pointPos1_("pointPos1", "Point2 Position", INVISIBLE, tgt::ivec3(-1), tgt::ivec3(INT_MAX-1),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , renderPoint2_("renderPoint2", "Render Point3", false)
    , pointColor2_("pointColor2", "Point3 Color", tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f))
    , pointPos2_("pointPos2", "Point3 Position", INVISIBLE, tgt::ivec3(-1), tgt::ivec3(INT_MAX-1),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , renderPoint3_("renderPoint3", "Render Point4", false)
    , pointColor3_("pointColor3", "Point4 Color", tgt::vec4(1.0f, 1.0f, 0.0f, 1.0f))
    , pointPos3_("pointPos3", "Point4 Position", INVISIBLE, tgt::ivec3(-1), tgt::ivec3(INT_MAX-1),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , shaderProp_("geometry.prg", "Shader", "trianglemesh.frag", "trianglemesh.vert", "trianglemesh.geom",
                  Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , currentGeometry_(0)
{
    //ports
    addPort(inport_);
    //properties
    addProperty(enable_);

    addProperty(pointRadius_);
        ON_PROPERTY_CHANGE(pointRadius_,SlicePointRenderer3D,updateGeometry);
        pointRadius_.setGroupID("general");
    addProperty(renderLines_);
        renderLines_.setGroupID("general");
    setPropertyGroupGuiName("general","General Settings");

    //points
    addProperty(renderPoint0_);
        renderPoint0_.setGroupID("point0");
    addProperty(pointColor0_);
        pointColor0_.setGroupID("point0");
    addProperty(pointPos0_);
        pointPos0_.setGroupID("point0");
    setPropertyGroupGuiName("point0","Point1 Settings");

    addProperty(renderPoint1_);
        renderPoint1_.setGroupID("point1");
    addProperty(pointColor1_);
        pointColor1_.setGroupID("point1");
    addProperty(pointPos1_);
        pointPos1_.setGroupID("point1");
    setPropertyGroupGuiName("point1","Point2 Settings");

    addProperty(renderPoint2_);
        renderPoint2_.setGroupID("point2");
    addProperty(pointColor2_);
        pointColor2_.setGroupID("point2");
    addProperty(pointPos2_);
        pointPos2_.setGroupID("point2");
    setPropertyGroupGuiName("point2","Point3 Settings");

    addProperty(renderPoint3_);
        renderPoint3_.setGroupID("point3");
    addProperty(pointColor3_);
        pointColor3_.setGroupID("point3");
    addProperty(pointPos3_);
        pointPos3_.setGroupID("point3");
    setPropertyGroupGuiName("point3","Point4 Settings");

    //geometry shader
    addProperty(shaderProp_);
}

void SlicePointRenderer3D::initialize() {
    GeometryRendererBase::initialize();
}

void SlicePointRenderer3D::deinitialize() {
    delete currentGeometry_;
    currentGeometry_ = 0;
    GeometryRendererBase::deinitialize();
}

void SlicePointRenderer3D::updateGeometry() {
    delete currentGeometry_;

    currentGeometry_ = new TriangleMeshGeometryUInt16IndexedSimple();
    currentGeometry_->addDiskGeometry(0.f,(float)pointRadius_.get());
}

tgt::Bounds SlicePointRenderer3D::getBoundingBox() const {
    if(inport_.isReady())
        return inport_.getData()->getBoundingBox().getBoundingBox(false);

    return GeometryRendererBase::getBoundingBox();
}

void SlicePointRenderer3D::process() {}

void SlicePointRenderer3D::render() {
    if (!inport_.isReady() || !enable_.get())
        return;

    //init geometry
    if(!currentGeometry_)
        updateGeometry();

    //init shader
    if(!shaderProp_.hasValidShader() || getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        shaderProp_.setHeader(generateHeader() + (currentGeometry_ ? currentGeometry_->getShaderDefines(): ""));
        shaderProp_.rebuild();
    }

    //shader must be valid
    if(!shaderProp_.hasValidShader()) {
        LERROR("Shader for geometry failed to compile");
        return;
    }

    //gl state changes
    glDepthRange(0.f,0.f); //pushes point to front
    glDepthFunc(GL_ALWAYS);

    // check which speheres should be rendered
    bool visiblePoint0, visiblePoint1, visiblePoint2, visiblePoint3;
    tgt::vec3 pointWolrdPos0, pointWolrdPos1, pointWolrdPos2, pointWolrdPos3;
    if((visiblePoint0 = renderPoint0_.get() && pointPos0_.get() != INVISIBLE) && pointColor0_.get().a > 0.f)
        pointWolrdPos0 = pointPos0_.get();
    if((visiblePoint1 = renderPoint1_.get() && pointPos1_.get() != INVISIBLE) && pointColor1_.get().a > 0.f)
        pointWolrdPos1 = pointPos1_.get();
    if((visiblePoint2 = renderPoint2_.get() && pointPos2_.get() != INVISIBLE) && pointColor2_.get().a > 0.f)
        pointWolrdPos2 = pointPos2_.get();
    if((visiblePoint3 = renderPoint3_.get() && pointPos3_.get() != INVISIBLE) && pointColor3_.get().a > 0.f)
        pointWolrdPos3 = pointPos3_.get();


    //render lines
    tgt::vec3 geomLlf(0,0,0);
    tgt::vec3 geomUrb(inport_.getData()->getDimensions());
    tgt::mat4 modelViewMat = camera_.getViewMatrix() * inport_.getData()->getVoxelToWorldMatrix();

    if(renderLines_.get()) {
        MatStack.pushMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadMatrix(modelViewMat);

        IMode.begin(tgt::ImmediateMode::LINES);
        if(visiblePoint0) renderLinesHelper(geomLlf,geomUrb,pointWolrdPos0,pointColor0_.get());
        if(visiblePoint1) renderLinesHelper(geomLlf,geomUrb,pointWolrdPos1,pointColor1_.get());
        if(visiblePoint2) renderLinesHelper(geomLlf,geomUrb,pointWolrdPos2,pointColor2_.get());
        if(visiblePoint3) renderLinesHelper(geomLlf,geomUrb,pointWolrdPos3,pointColor3_.get());
        IMode.end();
        LGL_ERROR;

        MatStack.popMatrix();
    }

    //render points
    tgt::ivec4 viewport;
    glGetIntegerv(GL_VIEWPORT,viewport.elem);
    tgt::mat4 pvMat = MatStack.getProjectionMatrix() * modelViewMat;

    tgt::Shader* prog = shaderProp_.getShader();
    prog->activate();
        prog->setUniform("projectionMatrix_", tgt::mat4::createOrtho(0.f,(float)viewport.z,0.f,(float)viewport.w,-1.f,1.f));
        prog->setUniform("viewMatrix_", tgt::mat4::createIdentity());

        if(visiblePoint0) renderDiskHelper(prog,viewport,pvMat,pointWolrdPos0,pointColor0_.get());
        if(visiblePoint1) renderDiskHelper(prog,viewport,pvMat,pointWolrdPos1,pointColor1_.get());
        if(visiblePoint2) renderDiskHelper(prog,viewport,pvMat,pointWolrdPos2,pointColor2_.get());
        if(visiblePoint3) renderDiskHelper(prog,viewport,pvMat,pointWolrdPos3,pointColor3_.get());

    prog->deactivate();
    LGL_ERROR;

    //set gl state back to default
    glDepthRange(0.f,1.f);
    glDepthFunc(GL_LESS);
    IMode.color(tgt::vec4::one);
    LGL_ERROR;
}

void SlicePointRenderer3D::renderDiskHelper(tgt::Shader* prog, const tgt::ivec4 viewport, const tgt::mat4& pvMat, const tgt::vec3& pointPos, const tgt::vec4 col) {
    tgt::vec4 clipSpace = pvMat * tgt::vec4(pointPos,1.f);
    if (clipSpace.w == 0.f)
        return;
    tgt::vec3 ndc = tgt::vec3(clipSpace.xyz() / clipSpace.w);
    tgt::vec3 pos = tgt::vec3(((ndc.xy() + tgt::vec2::one)/2.f) * tgt::vec2(viewport.zw()) + tgt::vec2(viewport.xy()),0.f);
    prog->setUniform("modelMatrix_", tgt::mat4::createTranslation(pos));
    prog->setUniform("solidColor_", tgt::vec4(col));
    currentGeometry_->render();
}


void SlicePointRenderer3D::renderLinesHelper(const tgt::vec3& geomLlf, const tgt::vec3& geomUrb, const tgt::vec3& pointPos, const tgt::vec4& col) {
    IMode.color(col);
    IMode.vertex(tgt::vec3(geomLlf.x,pointPos.y,pointPos.z));
    IMode.vertex(tgt::vec3(geomUrb.x,pointPos.y,pointPos.z));
    IMode.vertex(tgt::vec3(pointPos.x,geomLlf.y,pointPos.z));
    IMode.vertex(tgt::vec3(pointPos.x,geomUrb.y,pointPos.z));
    IMode.vertex(tgt::vec3(pointPos.x,pointPos.y,geomLlf.z));
    IMode.vertex(tgt::vec3(pointPos.x,pointPos.y,geomUrb.z));
}

void SlicePointRenderer3D::invalidate(int inv) {
    GeometryRendererBase::invalidate(inv);

    if (inport_.hasChanged() && inport_.hasData()) {
        tgt::ivec3 dim = tgt::ivec3(inport_.getData()->getDimensions());
        pointPos0_.setMaxValue(dim-tgt::ivec3::one);
        pointPos1_.setMaxValue(dim-tgt::ivec3::one);
        pointPos2_.setMaxValue(dim-tgt::ivec3::one);
        pointPos3_.setMaxValue(dim-tgt::ivec3::one);
    }
}

}

