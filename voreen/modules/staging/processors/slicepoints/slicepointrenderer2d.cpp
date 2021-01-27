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

#include "slicepointrenderer2d.h"

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/textureunit.h"

namespace voreen {

const std::string SlicePointRenderer2D::loggerCat_("voreen.SlicePointRenderer2D");

SlicePointRenderer2D::SlicePointRenderer2D()
    : ImageProcessor("image/background")
    , inport_(Port::INPORT,"inport","Slice Inport")
    , outport_(Port::OUTPORT,"outport","Slice Outport")
    , privatePort_(Port::OUTPORT,"private","Point Port")
    //properties
    , enable_("enable", "Enable", true)
    , pointRadius_("pointRadius","Point Radius",5,2,100)
    , renderPointInfo_("renderPointInfo","Point Info ?",true)
    // points
    , renderPoint0_("renderPoint0", "Render Point1", false)
    , pointColor0_("pointColor0", "Point1 Color", tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f))
    , pointPos0_("pointPos0", "Point1 Position", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(INT_MAX),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , renderPoint1_("renderPoint1", "Render Point2", false)
    , pointColor1_("pointColor1", "Point2 Color", tgt::vec4(0.0f, 1.0f, 0.0f, 1.0f))
    , pointPos1_("pointPos1", "Point2 Position", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(INT_MAX),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , renderPoint2_("renderPoint2", "Render Point3", false)
    , pointColor2_("pointColor2", "Point3 Color", tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f))
    , pointPos2_("pointPos2", "Point3 Position", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(INT_MAX),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , renderPoint3_("renderPoint3", "Render Point4", false)
    , pointColor3_("pointColor3", "Point4 Color", tgt::vec4(1.0f, 1.0f, 0.0f, 1.0f))
    , pointPos3_("pointPos3", "Point4 Position", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(INT_MAX),
                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    //to be linked
    , linkedMousePositionInSlice_("linkedMousePositionInSlice","\"Mouse Position\" from SliceViewer",tgt::ivec3(-1),tgt::ivec3(-1),tgt::ivec3(INT_MAX),
                                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , linkedSliceAlignment_("linkedSliceAlignment", "\"Slice Alignment\" from SliceViewer",Processor::INVALID_RESULT,false,Property::LOD_DEBUG)
    , linkedSliceIndex_("linkedSliceIndex", "\"Slice Number\" from SliceViewer", 0, 0, 10000,Processor::INVALID_RESULT,NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , linkedPickingMatrix_("linkedPickingMatrix", "\"Picking Matrix\" from SliceViewer", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f),
                            Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC,Property::LOD_DEBUG)
    // shader
    , shaderProp_("geometry.prg", "Shader", "trianglemesh.frag", "trianglemesh.vert", "trianglemesh.geom",
                  Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , point0EnableKeyEvent_(0), point1EnableKeyEvent_(0), point2EnableKeyEvent_(0), point3EnableKeyEvent_(0)
    , point0ActivateKeyEvent_(0), point1ActivateKeyEvent_(0), point2ActivateKeyEvent_(0), point3ActivateKeyEvent_(0)
    , activePoint_(-1), equalGeometry_(0), moreGeometry_(0), lessGeometry_(0), pointInfoFont_(0), fontLineWidth_(200.f)
{
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&privatePort_);

    //general
    addProperty(enable_);
    addProperty(pointRadius_);
        ON_PROPERTY_CHANGE(pointRadius_,SlicePointRenderer2D,updateGeometry);
        pointRadius_.setGroupID("general");
    addProperty(renderPointInfo_);
        renderPointInfo_.setGroupID("general");
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

    //to be linked
    addProperty(linkedMousePositionInSlice_);
        ON_PROPERTY_CHANGE(linkedMousePositionInSlice_,SlicePointRenderer2D,linkedMousePositionInSliceOnChange);
        linkedMousePositionInSlice_.setGroupID("linked");
    addProperty(linkedSliceAlignment_);
        linkedSliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
        linkedSliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
        linkedSliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
        linkedSliceAlignment_.setGroupID("linked");
    addProperty(linkedSliceIndex_);
        linkedSliceIndex_.setGroupID("linked");
    addProperty(linkedPickingMatrix_);
        linkedPickingMatrix_.setReadOnlyFlag(true);
        linkedPickingMatrix_.setGroupID("linked");
    setPropertyGroupGuiName("linked","To be linked");

    // shader
    addProperty(shaderProp_);
        shaderProp_.setVisibleFlag(false);
    // create font
    pointInfoFont_ = new tgt::Font(VoreenApplication::app()->getFontPath("VeraMono.ttf"),10,
                                   fontLineWidth_,tgt::Font::Centered);

    //events
    point0EnableKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point0EnableKey", "Point 0 Enable",
        this, &SlicePointRenderer2D::enablePoint0KeyEvent, tgt::KeyEvent::K_1, tgt::Event::CTRL, true);
    point1EnableKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point1EnableKey", "Point 1 Enable",
        this, &SlicePointRenderer2D::enablePoint1KeyEvent, tgt::KeyEvent::K_2, tgt::Event::CTRL, true);
    point2EnableKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point2EnableKey", "Point 2 Enable",
        this, &SlicePointRenderer2D::enablePoint2KeyEvent, tgt::KeyEvent::K_3, tgt::Event::CTRL, true);
    point3EnableKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point3EnableKey", "Point 3 Enable",
        this, &SlicePointRenderer2D::enablePoint3KeyEvent, tgt::KeyEvent::K_4, tgt::Event::CTRL, true);
    point0ActivateKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point0ActivateKey", "Point 0 Activate",
        this, &SlicePointRenderer2D::activatePoint0KeyEvent, tgt::KeyEvent::K_1, tgt::Event::MODIFIER_NONE, true);
    point1ActivateKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point1ActivateKey", "Point 1 Activate",
        this, &SlicePointRenderer2D::activatePoint1KeyEvent, tgt::KeyEvent::K_2, tgt::Event::MODIFIER_NONE, true);
    point2ActivateKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point2ActivateKey", "Point 2 Activate",
        this, &SlicePointRenderer2D::activatePoint2KeyEvent, tgt::KeyEvent::K_3, tgt::Event::MODIFIER_NONE, true);
    point3ActivateKeyEvent_ = new EventProperty<SlicePointRenderer2D>("keyEvent.point3ActivateKey", "Point 3 Activate",
        this, &SlicePointRenderer2D::activatePoint3KeyEvent, tgt::KeyEvent::K_4, tgt::Event::MODIFIER_NONE, true);

    addEventProperty(point0ActivateKeyEvent_);
    addEventProperty(point1ActivateKeyEvent_);
    addEventProperty(point2ActivateKeyEvent_);
    addEventProperty(point3ActivateKeyEvent_);
    addEventProperty(point0EnableKeyEvent_);
    addEventProperty(point1EnableKeyEvent_);
    addEventProperty(point2EnableKeyEvent_);
    addEventProperty(point3EnableKeyEvent_);
}

SlicePointRenderer2D::~SlicePointRenderer2D() {
    delete pointInfoFont_;
    delete point0EnableKeyEvent_;
    delete point1EnableKeyEvent_;
    delete point2EnableKeyEvent_;
    delete point3EnableKeyEvent_;
    delete point0ActivateKeyEvent_;
    delete point1ActivateKeyEvent_;
    delete point2ActivateKeyEvent_;
    delete point3ActivateKeyEvent_;
}

void SlicePointRenderer2D::deinitialize() {
    delete equalGeometry_;
    delete moreGeometry_;
    delete lessGeometry_;
    equalGeometry_ = 0;
    moreGeometry_ = 0;
    lessGeometry_ = 0;
    ImageProcessor::deinitialize();
}

std::string SlicePointRenderer2D::generateHeader(const tgt::GpuCapabilities::GlVersion* version /*= 0*/) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define BLEND_MODE_ADD\n";
    return header;
}

void SlicePointRenderer2D::process() {
    //init geometry
    if(!equalGeometry_ || !moreGeometry_ || !lessGeometry_)
        updateGeometry();

    //init shader
    if(!shaderProp_.hasValidShader() || getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        shaderProp_.setHeader(generateHeader() + (equalGeometry_ ? equalGeometry_->getShaderDefines(): "")); //all three geometries should be the same
        shaderProp_.rebuild();
    }

    //shader must be valid
    if(!shaderProp_.hasValidShader()) {
        LERROR("Shader for geometry failed to compile");
        return;
    }

    //render the points
    renderPoints();

    //combine both images
    blendPointsOverInport();
    LGL_ERROR;
}

void SlicePointRenderer2D::blendPointsOverInport() {
    //activate outport
    outport_.activateTarget();
        //clear old image
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // bind textures
        tgt::TextureUnit colorUnit0, depthUnit0, colorUnit1, depthUnit1;
        inport_.bindTextures(colorUnit0.getEnum(), depthUnit0.getEnum());
        privatePort_.bindTextures(colorUnit1.getEnum(), depthUnit1.getEnum());

        // activate shader and set parameteres
        program_->activate();
            setGlobalShaderParameters(program_);
            program_->setUniform("colorTex0_", colorUnit0.getUnitNumber());
            program_->setUniform("depthTex0_", depthUnit0.getUnitNumber());
            program_->setUniform("colorTex1_", colorUnit1.getUnitNumber());
            program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
            inport_.setTextureParameters(program_, "textureParameters0_");
            privatePort_.setTextureParameters(program_, "textureParameters1_");

        // blend images
        renderQuad();

        //deactivate shader
        program_->deactivate();
    //deactivate outport
    outport_.deactivateTarget();

    //clean up
    tgt::TextureUnit::setZeroUnit();
}

void SlicePointRenderer2D::renderPoints() {
    privatePort_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //return if disabled
    if(!enable_.get()) {
        privatePort_.deactivateTarget();
        return;
    }

    //render points
    tgt::mat4 invPicking; //from volume to viewport
    if(!linkedPickingMatrix_.get().invert(invPicking)) {
        tgtAssert(false,"Could not invert picking matrix.");
        return;
    }

    //determine, if point should be rendered
    bool visiblePoint0 = renderPoint0_.get() && pointColor0_.get().a > 0.f;
    bool visiblePoint1 = renderPoint1_.get() && pointColor1_.get().a > 0.f;
    bool visiblePoint2 = renderPoint2_.get() && pointColor2_.get().a > 0.f;
    bool visiblePoint3 = renderPoint3_.get() && pointColor3_.get().a > 0.f;

    //calc screenPos
    tgt::vec4 screenPos0, screenPos1, screenPos2, screenPos3;
    if(visiblePoint0) screenPos0 = invPicking * tgt::vec4(pointPos0_.get(),1.f);
    if(visiblePoint1) screenPos1 = invPicking * tgt::vec4(pointPos1_.get(),1.f);
    if(visiblePoint2) screenPos2 = invPicking * tgt::vec4(pointPos2_.get(),1.f);
    if(visiblePoint3) screenPos3 = invPicking * tgt::vec4(pointPos3_.get(),1.f);

    //render geometry
    tgt::Shader* prog = shaderProp_.getShader();
    prog->activate();
        prog->setUniform("projectionMatrix_", tgt::mat4::createOrtho(0.f, (float)privatePort_.getSize().x, 0.f, (float)privatePort_.getSize().y, -1.f, 1.f));
        prog->setUniform("viewMatrix_", tgt::mat4::createIdentity());

        if(visiblePoint0) renderPointGeometryHelper(prog, pointPos0_.get(), screenPos0, pointColor0_.get());
        if(visiblePoint1) renderPointGeometryHelper(prog, pointPos1_.get(), screenPos1, pointColor1_.get());
        if(visiblePoint2) renderPointGeometryHelper(prog, pointPos2_.get(), screenPos2, pointColor2_.get());
        if(visiblePoint3) renderPointGeometryHelper(prog, pointPos3_.get(), screenPos3, pointColor3_.get());
    prog->deactivate();


    //render info
    if(renderPointInfo_.get()) {
            if(visiblePoint0) renderPointInfoHelper(pointPos0_.get(), screenPos0, pointColor0_.get());
            if(visiblePoint1) renderPointInfoHelper(pointPos1_.get(), screenPos1, pointColor1_.get());
            if(visiblePoint2) renderPointInfoHelper(pointPos2_.get(), screenPos2, pointColor2_.get());
            if(visiblePoint3) renderPointInfoHelper(pointPos3_.get(), screenPos3, pointColor3_.get());
    }

    privatePort_.deactivateTarget();
    LGL_ERROR;
}

void SlicePointRenderer2D::renderPointGeometryHelper(tgt::Shader* prog, const tgt::ivec3& pointPos, const tgt::vec4& screenPos, const tgt::vec4& col) {
    prog->setUniform("modelMatrix_",tgt::mat4::createTranslation(tgt::vec3(screenPos.xy(),0.f)));
    prog->setUniform("solidColor_", col);
    // check if we are in the correct slice
    int dif = pointPos[static_cast<int>(linkedSliceAlignment_.getValue())] - linkedSliceIndex_.get();
    if(dif < 0) {
        lessGeometry_->render();
    } else {
        if(dif == 0)
            equalGeometry_->render();
        else
            moreGeometry_->render();
    }
}

void SlicePointRenderer2D::renderPointInfoHelper(const tgt::ivec3& pointPos, const tgt::vec4& screenPos, const tgt::vec4& col) {
    // check if we are in the correct slice
    int dif = pointPos[static_cast<int>(linkedSliceAlignment_.getValue())] - linkedSliceIndex_.get();

    //no font if we are in the current slice
    if(dif == 0) return;

    //adapt font
    int xpos = 0;
    if(screenPos.x < privatePort_.getSize().x/2) {
        xpos = screenPos.x-pointRadius_.get();
        pointInfoFont_->setTextAlignment(tgt::Font::MiddleLeft);
    } else {
        xpos = screenPos.x+pointRadius_.get()-fontLineWidth_;
        pointInfoFont_->setTextAlignment(tgt::Font::MiddleRight);
    }

    pointInfoFont_->setFontColor(col);
    pointInfoFont_->render(tgt::vec3(xpos,screenPos.y + ((dif < 0) ? 2.5 : -2.5 ) *pointRadius_.get(),0), itos(std::abs(dif)), privatePort_.getSize());
}

//-----------------------------------------------------------------------------------------
//                          event callback functions
//-----------------------------------------------------------------------------------------
void SlicePointRenderer2D::activatePoint0KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(0,&renderPoint0_);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}
void SlicePointRenderer2D::activatePoint1KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(1,&renderPoint1_);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}
void SlicePointRenderer2D::activatePoint2KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(2,&renderPoint2_);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}
void SlicePointRenderer2D::activatePoint3KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(3,&renderPoint3_);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}
void SlicePointRenderer2D::enablePoint0KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(0,&renderPoint0_,true);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}
void SlicePointRenderer2D::enablePoint1KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(1,&renderPoint1_,true);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}
void SlicePointRenderer2D::enablePoint2KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(2,&renderPoint2_,true);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}
void SlicePointRenderer2D::enablePoint3KeyEvent(tgt::KeyEvent* e) {
    //do not handle auto repeated keys
    if(e->autoRepeat() || !e->pressed()) {e->ignore(); return;}
    mainEventFunc(3,&renderPoint3_,true);
    //accept event
    e->accept();
    //invalidate rendering
    invalidate();
}

void SlicePointRenderer2D::mainEventFunc(int point, BoolProperty* enableProp, bool enable) {
    //switch enable if property has been passed
    if(enable) {
        enableProp->set(!enableProp->get());
        if(!enableProp->get())
            activePoint_ = -1;
    } else { //active points
        if(!enableProp->get()) return;
        if(activePoint_ == point)
            activePoint_ = -1;
        else {
            activePoint_ = point;
            //force point update
            linkedMousePositionInSliceOnChange();
        }
    }
}

//-----------------------------------------------------------------------------------------
//                          on chance callback functions
//-----------------------------------------------------------------------------------------
void SlicePointRenderer2D::linkedMousePositionInSliceOnChange() {
    //if mouse pos out of range do nothing
    if(tgt::hor(tgt::lessThan(linkedMousePositionInSlice_.get(),tgt::ivec3::zero)))
        return;
    //set mouse pos
    switch(activePoint_) {
    case 0:
        pointPos0_.set(linkedMousePositionInSlice_.get());
        break;
    case 1:
        pointPos1_.set(linkedMousePositionInSlice_.get());
        break;
    case 2:
        pointPos2_.set(linkedMousePositionInSlice_.get());
        break;
    case 3:
        pointPos3_.set(linkedMousePositionInSlice_.get());
        break;
    default:
        //nothing to do here
        break;
    }
}

void SlicePointRenderer2D::updateGeometry() {
    delete equalGeometry_;
    delete moreGeometry_;
    delete lessGeometry_;

    equalGeometry_ = new TriangleMeshGeometryUInt16IndexedSimple();
    equalGeometry_->addDiskGeometry(pointRadius_.get(),pointRadius_.get()+(static_cast<float>(pointRadius_.get())/2.f));

    moreGeometry_ = new TriangleMeshGeometryUInt16IndexedSimple();
    std::vector<VertexBase> vec(3);
    vec[0] = VertexBase(tgt::vec3(0.f,pointRadius_.get(),0.f));
    vec[1] = VertexBase(tgt::vec3((float)pointRadius_.get(), (float)-pointRadius_.get(), 0.f));
    vec[2] = VertexBase(tgt::vec3((float)-pointRadius_.get(), (float)-pointRadius_.get(), 0.f));
    moreGeometry_->setVertices(vec);
    moreGeometry_->addTriangle(tgt::ivec3(0,1,2));

    lessGeometry_ = new TriangleMeshGeometryUInt16IndexedSimple();
    vec[0] = VertexBase(tgt::vec3(0.f,-pointRadius_.get(),0.f));
    vec[1] = VertexBase(tgt::vec3((float)-pointRadius_.get(), (float)pointRadius_.get(), 0.f));
    vec[2] = VertexBase(tgt::vec3((float)pointRadius_.get(), (float)pointRadius_.get(), 0.f));
    lessGeometry_->setVertices(vec);
    lessGeometry_->addTriangle(tgt::ivec3(0,1,2));

    pointInfoFont_->setFontSize(2*(pointRadius_.get()+1)); //+1 for optical reasons
    fontLineWidth_ = static_cast<float>(std::min(pointInfoFont_->getFontSize()*2,300));
    pointInfoFont_->setLineWidth(fontLineWidth_);
}


} // namespace voreen
