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

#include "planewidgetprocessor.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include "tgt/glmath.h"

using tgt::vec3;
using tgt::vec4;

namespace voreen {

PlaneWidgetProcessor::PlaneWidgetProcessor()
    : GeometryRendererBase()
    , grabbed_(-1)
    , syncMovement_(false)
    , enable_("enable", "Enable", true)
    , showHandles_("showHandles", "Show Handles", true)
    , manipulatorScaling_("manipulatorScaling", "Manipulator scaling", 1.0f, 0.1f, 10.0f)
    , xColor_("xColor", "x Color", tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f))
    , yColor_("yColor", "y Color", tgt::vec4(0.0f, 1.0f, 0.0f, 1.0f))
    , zColor_("zColor", "z Color", tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f))
    , clipEnabledLeftX_("enableLeftX", "Enable left clipping plane (x)", true)
    , clipUseSphereManipulatorLeftX_("useSphereManipulatorLeftX", "Use sphere manipulator for left clipping plane (x)", false)
    , clipEnabledRightX_("enableRightX", "Enable right clipping plane (x)", true)
    , clipUseSphereManipulatorRightX_("useSphereManipulatorRightX", "Use sphere manipulator for right clipping plane (x)", false)
    , clipEnabledFrontY_("enableBottomY", "Enable front clipping plane (y)", true)
    , clipUseSphereManipulatorFrontY_("useSphereManipulatorFrontY", "Use sphere manipulator for front clipping plane (y)", false)
    , clipEnabledBackY_("enableTopY", "Enable back clipping plane (y)", true)
    , clipUseSphereManipulatorBackY_("useSphereManipulatorBackY", "Use sphere manipulator for back clipping plane (y)", false)
    , clipEnabledBottomZ_("enableFrontZ", "Enable bottom clipping plane (z)", true)
    , clipUseSphereManipulatorBottomZ_("useSphereManipulatorBottomZ", "Use sphere manipulator for bottom clipping plane (z)", false)
    , clipEnabledTopZ_("enableBackZ", "Enable top clipping plane (z)", true)
    , clipUseSphereManipulatorTopZ_("useSphereManipulatorTopZ", "Use sphere manipulator for top clipping plane (z)", false)
    , clipRegion_("clipRegion", "Clip Region")
    , width_("lineWidth", "Line Width", 1.0f, 1.0f, 10.0f)
    , showInnerBB_("showInnerBB", "Show inner box", false)
    , innerColor_("innerColor", "Inner box color", tgt::vec4(0.0f, 0.0f, 0.0f, 1.0f))
    , inport_(Port::INPORT, "volume", "Volume Input")
{
    addPort(inport_);
    inport_.onNewData(MemberFunctionCallback<PlaneWidgetProcessor>(this,&PlaneWidgetProcessor::inportNewData));

    addProperty(enable_);

    addProperty(showHandles_);
    //enable_.setGroupID("vis");

    moveEventProp_ = new EventProperty<PlaneWidgetProcessor>(
        "mouseEvent.clipplaneMovement", "Clipplane movement", this,
        &PlaneWidgetProcessor::planeMovement,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED,
        tgt::Event::MODIFIER_NONE, false);
    addEventProperty(moveEventProp_);
    syncMoveEventProp_ = new EventProperty<PlaneWidgetProcessor>(
        "mouseEvent.synchronizedMovement", "Synchronized movement", this,
        &PlaneWidgetProcessor::planeMovementSync,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED,
        tgt::Event::SHIFT, false);
    addEventProperty(syncMoveEventProp_);

    addProperty(width_);
        //width_.setGroupID("vis");
    addProperty(manipulatorScaling_);
        //manipulatorScaling_.setGroupID("vis");
    //setPropertyGroupGuiName("vis", "Visual Settings");

    addProperty(clipRegion_);

    addProperty(xColor_);
        xColor_.setGroupID("x");
    addProperty(clipEnabledRightX_);
        clipEnabledRightX_.setGroupID("x");
    addProperty(clipUseSphereManipulatorRightX_);
        clipUseSphereManipulatorRightX_.setGroupID("x");
    addProperty(clipEnabledLeftX_);
        clipEnabledLeftX_.setGroupID("x");
    addProperty(clipUseSphereManipulatorLeftX_);
        clipUseSphereManipulatorLeftX_.setGroupID("x");
    setPropertyGroupGuiName("x", "X Axis Clipping");

    addProperty(yColor_);
        yColor_.setGroupID("y");
    addProperty(clipEnabledFrontY_);
        clipEnabledFrontY_.setGroupID("y");
    addProperty(clipUseSphereManipulatorFrontY_);
        clipUseSphereManipulatorFrontY_.setGroupID("y");
    addProperty(clipEnabledBackY_);
        clipEnabledBackY_.setGroupID("y");
    addProperty(clipUseSphereManipulatorBackY_);
        clipUseSphereManipulatorBackY_.setGroupID("y");
    setPropertyGroupGuiName("y", "Y Axis Clipping");

    addProperty(zColor_);
        zColor_.setGroupID("z");
    addProperty(clipEnabledBottomZ_);
        clipEnabledBottomZ_.setGroupID("z");
    addProperty(clipUseSphereManipulatorBottomZ_);
        clipUseSphereManipulatorBottomZ_.setGroupID("z");
    addProperty(clipEnabledTopZ_);
        clipEnabledTopZ_.setGroupID("z");
    addProperty(clipUseSphereManipulatorTopZ_);
        clipUseSphereManipulatorTopZ_.setGroupID("z");
    setPropertyGroupGuiName("z", "Z Axis Clipping");

    addProperty(showInnerBB_);
        showInnerBB_.setGroupID("bb");
    addProperty(innerColor_);
        innerColor_.setGroupID("bb");
    setPropertyGroupGuiName("bb", "Inner Bounding Box");

    //X
    manipulators_.reserve(12);
    manipulators_.push_back(Manipulator(0, 0, &clipRegion_, false));
    manipulators_.push_back(Manipulator(0, 2, &clipRegion_, false));
    manipulators_.push_back(Manipulator(0, 0, &clipRegion_, true));
    manipulators_.push_back(Manipulator(0, 2, &clipRegion_, true));

    //Y
    manipulators_.push_back(Manipulator(1, 0, &clipRegion_, false));
    manipulators_.push_back(Manipulator(1, 2, &clipRegion_, false));
    manipulators_.push_back(Manipulator(1, 0, &clipRegion_, true));
    manipulators_.push_back(Manipulator(1, 2, &clipRegion_, true));

    //Z
    manipulators_.push_back(Manipulator(2, 0, &clipRegion_, false));
    manipulators_.push_back(Manipulator(2, 2, &clipRegion_, false));
    manipulators_.push_back(Manipulator(2, 0, &clipRegion_, true));
    manipulators_.push_back(Manipulator(2, 2, &clipRegion_, true));

    ON_PROPERTY_CHANGE(clipEnabledRightX_,PlaneWidgetProcessor,propertyVisibilityOnChange);
    ON_PROPERTY_CHANGE(clipEnabledLeftX_,PlaneWidgetProcessor,propertyVisibilityOnChange);
    ON_PROPERTY_CHANGE(clipEnabledFrontY_,PlaneWidgetProcessor,propertyVisibilityOnChange);
    ON_PROPERTY_CHANGE(clipEnabledBackY_,PlaneWidgetProcessor,propertyVisibilityOnChange);
    ON_PROPERTY_CHANGE(clipEnabledTopZ_,PlaneWidgetProcessor,propertyVisibilityOnChange);
    ON_PROPERTY_CHANGE(clipEnabledBottomZ_,PlaneWidgetProcessor,propertyVisibilityOnChange);
}

PlaneWidgetProcessor::~PlaneWidgetProcessor() {
    delete moveEventProp_;
    delete syncMoveEventProp_;
}

Processor* PlaneWidgetProcessor::create() const {
    return new PlaneWidgetProcessor();
}

tgt::Bounds PlaneWidgetProcessor::getBoundingBox() const {
    if (inport_.hasData())
        return inport_.getData()->getBoundingBox().getBoundingBox(false);

    return GeometryRendererBase::getBoundingBox();
}

std::string PlaneWidgetProcessor::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = GeometryRendererBase::generateHeader(version);
    header += arrowGeometry_.getShaderDefines();
    return header;
}

void PlaneWidgetProcessor::initialize() {
    GeometryRendererBase::initialize();

    // create static geometry
    sphereGeometry_.setSphereGeometry(0.025f, tgt::vec3::zero, tgt::vec4::one, 25);
    arrowGeometry_.addCylinderGeometry(0.0f, 0.025f, 0.05f, vec3(0.0f,0.0f,0.0f));
    arrowGeometry_.addCylinderGeometry(0.01f, 0.01f, 0.06f, vec3(0.0f,0.0f,0.05f));

    shader_ = ShdrMgr.loadSeparate("trianglemesh.vert", "trianglemesh.geom", "trianglemesh.frag", generateHeader(), false);
}

void PlaneWidgetProcessor::deinitialize() {
    ShdrMgr.dispose(shader_);
    GeometryRendererBase::deinitialize();
}

void PlaneWidgetProcessor::invalidate(int inv) {
    GeometryRendererBase::invalidate(inv);
}

void PlaneWidgetProcessor::inportNewData() {
   if (inport_.hasData()) {
        tgt::ivec3 oldRange = clipRegion_.getMaxValue();
        tgt::ivec3 newRange = tgt::ivec3(inport_.getData()->getDimensions())-tgt::ivec3::one;

        if (oldRange != newRange){
            clipRegion_.setMaxValue(newRange);
            clipRegion_.set(tgt::IntBounds(tgt::ivec3::zero, newRange));
        }
    }
    else {
        clipRegion_.setMaxValue(tgt::ivec3(100000));
    }
}

void PlaneWidgetProcessor::propertyVisibilityOnChange() {
    //x settings
    bool anyClipEnabled =
        clipEnabledLeftX_.get() || clipEnabledRightX_.get() ||
        clipEnabledFrontY_.get() || clipEnabledBackY_.get() ||
        clipEnabledBottomZ_.get() || clipEnabledTopZ_.get();

    clipRegion_.setVisibleFlag(anyClipEnabled);


    clipUseSphereManipulatorRightX_.setVisibleFlag(clipEnabledRightX_.get());
    clipUseSphereManipulatorLeftX_.setVisibleFlag(clipEnabledLeftX_.get());
    //y settings
    clipUseSphereManipulatorFrontY_.setVisibleFlag(clipEnabledFrontY_.get());
    clipUseSphereManipulatorBackY_.setVisibleFlag(clipEnabledBackY_.get());
    //z settings
    clipUseSphereManipulatorTopZ_.setVisibleFlag(clipEnabledTopZ_.get());
    clipUseSphereManipulatorBottomZ_.setVisibleFlag(clipEnabledBottomZ_.get());
}

void PlaneWidgetProcessor::setIDManager(IDManager* idm) {
    if(idm == idManager_)
        return;

    idManager_ = idm;

    if(idManager_) {
        for (size_t i = 0; i < manipulators_.size(); ++i)
            idManager_->registerObject(&manipulators_[i]);
    }
}

void PlaneWidgetProcessor::planeMovement(tgt::MouseEvent* e) {
    e->ignore();
    if (!idManager_)
        return;

    if (e->action() & tgt::MouseEvent::PRESSED) {
        const void* obj = idManager_->getObjectAtPos(tgt::ivec2(e->coord().x, e->viewport().y - e->coord().y));
        for (size_t i = 0; i < manipulators_.size(); ++i) {
            if(obj == &manipulators_[i]) {
                grabbed_ = static_cast<int>(i);
                break;
            }
        }
        if(grabbed_ != -1) {
            e->accept();
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            invalidate();
        }
    }

    if (e->action() & tgt::MouseEvent::MOTION) {
        if(grabbed_ != -1) {
            e->accept();

            MatStack.pushMatrix();
            MatStack.multMatrix(inport_.getData()->getPhysicalToWorldMatrix());

            //adjacent vertices
            tgt::vec3 vertex1 = getCorner(manipulators_[grabbed_].axis_, manipulators_[grabbed_].corner_, 1.0f);
            tgt::vec3 vertex2 = getCorner(manipulators_[grabbed_].axis_, manipulators_[grabbed_].corner_, 0.0f);

            //convert coordinates of both points in windowcoordinates
            tgt::mat4 viewMat = inport_.getData()->getPhysicalToWorldMatrix();
            tgt::vec3 vertex1Projected = getWindowPos(vertex1, viewMat);
            tgt::vec3 vertex2Projected = getWindowPos(vertex2, viewMat);

            //calculate projection of mouseposition to line between both vertexs
            tgt::vec2 mousePos = tgt::vec2(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y));
            tgt::vec3 direction = vertex2Projected-vertex1Projected;
            float t = tgt::dot(mousePos-vertex1Projected.xy(), direction.xy()) / tgt::lengthSq(direction.xy());
            if (t < 0.f)
                t = 0.f;
            if (t > 1.f)
                t = 1.f;

            Manipulator *manipulator = &manipulators_[grabbed_];
            tgt::IntBounds bounds = manipulator->prop_->get();

            tgt::ivec3 vecbounds[2];
            vecbounds[0] = bounds.getURB();
            vecbounds[1] = bounds.getLLF();
            int side = manipulator->inverted_ ? 0 : 1;
            int oppositeSide = 1-side;
            int axis = manipulator->axis_;

            //save difference before move for synced movement:
            float diff = static_cast<float>(vecbounds[side][axis]-vecbounds[oppositeSide][axis]);


            //update property value:
            vecbounds[side][axis] = static_cast<int>((1.0f-t)*manipulator->prop_->getMaxValue()[axis]);

            if(syncMovement_) {
                float temp = vecbounds[side][axis] - diff;
                vecbounds[oppositeSide][axis] = static_cast<int>(temp);
            }

            // we need to make sure that interval has positive volume
            if (side == 0){
                vecbounds[oppositeSide][axis] = std::min(vecbounds[side][axis], vecbounds[oppositeSide][axis]);
            }else{
                vecbounds[oppositeSide][axis] = std::max(vecbounds[side][axis], vecbounds[oppositeSide][axis]);
            }

            bounds = tgt::IntBounds(vecbounds[0], vecbounds[1]);
            manipulator->prop_->set(bounds);

            MatStack.popMatrix();

            invalidate();
        }
    }

    if (e->action() & tgt::MouseEvent::RELEASED) {
        if(grabbed_ != -1) {
            e->accept();
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            invalidate();
            grabbed_ = -1;
        }
        syncMovement_ = false;
    }
}

void PlaneWidgetProcessor::planeMovementSync(tgt::MouseEvent* e) {
    e->ignore();
    if (!idManager_)
        return;

    if (e->action() & tgt::MouseEvent::PRESSED) {
        syncMovement_ = true;
        planeMovement(e);
        if(grabbed_ == -1)
            syncMovement_ = false;
    }
    if (e->action() & tgt::MouseEvent::MOTION) {
        planeMovement(e);
    }
    if (e->action() & tgt::MouseEvent::RELEASED) {
        planeMovement(e);
        syncMovement_ = false;
    }
}

void PlaneWidgetProcessor::renderPicking() {
    if (!enable_.get() || !showHandles_.get())
        return;

    if (!idManager_)
        return;

    MatStack.pushMatrix();
    MatStack.multMatrix(inport_.getData()->getPhysicalToWorldMatrix());

    if(clipEnabledLeftX_.get()) {
        paintManipulator(manipulators_[0], clipUseSphereManipulatorLeftX_.get(), idManager_->getColorFromObjectFloat(&manipulators_[0]), false);
        paintManipulator(manipulators_[1], clipUseSphereManipulatorLeftX_.get(), idManager_->getColorFromObjectFloat(&manipulators_[1]), false);
    }

    if(clipEnabledRightX_.get()) {
        paintManipulator(manipulators_[2], clipUseSphereManipulatorRightX_.get(), idManager_->getColorFromObjectFloat(&manipulators_[2]), false);
        paintManipulator(manipulators_[3], clipUseSphereManipulatorRightX_.get(), idManager_->getColorFromObjectFloat(&manipulators_[3]), false);
    }

    if(clipEnabledFrontY_.get()) {
        paintManipulator(manipulators_[4], clipUseSphereManipulatorFrontY_.get(), idManager_->getColorFromObjectFloat(&manipulators_[4]), false);
        paintManipulator(manipulators_[5], clipUseSphereManipulatorFrontY_.get(), idManager_->getColorFromObjectFloat(&manipulators_[5]), false);
    }

    if(clipEnabledBackY_.get()) {
        paintManipulator(manipulators_[6], clipUseSphereManipulatorBackY_.get(), idManager_->getColorFromObjectFloat(&manipulators_[6]), false);
        paintManipulator(manipulators_[7], clipUseSphereManipulatorBackY_.get(), idManager_->getColorFromObjectFloat(&manipulators_[7]), false);
    }

    if(clipEnabledTopZ_.get()) {
        paintManipulator(manipulators_[8], clipUseSphereManipulatorTopZ_.get(), idManager_->getColorFromObjectFloat(&manipulators_[8]), false);
        paintManipulator(manipulators_[9], clipUseSphereManipulatorTopZ_.get(), idManager_->getColorFromObjectFloat(&manipulators_[9]), false);
    }

    if(clipEnabledBottomZ_.get()) {
        paintManipulator(manipulators_[10], clipUseSphereManipulatorBottomZ_.get(), idManager_->getColorFromObjectFloat(&manipulators_[10]), false);
        paintManipulator(manipulators_[11], clipUseSphereManipulatorBottomZ_.get(), idManager_->getColorFromObjectFloat(&manipulators_[11]), false);
    }

    MatStack.popMatrix();

    LGL_ERROR;
}

void PlaneWidgetProcessor::render() {
    if (!enable_.get())
        return;


    MatStack.pushMatrix();
    MatStack.multMatrix(inport_.getData()->getPhysicalToWorldMatrix());

    glLineWidth(width_.get());

    // llf and urb parameter values are computed for normalized t parameters like using voxelToTexture matrix of the volume to be consistent with bounding boxes and clipping in proxy geometry
    tgt::vec3 llf = inport_.getData()->getVoxelToTextureMatrix() * (tgt::vec3(clipRegion_.get().getLLF()) - tgt::vec3(0.5f));
    tgt::vec3 urb = inport_.getData()->getVoxelToTextureMatrix() * (tgt::vec3(clipRegion_.get().getURB()) + tgt::vec3(0.5f));

    if ( ((grabbed_ == 0) || (grabbed_ == 1))
      || (syncMovement_ && ((grabbed_ == 2) || (grabbed_ == 3)))) {
        float t = llf.x;
        glColor4fv(xColor_.get().elem);
        glBegin(GL_LINE_LOOP);
        glVertex3fv(getCorner(0, 0, t).elem);
        glVertex3fv(getCorner(0, 1, t).elem);
        glVertex3fv(getCorner(0, 2, t).elem);
        glVertex3fv(getCorner(0, 3, t).elem);
        glEnd();
    }
    if ( ((grabbed_ == 2) || (grabbed_ == 3))
      || (syncMovement_ && ((grabbed_ == 0) || (grabbed_ == 1))) ) {
        float t = urb.x;
        glColor4fv(xColor_.get().elem);
        glBegin(GL_LINE_LOOP);
        glVertex3fv(getCorner(0, 0, t).elem);
        glVertex3fv(getCorner(0, 1, t).elem);
        glVertex3fv(getCorner(0, 2, t).elem);
        glVertex3fv(getCorner(0, 3, t).elem);
        glEnd();
    }
    if ( ((grabbed_ == 4) || (grabbed_ == 5))
      || (syncMovement_ && ((grabbed_ == 6) || (grabbed_ == 7))) ){
        float t = llf.y;
        glColor4fv(yColor_.get().elem);
        glBegin(GL_LINE_LOOP);
        glVertex3fv(getCorner(1, 0, t).elem);
        glVertex3fv(getCorner(1, 1, t).elem);
        glVertex3fv(getCorner(1, 2, t).elem);
        glVertex3fv(getCorner(1, 3, t).elem);
        glEnd();
    }
    if ( ((grabbed_ == 6) || (grabbed_ == 7))
      || (syncMovement_ && ((grabbed_ == 4) || (grabbed_ == 5))) ) {
        float t = urb.y;
        glColor4fv(yColor_.get().elem);
        glBegin(GL_LINE_LOOP);
        glVertex3fv(getCorner(1, 0, t).elem);
        glVertex3fv(getCorner(1, 1, t).elem);
        glVertex3fv(getCorner(1, 2, t).elem);
        glVertex3fv(getCorner(1, 3, t).elem);
        glEnd();
    }
    if ( ((grabbed_ == 8) || (grabbed_ == 9))
      || (syncMovement_ && ((grabbed_ == 10) || (grabbed_ == 11))) ) {
        float t = llf.z;
        glColor4fv(zColor_.get().elem);
        glBegin(GL_LINE_LOOP);
        glVertex3fv(getCorner(2, 0, t).elem);
        glVertex3fv(getCorner(2, 1, t).elem);
        glVertex3fv(getCorner(2, 2, t).elem);
        glVertex3fv(getCorner(2, 3, t).elem);
        glEnd();
    }
    if ( ((grabbed_ == 10) || (grabbed_ == 11))
      || (syncMovement_ && ((grabbed_ == 8) || (grabbed_ == 9))) ) {
        float t = urb.z;
        glColor4fv(zColor_.get().elem);
        glBegin(GL_LINE_LOOP);
        glVertex3fv(getCorner(2, 0, t).elem);
        glVertex3fv(getCorner(2, 1, t).elem);
        glVertex3fv(getCorner(2, 2, t).elem);
        glVertex3fv(getCorner(2, 3, t).elem);
        glEnd();
    }

    if (showInnerBB_.get()) {
        //inner boundingbox:
        vec3 min = llf * (getUrb() - getLlf()) + getLlf();
        vec3 max = urb * (getUrb() - getLlf()) + getLlf();

        glColor4fv(innerColor_.get().elem);
        glBegin(GL_LINES);

        //---------------------------

        glVertex3f(min.x, min.y, min.z);
        glVertex3f(min.x, min.y, max.z);

        glVertex3f(min.x, max.y, min.z);
        glVertex3f(min.x, max.y, max.z);

        glVertex3f(max.x, min.y, min.z);
        glVertex3f(max.x, min.y, max.z);

        glVertex3f(max.x, max.y, min.z);
        glVertex3f(max.x, max.y, max.z);

        //---------------------------

        glVertex3f(min.x, min.y, min.z);
        glVertex3f(max.x, min.y, min.z);

        glVertex3f(min.x, max.y, min.z);
        glVertex3f(max.x, max.y, min.z);

        glVertex3f(min.x, min.y, max.z);
        glVertex3f(max.x, min.y, max.z);

        glVertex3f(min.x, max.y, max.z);
        glVertex3f(max.x, max.y, max.z);

        //---------------------------

        glVertex3f(min.x, min.y, min.z);
        glVertex3f(min.x, max.y, min.z);

        glVertex3f(max.x, min.y, min.z);
        glVertex3f(max.x, max.y, min.z);

        glVertex3f(min.x, min.y, max.z);
        glVertex3f(min.x, max.y, max.z);

        glVertex3f(max.x, min.y, max.z);
        glVertex3f(max.x, max.y, max.z);

        //---------------------------

        glEnd();

    }
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glLineWidth(1.0f);

    //---------------------------------------------------------------

    //render manipulators:
    if(showHandles_.get()) {
        //X
        if(clipEnabledLeftX_.get()) {
            paintManipulator(manipulators_[0], clipUseSphereManipulatorLeftX_.get(), xColor_.get(), true);
            paintManipulator(manipulators_[1], clipUseSphereManipulatorLeftX_.get(), xColor_.get(), true);
        }
        if(clipEnabledRightX_.get()) {
            paintManipulator(manipulators_[2], clipUseSphereManipulatorRightX_.get(), xColor_.get(), true);
            paintManipulator(manipulators_[3], clipUseSphereManipulatorRightX_.get(), xColor_.get(), true);
        }

        //Y
        if(clipEnabledFrontY_.get()) {
            paintManipulator(manipulators_[4], clipUseSphereManipulatorFrontY_.get(), yColor_.get(), true);
            paintManipulator(manipulators_[5], clipUseSphereManipulatorFrontY_.get(), yColor_.get(), true);
        }
        if(clipEnabledBackY_.get()) {
            paintManipulator(manipulators_[6], clipUseSphereManipulatorBackY_.get(), yColor_.get(), true);
            paintManipulator(manipulators_[7], clipUseSphereManipulatorBackY_.get(), yColor_.get(), true);
        }

        //Z
        if(clipEnabledTopZ_.get()) {
            paintManipulator(manipulators_[8], clipUseSphereManipulatorTopZ_.get(), zColor_.get(), true);
            paintManipulator(manipulators_[9], clipUseSphereManipulatorTopZ_.get(), zColor_.get(), true);
        }
        if(clipEnabledBottomZ_.get()) {
            paintManipulator(manipulators_[10], clipUseSphereManipulatorBottomZ_.get(), zColor_.get(), true);
            paintManipulator(manipulators_[11], clipUseSphereManipulatorBottomZ_.get(), zColor_.get(), true);
        }
    }

    MatStack.popMatrix();

    LGL_ERROR;
}

void PlaneWidgetProcessor::paintManipulator(const Manipulator& a, bool useSphereManipulator, vec4 color, bool lighting) {
    vec3 rot;
    float angle = 0.f;
    switch(a.axis_) {
        case 0: rot = vec3(0.0f, 1.0f, 0.0f);
                angle = -90.0f;
                if(a.inverted_)
                    angle = 90.0f;
                break;
        case 1: rot = vec3(1.0f, 0.0f, 0.0f);
                angle = 90.0f;
                if(a.inverted_)
                    angle = -90.0f;
                break;
        case 2: rot = vec3(1.0f, 1.0f, 0.0f);
                angle = 0.0f;
                if(!a.inverted_)
                    angle = 180.0f;
                break;
    }

    float scalingFactor = 0.3f * tgt::length(getUrb() - getLlf()) * manipulatorScaling_.get();
    int axis = a.axis_;
    float cornerPos;
    // cornerPos is computed for normalized t parameter like using voxelToTexture matrix of the volume to be consistent with bounding boxes and clipping in proxy geometry
    if (a.inverted_){
        cornerPos = static_cast<float>(a.prop_->get().getURB()[axis]) + 0.5f; // voxel space
    }
    else{
        cornerPos = static_cast<float>(a.prop_->get().getLLF()[axis]) - 0.5f; // voxel space
    }
    // convert from voxel space to normalized texture space t parameter
    cornerPos += 0.5f;
    cornerPos /= (static_cast<float>(a.prop_->getMaxValue()[axis]) + 1.f);

    paintManipulator(getCorner(a.axis_, a.corner_, cornerPos), rot, angle, scalingFactor, useSphereManipulator, color, lighting);
}

void PlaneWidgetProcessor::paintManipulator(tgt::vec3 translation, tgt::vec3 rotationAxis, float rotationAngle, float scaleFactor, bool useSphereManipulator, vec4 color, bool lighting)
{
    MatStack.pushMatrix();
    MatStack.translate(translation.x, translation.y, translation.z);
    MatStack.rotate(rotationAngle, rotationAxis.x, rotationAxis.y, rotationAxis.z);
    MatStack.scale(scaleFactor, scaleFactor, scaleFactor);

    // select the appropriate geometry
    //TriangleMeshGeometryUInt16IndexedNormal* geometry = NULL;
    Geometry* geometry = NULL;
    if(useSphereManipulator) {
        geometry = &sphereGeometry_;
    }
    else {
        geometry = &arrowGeometry_;
    }

    shader_->activate();

    shader_->setUniform("solidColor_", color);

    // set transformation matrices
    shader_->setUniform("modelMatrix_", tgt::mat4::identity);
    shader_->setUniform("viewMatrix_", MatStack.getModelViewMatrix());
    shader_->setUniform("projectionMatrix_", MatStack.getProjectionMatrix());
    tgt::mat3 modelViewMat3 = MatStack.getModelViewMatrix().getRotationalPartMat3();
    tgt::mat3 normalMatrix;
    if(!modelViewMat3.invert(normalMatrix)) {
        LWARNING("Could not generate normal matrix out of current view / model matrix, using identity.");
        normalMatrix = tgt::mat3::identity;
    }
    normalMatrix = transpose(normalMatrix);
    shader_->setUniform("normalMatrix_", normalMatrix);


    // set lighting uniforms
    shader_->setUniform("lightSource_.ambientColor_", vec3(0.5f, 0.5f, 0.5f));
    shader_->setUniform("lightSource_.diffuseColor_", vec3(1.0f, 1.0f, 1.0f));
    shader_->setUniform("lightSource_.specularColor_", vec3(1.0f, 1.0f, 1.0f));
    shader_->setUniform("shininess_", 25.0f);
    shader_->setUniform("lightPositionEye_", vec3(0.f, 0.f, 1.f));
    shader_->setUniform("enableLighting_", lighting);

    // disable clipping
    shader_->setUniform("enableClipping_", false);
    shader_->setUniform("plane_", tgt::vec4(0.0f, 0.0f, 0.0f, 0.0f));

    // draw the generated geometry
    geometry->render();

    shader_->deactivate();
    MatStack.popMatrix();
}

tgt::vec3 PlaneWidgetProcessor::getLlf() {
    return inport_.getData()->getLLF();
}

tgt::vec3 PlaneWidgetProcessor::getUrb() {
    return inport_.getData()->getURB();
}

vec3 PlaneWidgetProcessor::getCorner(int axis, int num, float t) {
    if(axis == 0) {
        float xSlice = (t * (getUrb().x - getLlf().x)) + getLlf().x;
        switch(num) {
            case 0:
                return vec3(xSlice, getUrb().y, getUrb().z);
            case 1:
                return vec3(xSlice, getLlf().y, getUrb().z);
            case 2:
                return vec3(xSlice, getLlf().y, getLlf().z);
            case 3:
                return vec3(xSlice, getUrb().y, getLlf().z);
        }
    }
    else if(axis == 1) {
        float ySlice = (t * (getUrb().y - getLlf().y)) + getLlf().y;
        switch(num) {
            case 0:
                return vec3(getLlf().x, ySlice, getLlf().z);
            case 1:
                return vec3(getLlf().x, ySlice, getUrb().z);
            case 2:
                return vec3(getUrb().x, ySlice, getUrb().z);
            case 3:
                return vec3(getUrb().x, ySlice, getLlf().z);
        }
    }
    else if(axis == 2) {
        float zSlice = (t * (getUrb().z - getLlf().z)) + getLlf().z;
        switch(num) {
            case 0:
                return vec3(getLlf().x, getUrb().y, zSlice);
            case 1:
                return vec3(getLlf().x, getLlf().y, zSlice);
            case 2:
                return vec3(getUrb().x, getLlf().y, zSlice);
            case 3:
                return vec3(getUrb().x, getUrb().y, zSlice);
        }
    }
    return vec3(0.0f);
}


} //namespace voreen
