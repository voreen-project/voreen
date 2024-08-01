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

#include "regionofinterest2d.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"

using tgt::TextureUnit;

namespace voreen {

using tgt::ivec2;
using tgt::mat4;
using tgt::MouseEvent;
using tgt::vertex;
using tgt::vec2;
using tgt::vec3;
using tgt::vec4;
using std::vector;

//---------------------------------------------------------------------------

namespace {

// GLU tesselation callbacks

#ifndef CALLBACK
#define CALLBACK
#endif

    void CALLBACK tessBegin(GLenum primitive) {
        glBegin(primitive);
    }

    void CALLBACK tessEnd() {
        glEnd();
    }

    void CALLBACK tessVertex2dv(GLdouble* data) {
        glVertex2dv(data);
    }

    void CALLBACK tessCombine(GLdouble coords[3],
        GLdouble* /*vertex_data[4]*/, GLfloat /*weight[4]*/,
        GLdouble** dataOut)
    {
        GLdouble* vertex = new GLdouble[3];
        vertex[0] = coords[0];
        vertex[1] = coords[1];
        vertex[2] = coords[2];

        *dataOut = vertex;
    }
}

//---------------------------------------------------------------------------


const std::string RegionOfInterest2D::RegionOfInterestGeometry::loggerCat_("voreen.RegionOfInterest2D.RegionOfInterestGeometry");

RegionOfInterest2D::RegionOfInterestGeometry::RegionOfInterestGeometry()
    : Geometry()
    , layer_(0)
    , transformMatrix_(mat4::createIdentity())
    , geometryMode_(LINE_STRIP)
    , roiColor_(tgt::Color(1.f, 1.f, 1.f, 1.f))
    , maskColor_(tgt::Color(1.f, 1.f, 1.f, 1.f))
    , size_(1.f)
    , antialiasing_(false)
{}

RegionOfInterest2D::RegionOfInterestGeometry::~RegionOfInterestGeometry() {
}

void RegionOfInterest2D::RegionOfInterestGeometry::renderBoundingBox() const {
    // apply roi transformation matrix
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(transformMatrix_);

    vec2 ll = vec2(boundingBox_.x, boundingBox_.y);
    vec2 ur = vec2(boundingBox_.z, boundingBox_.w);
    vec2 lr = vec2(ur.x, ll.y);
    vec2 ul = vec2(ll.x, ur.y);

    IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(ll);
        IMode.vertex(lr);
        IMode.vertex(ur);
        IMode.vertex(ul);
    IMode.end();

    MatStack.popMatrix();

    LGL_ERROR;
}

void RegionOfInterest2D::RegionOfInterestGeometry::computeBoundingBox() {
    // compute bounding box
    vec2 ll = vec2(FLT_MAX);
    vec2 ur = vec2(-1.f);
    for (size_t i=0; i < points_.size(); ++i) {
        ll = tgt::min(ll, points_.at(i));
        ur = tgt::max(ur, points_.at(i));
    }
    boundingBox_ = vec4(ll.x, ll.y, ur.x, ur.y);
}

void RegionOfInterest2D::RegionOfInterestGeometry::fillArea(tgt::Texture* mask) {
    fillPoints_.clear();
    mask->downloadTexture();
    for(int i = 0; i< mask->getWidth(); i++) {
        for(int j = 0; j < mask->getHeight(); j++) {
            if(mask->texelAsFloat(i, j).r > 0.f)
                fillPoints_.push_back(tgt::vec2((float)i, (float)j));
        }
    }
}

void RegionOfInterest2D::RegionOfInterestGeometry::render() const {
    // apply roi transformation matrix
    MatStack.pushMatrix();
    MatStack.multMatrix(transformMatrix_);

    // Use OpenGL anti-aliasing?
    if (antialiasing_) {
        //glEnable(GL_POINT_SMOOTH); //TODO: reenable via shaders? probably not worth it
        glEnable(GL_LINE_SMOOTH);
    }
    else {
        //glDisable(GL_POINT_SMOOTH);
        glDisable(GL_LINE_SMOOTH);
    }

    LGL_ERROR;

    RegionOfInterest2D::RegionOfInterestGeometry::GeometryMode mode = geometryMode_;
    if (mode == POINTS || mode == LINES || mode == LINE_STRIP || mode == LINE_LOOP || mode == RECTANGLE) {
        // determine GL geometry mode
        tgt::ImmediateMode::PolygonMode glMode;
        switch (mode) {
        case POINTS:
            glMode = tgt::ImmediateMode::POINTS;
            break;
        case LINES:
            glMode = tgt::ImmediateMode::LINES;
            break;
        case LINE_STRIP:
            glMode = tgt::ImmediateMode::LINE_STRIP;
            break;
        case LINE_LOOP:
            glMode = tgt::ImmediateMode::LINE_LOOP;
            break;
        case RECTANGLE:
            glMode = tgt::ImmediateMode::QUADS;
            break;
        default:
            glMode = tgt::ImmediateMode::POINTS;
            break;
        }

        // render curve with specified line / point size
        if (geometryMode_ == POINTS)
            glPointSize(size_);
        else
            glLineWidth(size_);

        IMode.begin(glMode);
        for (size_t i = 0; i < points_.size(); ++i)
            IMode.vertex(points_[i]);
        IMode.end();
        LGL_ERROR;

    }
    //not supported in core profile. Shaders are not working with glu.. command
    /*else if (mode == SPLINE) {
        // transform ROI points to spline points
        vector<vec2> roiPoints = points_;
        vector<vec3> splinePoints;
        float layer = static_cast<float>(layer_);
        vec3 last = vec3(roiPoints.front(), layer);
        last.z = 0.f;
        splinePoints.push_back(last);
        splinePoints.push_back(vec3(last));
        for (size_t i = 1; i < points_.size(); ++i) {
            vec3 cur = vec3(roiPoints.at(i), layer);
            cur.z = 0.f;
            if (length(cur - last) > 0.f) {
                splinePoints.push_back(cur);
                last = cur;
            }
        }
        last = vec3(roiPoints.back(), layer);
        last.z = 0.f;
        splinePoints.push_back(last);
        splinePoints.push_back(last);

        // Create glu nurbs renderer
        GLUnurbsObj* theNurb = gluNewNurbsRenderer();

        // Create NURBS knots and control points
        const int order = 4;
        int nknots = static_cast<int>(splinePoints.size()) + order;
        GLfloat* knots = new GLfloat[nknots];
        for (int i = 0; i < nknots; ++i)
            knots[i] = GLfloat(i);
        GLfloat* ctlarray = new GLfloat[(nknots-order)*3];
        for (int i = 0; i < nknots-order; ++i) {
            ctlarray[i*3] = splinePoints.at(i).x;
            ctlarray[i*3+1] = splinePoints.at(i).y;
            ctlarray[i*3+2] = splinePoints.at(i).z;
        }

        // render NURBS curve with specified size
        glLineWidth(size_);
        gluBeginCurve(theNurb);
        gluNurbsCurve(theNurb, nknots, knots, 3, ctlarray, 4, GL_MAP1_VERTEX_3);
        gluEndCurve(theNurb);

        delete[] knots;
        delete[] ctlarray;
        gluDeleteNurbsRenderer(theNurb);

        LGL_ERROR;
    }

    else if (mode == PATCH || mode == CIRCLE_PATCH) {
        // initialize GLU tesselator
        GLUtesselator* tobj = gluNewTess();
        gluTessCallback(tobj, GLU_TESS_BEGIN, (void(CALLBACK*)() ) &tessBegin);
        gluTessCallback(tobj, GLU_TESS_END, (void(CALLBACK*)() ) &tessEnd);
        gluTessCallback(tobj, GLU_TESS_VERTEX, (void(CALLBACK*)() ) &tessVertex2dv);
        gluTessCallback(tobj, GLU_TESS_COMBINE, (void(CALLBACK*)() ) &tessCombine);
        gluTessProperty(tobj, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO);

        // tesselate freeform patch and write result to OpenGL feedback buffer
        gluTessBeginPolygon(tobj, NULL);
        gluTessBeginContour(tobj);
        for (size_t i=0; i < points_.size(); ++i) {
            GLdouble* vertex = new GLdouble[3];
            vertex[0] = static_cast<double>(points_[i].x);
            vertex[1] = static_cast<double>(points_[i].y);
            vertex[2] = 0.0;
            gluTessVertex(tobj, vertex, vertex);
        }
        gluTessEndContour(tobj);
        gluTessEndPolygon(tobj);
        gluDeleteTess(tobj);

        LGL_ERROR;
    }*/
    else
        LERROR("Unknown ROI geometry mode. (PATCH, CIRCLE-PATCH and SPLINE are deprecated)");

    // revert roi transformation
    MatStack.popMatrix();

    // Restore gl state
    //glDisable(GL_POINT_SMOOTH);
    glDisable(GL_LINE_SMOOTH);
    glLineWidth(1.0f);
    glPointSize(1.0f);
    LGL_ERROR;
}

void RegionOfInterest2D::RegionOfInterestGeometry::serialize(Serializer& s) const {
    s.serialize("points", points_);
    s.serialize("layer", layer_);
    s.serialize("viewport", viewportSize_);
    s.serialize("trafoMatrix", transformMatrix_);
    s.serialize("geometryMode", geometryMode_);
    s.serialize("roiColor", roiColor_);
    s.serialize("maskColor", maskColor_);
    s.serialize("size", size_);
    s.serialize("antialiasing", antialiasing_);
}

void RegionOfInterest2D::RegionOfInterestGeometry::deserialize(Deserializer& s) {

    try {
        s.deserialize("points", points_);
        s.deserialize("viewport", viewportSize_);
        int mode;
        s.deserialize("geometryMode", mode);
        geometryMode_ = static_cast<RegionOfInterestGeometry::GeometryMode>(mode);
        s.deserialize("roiColor", roiColor_);
        s.deserialize("maskColor", maskColor_);
        s.deserialize("size", size_);
        s.deserialize("antialiasing", antialiasing_);
        s.deserialize("layer", layer_);
        s.deserialize("trafoMatrix", transformMatrix_);
    }
    catch (SerializationNoSuchDataException&) {
        // superfluous or missing XML data item, just ignore
        s.removeLastError();
    }

    computeBoundingBox();
}

Geometry* RegionOfInterest2D::RegionOfInterestGeometry::create() const {
    return new RegionOfInterestGeometry();
}

//---------------------------------------------------------------------------

RegionOfInterest2D::RegionOfInterest2D()
    : ImageProcessor("copyimage")
    , inport_(Port::INPORT, "inport", "Image Input")
    , outportMask_(Port::OUTPORT, "image.mask", "Mask Image Output")
    , maskPrivate_(Port::OUTPORT, "image.mask_private")
    , outportROI_(Port::OUTPORT, "image.roi", "ROI Image Output")
    , outportOverlay_(Port::OUTPORT, "image.overlay", "ROI Image Overlay Output")
    , pickingPort_(Port::OUTPORT, "pickingTarget")
    , outportROIGeometryLayer_(Port::OUTPORT, "geometry.roiGeometryLayer", "ROI Geometry Layer")
    , outportROIGeometry_(Port::OUTPORT, "geometry.roiGeometry", "ROI Geometry Output")
    , geometryMode_("geometryMode", "Geometry")
    , roiColor_("roiColor", "ROI Color", tgt::Color(1.f, 1.f, 0.0f, 0.5f))
    , maskColor_("maskColor", "Mask Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , maskBackgroundColor_("maskBackgroundColor", "Mask Background", tgt::Color(0.f, 0.f, 0.f, 0.f))
    , size_("size", "Line/Point size", 3.f, 1.f, 35.f)
    , antialiasing_("antialiasing", "Antialiasing", false)
    , currentLayer_("roiLayer", "ROI Layer", 0, 0, 10000)
    , renderCurrentLayerOnly_("renderCurrentLayerOnly", "Render current layer only", false)
    , boundingBoxColor_("boundingBoxColor", "Bounding Box Color", tgt::Color(0.8f, 0.8f, 0.8f, 0.8f))
    , minRoiPointDistance_("minRoiPointDistance", "Min ROI Point Distance", 2.f, 0.f, 20.f)
    , fillRegions_("fillRegions", "Inner points are ROI points", false)
    , roiTransformMatrix_("transformationMatrix", "Transformation Matrix", tgt::mat4::createIdentity(), tgt::mat4(-1e6), tgt::mat4(1e6))
    , loadROIs_("loadROIs", "Load ROIs", "Select ROI file", VoreenApplication::app()->getUserDataPath(),
        "Voreen ROI files (*.vrf)", FileDialogProperty::OPEN_FILE, Processor::VALID)
    , saveROIs_("saveROIs", "Save ROIs", "Select ROI file", VoreenApplication::app()->getUserDataPath(),
        "Voreen ROI files (*.vrf)", FileDialogProperty::SAVE_FILE, Processor::VALID)
    , clearROIs_("clearROIs", "Clear ROIs")
    , removeLastROI_("deleteLastRoi", "Remove Last ROI")
    , removeLastROIInCurSlice_("deleteLastRoiInSlice", "Remove Last ROI in Slice")
    , numberOfROIs_("numberOfRois", "Number of ROIs", 0, 0, 100000, Processor::VALID)
    , mouseEventAddROI_("mouseEvent.addROI", "Add ROI", this, &RegionOfInterest2D::addROIEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        static_cast<tgt::MouseEvent::MouseAction>(tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED),
        tgt::Event::CTRL)
    , mouseEventRemoveLastROI_("mouseEvent.removeLastROI", "Remove Last ROI", this, &RegionOfInterest2D::removeLastROI,
        tgt::MouseEvent::MOUSE_BUTTON_RIGHT, tgt::MouseEvent::PRESSED, tgt::Event::CTRL)
    , mouseEventRemoveLastROIInCurSlice_("mouseEvent.removeLastROIInCurSlice", "Remove Last ROI in Cur Slice", this,
        &RegionOfInterest2D::removeROIInCurrentSlice, tgt::MouseEvent::MOUSE_BUTTON_MIDDLE, tgt::MouseEvent::PRESSED, tgt::Event::CTRL)
    , mouseEventRemoveClosestROIInCurSlice_("mouseEvent.removeClosestROIInCurSlice_", "Remove closest ROI in Cur Slice", this,
        &RegionOfInterest2D::removeClosestROIInCurrentSlice, tgt::MouseEvent::MOUSE_BUTTON_4, tgt::MouseEvent::PRESSED, tgt::Event::CTRL)
    , mouseEventClearROIs_("mouseEvent.clearROIs", "Clear ROIs", this, &RegionOfInterest2D::clearROIs,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::DOUBLECLICK, tgt::Event::CTRL)
    , mouseEventShiftROI_("mouseEvent.shiftROI", "Shift ROI", this, &RegionOfInterest2D::shiftROI,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        static_cast<tgt::MouseEvent::MouseAction>(tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED),
        tgt::Event::ALT)
    , mouseEventScaleROI_("mouseEvent.scaleROI", "Scale ROI", this, &RegionOfInterest2D::scaleROI,
        tgt::MouseEvent::MOUSE_BUTTON_RIGHT,
        static_cast<tgt::MouseEvent::MouseAction>(tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED),
        tgt::Event::ALT)
    , blendShader_(0)
    , selectedROI_(-1)
    , lastMousePos_(-1)
    , geometryHasChanged_(true)
{
    addPort(inport_);
    addPort(outportMask_);
    addPort(outportROI_);
    addPort(outportOverlay_);
    addPrivateRenderPort(pickingPort_);
    addPort(outportROIGeometryLayer_);
    addPort(outportROIGeometry_);

    addPrivateRenderPort(maskPrivate_);

    geometryMode_.addOption("points",    "Points",          RegionOfInterestGeometry::POINTS);
    geometryMode_.addOption("lines",     "Lines",           RegionOfInterestGeometry::LINES);
    geometryMode_.addOption("lineStrip", "Line-Strip",      RegionOfInterestGeometry::LINE_STRIP);
    geometryMode_.addOption("lineLoop",  "Line-Loop",       RegionOfInterestGeometry::LINE_LOOP);
    //geometryMode_.addOption("spline",    "Spline",          RegionOfInterestGeometry::SPLINE);
    geometryMode_.addOption("rectangle", "Rectangle",       RegionOfInterestGeometry::RECTANGLE);
    //geometryMode_.addOption("circle",    "Circle",          RegionOfInterestGeometry::CIRCLE_PATCH);
    //geometryMode_.addOption("patch",     "Free-form patch", RegionOfInterestGeometry::PATCH);
    geometryMode_.set("points");

    removeLastROI_.onChange(MemberFunctionOneParameterCallback<RegionOfInterest2D, tgt::MouseEvent*>(this, &RegionOfInterest2D::removeLastROI, 0));
    removeLastROIInCurSlice_.onChange(MemberFunctionOneParameterCallback<RegionOfInterest2D, tgt::MouseEvent*>(this,
                                      &RegionOfInterest2D::removeROIInCurrentSlice, 0));
    clearROIs_.onChange(MemberFunctionOneParameterCallback<RegionOfInterest2D, tgt::MouseEvent*>(this, &RegionOfInterest2D::clearROIs, 0));
    saveROIs_.onChange(MemberFunctionCallback<RegionOfInterest2D>(this, &RegionOfInterest2D::saveFileIssued));
    loadROIs_.onChange(MemberFunctionCallback<RegionOfInterest2D>(this, &RegionOfInterest2D::loadFileIssued));

    mouseEventShiftROI_.setEnabled(false);
    mouseEventScaleROI_.setEnabled(false);

    minRoiPointDistance_.setStepping(1.f);
    numberOfROIs_.setReadOnlyFlag(true);

    addProperty(geometryMode_);
    addProperty(roiColor_);
    addProperty(maskColor_);
    addProperty(maskBackgroundColor_);
    addProperty(size_);
    addProperty(antialiasing_);
    addProperty(currentLayer_);
    addProperty(renderCurrentLayerOnly_);
    addProperty(boundingBoxColor_);
    addProperty(minRoiPointDistance_);
    addProperty(fillRegions_);
    addProperty(roiTransformMatrix_);
    addProperty(loadROIs_);
    addProperty(saveROIs_);
    addProperty(clearROIs_);
    addProperty(removeLastROI_);
    addProperty(removeLastROIInCurSlice_);
    addProperty(numberOfROIs_);

    addEventProperty(mouseEventAddROI_);
    addEventProperty(mouseEventRemoveLastROI_);
    addEventProperty(mouseEventRemoveLastROIInCurSlice_);
    addEventProperty(mouseEventRemoveClosestROIInCurSlice_);
    addEventProperty(mouseEventClearROIs_);
    addEventProperty(mouseEventShiftROI_);
    addEventProperty(mouseEventScaleROI_);
}

bool RegionOfInterest2D::isReady() const {
    bool ready = false;
    ready |= outportMask_.isConnected();
    ready |= outportOverlay_.isConnected();
    ready |= outportROI_.isConnected();
    ready |= outportROIGeometryLayer_.isConnected();
    ready |= outportROIGeometry_.isConnected();
    ready &= inport_.isReady();
    return ready;
}

void RegionOfInterest2D::initialize() {
    ImageProcessor::initialize();

    blendShader_ = ShdrMgr.loadSeparate("immediatemode.vert", "blendwithimage.frag",
        generateHeader(), false);

    idManager_.setRenderTarget(pickingPort_.getRenderTarget());
    idManager_.initializeTarget();

    numberOfROIs_.set(static_cast<int>(regionsOfInterest_.size()));
}

void RegionOfInterest2D::deinitialize() {
    ShdrMgr.dispose(blendShader_);
    blendShader_ = 0;

    ImageProcessor::deinitialize();
}

void RegionOfInterest2D::process() {
    // mask
    if (outportMask_.isReady()) {
        // clear output
        outportMask_.activateTarget();
        tgt::Color bgCol = maskBackgroundColor_.get();
        glClearColor(bgCol.r, bgCol.g, bgCol.b, bgCol.a);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.f, 0.f, 0.f, 0.f);

        // render ROIs with mask color
        renderROIs(MASK_MODE);

        outportMask_.deactivateTarget();
    }

    LGL_ERROR;

    // rendered ROI
    if (outportROI_.isReady()) {
        // clear output buffer
        outportROI_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // render ROIs with roi color, including alpha
        renderROIs(COLOR_MODE);

        // bounding box of selected roi
        if (selectedROI_ > -1) {
            tgtAssert(selectedROI_ >= 0 && selectedROI_ < (int)regionsOfInterest_.size(), "Invalid selected roi");
            if (!renderCurrentLayerOnly_.get() || (regionsOfInterest_.at(static_cast<size_t>(selectedROI_))->layer_ == currentLayer_.get()) )
                renderBoundingBox(regionsOfInterest_.at(static_cast<size_t>(selectedROI_)));
            LGL_ERROR;
        }

        outportROI_.deactivateTarget();
    }

    LGL_ERROR;

    // image with ROI overlay
    if (outportOverlay_.isReady()) {
        outportOverlay_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // copy input image to output buffer
        TextureUnit colorUnit, depthUnit;
        inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());
        program_->activate();
        setGlobalShaderParameters(program_);
        program_->setUniform("colorTex_", colorUnit.getUnitNumber());
        program_->setUniform("depthTex_", depthUnit.getUnitNumber());
        inport_.setTextureParameters(program_, "texParams_");
        renderQuad();
        program_->deactivate();
        LGL_ERROR;
        // render ROIs with ROI color and blend them over the input image
        inport_.bindColorTexture(colorUnit.getEnum());
        blendShader_->activate();
        //setGlobalShaderParameters(blendShader_); done in IMode.end()
        blendShader_->setUniform("colorTex_", colorUnit.getUnitNumber());
        inport_.setTextureParameters(blendShader_, "texParams_");
        renderROIs(COLOR_MODE);
        blendShader_->deactivate();
        LGL_ERROR;
        // bounding box of selected roi
        if (selectedROI_ > -1) {
            tgtAssert(selectedROI_ >= 0 && selectedROI_ < (int)regionsOfInterest_.size(), "Invalid selected roi");
            if (!renderCurrentLayerOnly_.get() || (regionsOfInterest_.at(static_cast<size_t>(selectedROI_))->layer_ == currentLayer_.get()) )
                renderBoundingBox(regionsOfInterest_.at(static_cast<size_t>(selectedROI_)));
        }

        outportOverlay_.deactivateTarget();
    }

    // picking
    renderROIs(PICKING_MODE);

    updateGeometry();

    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void RegionOfInterest2D::renderROIs(RoiMode roiMode) {
    if (!regionsOfInterest_.empty()) {
        glDepthFunc(GL_ALWAYS);

        // set up transformation matrices for specifying
        // ROI points in screen coordinates
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.pushMatrix();
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();

        // viewport transformation
        MatStack.loadIdentity();
        MatStack.translate(-1.f, -1.f, 0.f);
        ivec2 viewportDim = inport_.getSize();
        MatStack.scale(2.f/viewportDim.x , 2.f / viewportDim.y, 1.f);

        // inverse transformation matrix taking us back to screen coords
        mat4 transformMatrixInverse;
        roiTransformMatrix_.get().invert(transformMatrixInverse);
        MatStack.multMatrix(transformMatrixInverse);
        LGL_ERROR;

        if (roiMode == PICKING_MODE) {
            idManager_.clearRegisteredObjects();
            idManager_.activateTarget(getID());
            idManager_.clearTarget();
            LGL_ERROR;
        }

        //
        // Iterate over all ROIs
        //
        for (size_t k=0; k<regionsOfInterest_.size(); k++) {
            if (renderCurrentLayerOnly_.get() && (regionsOfInterest_[k]->layer_ != currentLayer_.get()))
                continue;

            if (regionsOfInterest_[k]->points_.empty())
                continue;

            if (roiMode == MASK_MODE)
                IMode.color(regionsOfInterest_[k]->maskColor_);
            else if (roiMode == COLOR_MODE)
                IMode.color(regionsOfInterest_[k]->roiColor_);
            else if (roiMode == PICKING_MODE) {
                idManager_.registerObject(regionsOfInterest_[k]);
                idManager_.setGLColor(regionsOfInterest_[k]);
                LGL_ERROR;
            }

            regionsOfInterest_[k]->render();
            LGL_ERROR;
        }

        if (roiMode == PICKING_MODE) {
            idManager_.deactivateTarget();
            LGL_ERROR;
        }

        // restore attributes and matrices
        glDepthFunc(GL_LESS);
        glColor4f(1.f, 1.f, 1.f, 1.f);
        IMode.color(tgt::vec4::one);

        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.popMatrix();

        LGL_ERROR;
    }
}

void RegionOfInterest2D::renderROIMask(RegionOfInterestGeometry* roi) {
    if (roi->points_.empty())
        return;

    glDepthFunc(GL_ALWAYS);

    // set up transformation matrices for specifying
    // ROI points in screen coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();

    // viewport transformation
    MatStack.loadIdentity();
    MatStack.translate(-1.f, -1.f, 0.f);
    ivec2 viewportDim = inport_.getSize();
    MatStack.scale(2.f/viewportDim.x , 2.f / viewportDim.y, 1.f);

    // inverse transformation matrix taking us back to screen coords
    mat4 transformMatrixInverse;
    roiTransformMatrix_.get().invert(transformMatrixInverse);
    MatStack.multMatrix(transformMatrixInverse);
    LGL_ERROR;

    IMode.color(roi->maskColor_);

    roi->render();
    LGL_ERROR;

    // restore attributes and matrices
    glDepthFunc(GL_LESS);
    glColor4f(1.f, 1.f, 1.f, 1.f);
    IMode.color(tgt::vec4::one);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    LGL_ERROR;
}

void RegionOfInterest2D::renderBoundingBox(RegionOfInterestGeometry* roi) const {
    glDepthFunc(GL_ALWAYS);

    // set up transformation matrices for specifying
    // ROI points in screen coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();

    // viewport transformation
    MatStack.loadIdentity();
    MatStack.translate(-1.f, -1.f, 0.f);
    ivec2 viewportDim = inport_.getSize();
    MatStack.scale(2.f/viewportDim.x , 2.f / viewportDim.y, 1.f);

    // inverse transformation matrix taking us back to screen coords
    mat4 transformMatrixInverse;
    roiTransformMatrix_.get().invert(transformMatrixInverse);
    MatStack.multMatrix(transformMatrixInverse);

    // render bounding box
    IMode.color(boundingBoxColor_.get());
    roi->renderBoundingBox();

    // restore attributes and matrices
    glDepthFunc(GL_LESS);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    IMode.color(tgt::vec4::one);

    LGL_ERROR;
}

void RegionOfInterest2D::addROIEvent(tgt::MouseEvent* e) {
    tgtAssert(e, "No event");

    if (!inport_.isReady())
        return;

    tgt::ivec2 mousePos(e->x(), inport_.getSize().y - e->y());

    enum Actions {
        NO_ACTION = 0,
        ADD_ROI = 1,
        ADD_ROI_POINT = 2,
        FINISH_ROI = 4
    };

    int actions = NO_ACTION;

    //
    // determine actions
    //
    if (e->action() == MouseEvent::PRESSED) {
        actions |= ADD_ROI;
        actions |= ADD_ROI_POINT;
    }
    if (e->action() == MouseEvent::RELEASED) {
        actions |= FINISH_ROI;
    }
    if (e->action() == MouseEvent::MOTION) {
        actions |= ADD_ROI_POINT;
    }

    //
    // perform actions
    //
    if (actions & ADD_ROI) {
        RegionOfInterestGeometry* roi = new RegionOfInterestGeometry();
        roi->geometryMode_ = geometryMode_.getValue();
        roi->layer_ = currentLayer_.get();
        roi->viewportSize_ = inport_.getSize();
        roi->transformMatrix_ = roiTransformMatrix_.get();
        roi->roiColor_ = roiColor_.get();
        roi->maskColor_ = maskColor_.get();
        roi->size_ = size_.get();
        roi->antialiasing_ = antialiasing_.get();

        // render patch contour during interaction
        /*if (roi->geometryMode_ == RegionOfInterestGeometry::PATCH)
            roi->geometryMode_ = RegionOfInterestGeometry::LINE_STRIP;*/

        regionsOfInterest_.push_back(roi);

        lastMousePos_ = mousePos;
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
    }

    if (actions & ADD_ROI_POINT) {
        if (!regionsOfInterest_.empty()) {
            vec2 point((float)e->coord().x, (float)e->viewport().y - e->coord().y);
            /*if (regionsOfInterest_.back()->geometryMode_ == RegionOfInterestGeometry::CIRCLE_PATCH) {
                float radius = tgt::length(tgt::vec2(mousePos - lastMousePos_));
                if (radius > 0.f)
                    regionsOfInterest_.back()->points_ = generateCircle(lastMousePos_, radius);
            }
            else*/ if (regionsOfInterest_.back()->geometryMode_ == RegionOfInterestGeometry::RECTANGLE) {
                regionsOfInterest_.back()->points_.clear();
                if (tgt::length(tgt::vec2(mousePos - lastMousePos_)) > 0.f) {
                    regionsOfInterest_.back()->points_ = generateRectangle(tgt::vec2(lastMousePos_), tgt::vec2(mousePos));
                }
            }
            else {
                if (regionsOfInterest_.back()->points_.empty() || tgt::distance(point, regionsOfInterest_.back()->points_.back()) >= minRoiPointDistance_.get()) {
                    regionsOfInterest_.back()->points_.push_back(point);
                }
            }
        }
    }

    if (actions & FINISH_ROI) {
        if (!regionsOfInterest_.empty()) {
            // now render patch as patch
            /*if (geometryMode_.getValue() == RegionOfInterestGeometry::PATCH ||
                geometryMode_.getValue() == RegionOfInterestGeometry::CIRCLE_PATCH) {
                regionsOfInterest_.back()->geometryMode_ = RegionOfInterestGeometry::PATCH;
            }*/

            // compute bounding box of roi points
            regionsOfInterest_.back()->computeBoundingBox();

        }
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
        lastMousePos_ = tgt::ivec2(-1);
        geometryHasChanged_ = true;
        updateGeometry();
        numberOfROIs_.set(static_cast<int>(regionsOfInterest_.size()));
    }

    e->accept();
    invalidate();
}

void RegionOfInterest2D::removeLastROI(tgt::MouseEvent* e) {
    if (!regionsOfInterest_.empty()) {
        delete regionsOfInterest_.back();
        regionsOfInterest_.pop_back();
        if (selectedROI_ >= static_cast<int>(regionsOfInterest_.size()))
            selectedROI_ = -1;
        geometryHasChanged_ = true;
        updateGeometry();
        invalidate();
    }
    numberOfROIs_.set(static_cast<int>(regionsOfInterest_.size()));

    if (e)
        e->accept();
}

void RegionOfInterest2D::removeROIInCurrentSlice(tgt::MouseEvent* e) {
    vector<vector<RegionOfInterestGeometry*>::iterator> currentLayerGeoms;
    for (size_t i = 0; i < regionsOfInterest_.size(); ++i) {

        RegionOfInterestGeometry* roi = regionsOfInterest_[i];
        if (roi->layer_ == currentLayer_.get())
            currentLayerGeoms.push_back(regionsOfInterest_.begin() + i);
    }

    if (!currentLayerGeoms.empty()) {
        delete *(currentLayerGeoms.back());
        regionsOfInterest_.erase(currentLayerGeoms.back());

        if (selectedROI_ >= static_cast<int>(regionsOfInterest_.size()))
            selectedROI_ = -1;
        geometryHasChanged_ = true;
        updateGeometry();
        invalidate();
    }
    numberOfROIs_.set(static_cast<int>(regionsOfInterest_.size()));

    if (e)
        e->accept();
}

void RegionOfInterest2D::removeClosestROIInCurrentSlice(tgt::MouseEvent* e) {
    vector<RegionOfInterestGeometry*>::iterator closest = regionsOfInterest_.end();
    float distanceSq = FLT_MAX;

    tgt::vec2 volMousePos(e->x(), inport_.getSize().y - e->y());

    mat4 voxelToCurrentView;
    bool inv = roiTransformMatrix_.get().invert(voxelToCurrentView);
    tgtAssert(inv, "Matrix inversion failed");

    for (vector<RegionOfInterestGeometry*>::iterator it = regionsOfInterest_.begin(); it != regionsOfInterest_.end(); ++it) {

        RegionOfInterestGeometry* roi = *it;

        float localDistanceSq = FLT_MAX;

        if (roi->layer_ == currentLayer_.get()) {
            mat4 roiToCurrentView = voxelToCurrentView * roi->transformMatrix_;

            for(tgt::vec2 p : roi->points_) {
                tgt::vec2 pView = roiToCurrentView.transform(tgt::vec3(p, 0)).xy();
                localDistanceSq = std::min(localDistanceSq, tgt::distanceSq(pView, volMousePos));
            }

            if(localDistanceSq < distanceSq) {
                distanceSq = localDistanceSq;
                closest = it;
            }
        }
    }

    if (closest != regionsOfInterest_.end()) {
        delete *closest;
        regionsOfInterest_.erase(closest);

        if (selectedROI_ >= static_cast<int>(regionsOfInterest_.size()))
            selectedROI_ = -1;
        geometryHasChanged_ = true;
        updateGeometry();
        invalidate();
    }
    numberOfROIs_.set(static_cast<int>(regionsOfInterest_.size()));

    if (e)
        e->accept();
}

void RegionOfInterest2D::clearROIs(tgt::MouseEvent* e) {
    if (!regionsOfInterest_.empty()) {
        for (size_t i=0; i<regionsOfInterest_.size(); i++)
            delete regionsOfInterest_[i];
        regionsOfInterest_.clear();
        selectedROI_ = -1;
        geometryHasChanged_ = true;
        updateGeometry();
        invalidate();
    }
    numberOfROIs_.set(static_cast<int>(regionsOfInterest_.size()));

    if (e)
        e->accept();
}

void RegionOfInterest2D::shiftROI(tgt::MouseEvent* e){
    tgtAssert(e, "No mouse event");
    e->accept();

    if (!inport_.isReady())
        return;

    tgt::ivec2 mousePos(e->x(), inport_.getSize().y - e->y());

    if (e->action() == tgt::MouseEvent::PRESSED) {
        selectedROI_ = -1;
        lastMousePos_ = tgt::ivec2(-1);

        for (size_t i = 0; i < regionsOfInterest_.size() && (selectedROI_ == -1); ++i) {
            if (idManager_.isHit(mousePos, regionsOfInterest_[i])) {
                selectedROI_ = static_cast<int>(i);
                lastMousePos_ = mousePos;
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            }
        }

        geometryHasChanged_ = true;
        invalidate();
    }
    else if (e->action() == tgt::MouseEvent::MOTION) {
        if (selectedROI_ > -1 && lastMousePos_ != tgt::ivec2(-1)) {
            tgtAssert(selectedROI_ >= 0 && selectedROI_ < (int)regionsOfInterest_.size(), "Invalid selected roi");
            RegionOfInterestGeometry* roi = regionsOfInterest_.at(static_cast<size_t>(selectedROI_));
            vec2 offset = tgt::vec2(mousePos - lastMousePos_);
            if (offset != vec2(0.f)) {
                for (size_t i=0; i<roi->points_.size(); i++) {
                    roi->points_[i] += offset;
                }
                roi->boundingBox_ += tgt::vec4(offset.x, offset.y, offset.x, offset.y);

                lastMousePos_ = mousePos;
                geometryHasChanged_ = true;
                invalidate();
            }
        }
    }
    else if (e->action() == tgt::MouseEvent::RELEASED) {
        selectedROI_ = -1;
        lastMousePos_ = ivec2(-1);
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
        geometryHasChanged_ = true;
        invalidate();
    }
}

void RegionOfInterest2D::scaleROI(tgt::MouseEvent* e) {
    tgtAssert(e, "No mouse event");
    e->accept();

    if (!inport_.isReady())
        return;

    tgt::ivec2 mousePos(e->x(), inport_.getSize().y - e->y());

    if (e->action() == tgt::MouseEvent::PRESSED) {
        selectedROI_ = -1;
        lastMousePos_ = tgt::ivec2(-1);

        for (size_t i=0; i<regionsOfInterest_.size() && (selectedROI_ == -1); i++) {
            if (idManager_.isHit(mousePos, regionsOfInterest_[i])) {
                LINFO("Hit: " << i);
                selectedROI_ = static_cast<int>(i);
                lastMousePos_ = mousePos;
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            }
        }

        geometryHasChanged_ = true;
        invalidate();
    }
    else if (e->action() == tgt::MouseEvent::MOTION) {
        if (selectedROI_ > -1 && lastMousePos_ != tgt::ivec2(-1)) {
            tgtAssert(selectedROI_ >= 0 && selectedROI_ < (int)regionsOfInterest_.size(), "Invalid selected roi");
            RegionOfInterestGeometry* roi = regionsOfInterest_.at(static_cast<size_t>(selectedROI_));
            vec2 offset = tgt::vec2(mousePos - lastMousePos_);
            if (offset != tgt::vec2(0.f)) {
                vec2 bbCenter = vec2(roi->boundingBox_.z + roi->boundingBox_.x, roi->boundingBox_.w + roi->boundingBox_.y) / 2.f;
                vec2 scale = (1.f + (offset / vec2(inport_.getSize())) * 6.f);
                for (size_t i=0; i<roi->points_.size(); i++) {
                    roi->points_[i] = (roi->points_[i] - bbCenter) * scale + bbCenter;
                }
                roi->boundingBox_ = (roi->boundingBox_ - tgt::vec4(bbCenter.x, bbCenter.y, bbCenter.x, bbCenter.y)) * tgt::vec4(scale.x, scale.y, scale.x, scale.y)
                    + tgt::vec4(bbCenter.x, bbCenter.y, bbCenter.x, bbCenter.y);

                lastMousePos_ = mousePos;
                geometryHasChanged_ = true;
                invalidate();
            }
        }
    }
    else if (e->action() == tgt::MouseEvent::RELEASED) {
        selectedROI_ = -1;
        lastMousePos_ = ivec2(-1);
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
        geometryHasChanged_ = true;
        invalidate();
    }
}

void RegionOfInterest2D::serialize(Serializer& s) const {
    ImageProcessor::serialize(s);
    s.serialize("RegionsOfInterest", regionsOfInterest_, "Region");
}

void RegionOfInterest2D::deserialize(Deserializer& s) {
    ImageProcessor::deserialize(s);
    try {
        clearROIs();
        s.deserialize("RegionsOfInterest", regionsOfInterest_, "Region");
    }
    catch (SerializationNoSuchDataException&) {
        // just no ROI data present => ignore
        s.removeLastError();
    }
}

void RegionOfInterest2D::saveROIs(const std::string& filename) const {
    if (filename.empty()) {
        LERROR("ROI filename is empty");
        return;
    }

    LINFO("Saving regions of interest to file: " << filename);

    // open file stream
    std::ofstream stream(filename.c_str(), std::ios_base::out);
    if (stream.fail()) {
        LWARNING("Unable to open file " << filename << " for writing.");
        return;
    }

    // serialize to stream
    bool success = true;
    try {
        XmlSerializer s(filename);
        Serializer serializer(s);
        serializer.serialize("RegionsOfInterest", regionsOfInterest_, "Region");

        s.write(stream);
        if (stream.bad()) {
            LWARNING("Unable to write to file: " << filename);
            success = false;
        }
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNING("SerializationException: " << e.what());
        stream.close();
        success = false;
    }

    // log result
    if (!success)
        LWARNING("Saving regions of interest failed.");
}

void RegionOfInterest2D::loadROIs(const std::string& filename) {
    if (filename.empty()) {
        LERROR("ROI filename is empty");
        return;
    }

    LINFO("Loading regions of interest from file: " << filename);

    // open file stream
    std::ifstream stream(filename.c_str(), std::ios_base::in);
    if (stream.fail()) {
        LWARNING("Unable to open file " << filename << " for reading.");
        return;
    }

    // deserialize from stream
    bool success = true;
    try {
        XmlDeserializer d(filename);
        d.read(stream);
        Deserializer deserializer(d);
        deserializer.deserialize("RegionsOfInterest", regionsOfInterest_, "Region");
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNING("SerializationException: " << e.what());
        stream.close();
        success = false;
    }

    if (!success)
        LWARNING("Loading regions of interest failed.");
    numberOfROIs_.set(static_cast<int>(regionsOfInterest_.size()));

    geometryHasChanged_ = true;
    invalidate();
}

void RegionOfInterest2D::saveFileIssued() {
    if (!saveROIs_.get().empty()) {
        saveROIs(saveROIs_.get());
        saveROIs_.set("");
    }
}

void RegionOfInterest2D::loadFileIssued() {
    if (!loadROIs_.get().empty()) {
        loadROIs(loadROIs_.get());
        loadROIs_.set("");
    }
}

vector<vec2> RegionOfInterest2D::generateCircle(vec2 center, float radius) const {
    vector<vec2> result;
    int numPoints = tgt::ifloor(2.f*radius*tgt::PIf);
    float step = 2.f*tgt::PIf / numPoints;
    for (float t = 0.f; t <= 2.f*tgt::PIf; t += step) {
        vec2 nextPoint = radius * vec2(cos(t), sin(t)) + center;
        if (result.empty() || tgt::length(result.back() - nextPoint) >= minRoiPointDistance_.get())
            result.push_back(nextPoint);
    }
    return result;
}

vector<vec2> RegionOfInterest2D::generateRectangle(vec2 first, tgt::vec2 second) const {
    vector<vec2> result;
    result.push_back(first);
    result.push_back(vec2(first.x, second.y));
    result.push_back(second);
    result.push_back(vec2(second.x, first.y));
    return result;
}

void RegionOfInterest2D::updateGeometry() {
    // only update if the geometry has changed
    // (the second term fixes problems when the input volume changed and the network evaluator cleared the output volume)
    if(!geometryHasChanged_ && ( (!outportROIGeometry_.isConnected()) || (outportROIGeometry_.hasData()) ) )
        return;

    // fill regions where necessary
    /*for (size_t i = 0; i < regionsOfInterest_.size(); ++i) {
        RegionOfInterestGeometry* roi = regionsOfInterest_[i];
        if(fillRegions_.get() && roi->geometryMode_ == RegionOfInterestGeometry::PATCH) {
            maskPrivate_.resize(roi->viewportSize_);
            maskPrivate_.activateTarget();
            maskPrivate_.clearTarget();

            renderROIMask(roi);

            maskPrivate_.deactivateTarget();
            roi->fillArea(maskPrivate_.getColorTexture());
        }
    }*/

    // roi geometry (layer only)
    if (outportROIGeometryLayer_.isConnected()) {
        PointSegmentListGeometryVec3* segmentListGeom = new PointSegmentListGeometryVec3();
        vector<vector<vec3> > segmentList;
        for (size_t i = 0; i < regionsOfInterest_.size(); ++i) {

            RegionOfInterestGeometry* roi = regionsOfInterest_[i];
            if (roi->layer_ != currentLayer_.get())
                continue;

            /*if(fillRegions_.get() && roi->geometryMode_ == RegionOfInterestGeometry::PATCH) {
                maskPrivate_.resize(roi->viewportSize_);
                maskPrivate_.activateTarget();
                maskPrivate_.clearTarget();

                renderROIs(MASK_MODE);

                maskPrivate_.deactivateTarget();
                roi->fillArea(maskPrivate_.getColorTexture());
            }*/

            vector<vec3> pointList;
            vector<vec2> tmpList = fillRegions_.get() ? roi->fillPoints_ : roi->points_;
            for (size_t p = 0; p < tmpList.size(); ++p)
                pointList.push_back((roi->transformMatrix_ * vec4(tmpList[p], static_cast<float>(roi->layer_), 1.f)).xyz());
            segmentList.push_back(pointList);
        }
        segmentListGeom->setData(segmentList);
        outportROIGeometryLayer_.setData(segmentListGeom, true);
    }

    // roi geometry (all rois)
    if (outportROIGeometry_.isConnected()) {
        PointSegmentListGeometryVec3* segmentListGeom = new PointSegmentListGeometryVec3();
        vector<vector<vec3> > segmentList;
        for (size_t i = 0; i < regionsOfInterest_.size(); ++i) {
            RegionOfInterestGeometry* roi = regionsOfInterest_[i];
            vector<vec3> pointList;
            vector<vec2> tmpList = fillRegions_.get() ? roi->fillPoints_ : roi->points_;
            for (size_t p = 0; p < tmpList.size(); ++p)
                pointList.push_back((roi->transformMatrix_ * vec4(tmpList[p], static_cast<float>(roi->layer_), 1.f)).xyz());
            segmentList.push_back(pointList);
        }
        segmentListGeom->setData(segmentList);
        outportROIGeometry_.setData(segmentListGeom, true);
    }

    geometryHasChanged_ = false;
}

} // voreen namespace

