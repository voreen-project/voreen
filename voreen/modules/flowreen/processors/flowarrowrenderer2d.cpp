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

#include "flowarrowrenderer2d.h"

#include "tgt/tgt_math.h"
#include <sstream>

#include "tgt/gpucapabilities.h"
#include "tgt/glmath.h"
#include "tgt/font.h"
#include "tgt/tgt_gl.h"
#include "tgt/textureunit.h"

#include "voreen/core/datastructures/octree/volumeoctreebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/slice/slicehelper.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

using tgt::TextureUnit;

namespace voreen {

const std::string FlowArrowRenderer2D::loggerCat_("voreen.bovenkamp.FlowArrowRenderer2D");
const std::string FlowArrowRenderer2D::fontName_("Vera.ttf");

FlowArrowRenderer2D::FlowArrowRenderer2D()
    : VolumeRenderer()
    , velocityVolumeInport_(Port::INPORT, "velocityVolumeInport", "Velocity Input")
    , magnitudeVolumeInport_(Port::INPORT, "magnitudeVolumeInport", "Magnitude Input")
    , renderOutport_(Port::OUTPORT, "renderOutport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)

    , mappingOption_("mappingOption","Derive arrow color...")
    , arrowSize_("arrowSize","Arrow Size Factor",1.f,0.01f,10.f)
    //magnitude setting
    , magnitudeTF_("magnitudeTF", "Magnitude")
    //arrow setting
    , arrowTF_("arrowTF", "Arrow")

    , sliceAlignment_("sliceAlignmentProp", "Slice Alignment")
    , sliceIndex_("sliceIndex", "Slice Number ", 0, 0, 10000)
    , mouseXCoord_("sliceX", "Slice Mouse Pos X", 0, 0, 10000)
    , mouseYCoord_("sliceY", "Slice Mouse Pos Y", 0, 0, 10000)
    , mouseZCoord_("sliceZ", "Slice Mouse Pos Z", 0, 0, 10000)
    , renderSliceBoundaries_("renderSliceBoundaries", "Render Slice Boundaries", true)
    , boundaryColor_("boundaryColor", "Boundary Color", tgt::Color(1.f, 1.f, 1.f, 1.f))

    , showCursorInfos_("showCursorInformation", "Show Cursor Info")
    , showSliceNumber_("showSliceNumber", "Show Slice Number", true)
    , showScaleLegend_("renderLegend", "Show Scale Legend", false)
    , fontSize_("fontSize", "Font Size", 14, 8, 48)
    , voxelOffset_("voxelOffset", "Voxel Offset", tgt::vec2(0.f), tgt::vec2(-10000.f), tgt::vec2(10000.f))
    , zoomFactor_("zoomFactor", "Zoom Factor", 1.f, 0.01f, 1.f)
    , pickingMatrix_("pickingMatrix", "Picking Matrix", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f), Processor::VALID)
    , mwheelCycleHandler_("mouseWheelHandler", "Slice Cycling", &sliceIndex_)
    , mwheelZoomHandler_("zoomHandler", "Slice Zoom", &zoomFactor_, tgt::MouseEvent::CTRL)
    , magnitudeShader_(0), velocityShader_(0)
    , voxelPosPermutation_(0, 1, 2)
    , sliceLowerLeft_(1.f)
    , sliceSize_(0.f)
    , mousePosition_(-1, -1)
    , mouseIsPressed_(false)
    , lastPickingPosition_(-1, -1, -1)
{
    magnitudeVolumeInport_.addCondition(new PortConditionVolumeTypeGL());
    magnitudeVolumeInport_.showTextureAccessProperties(true);
    addPort(magnitudeVolumeInport_);
    addPort(velocityVolumeInport_);
    addPort(renderOutport_);

    // interaction
    mouseEventShift_ = new EventProperty<FlowArrowRenderer2D>("mouseEvent.Shift", "Slice Shift",
        this, &FlowArrowRenderer2D::shiftEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::MOTION, tgt::Event::CTRL);

    mouseEventMove_ = new EventProperty<FlowArrowRenderer2D>("mouseEvent.cursorPositionMove", "Cursor Position Move",
        this, &FlowArrowRenderer2D::mouseLocalization,
        tgt::MouseEvent::MOUSE_BUTTON_NONE,
        tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE,
        true);

    mouseEventPress_ = new EventProperty<FlowArrowRenderer2D>("mouseEvent.cursorPositionPress", "Cursor Position Press",
        this, &FlowArrowRenderer2D::mouseLocalization,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_MIDDLE),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED | tgt::MouseEvent::WHEEL | tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE,
        true);

    addEventProperty(mouseEventPress_);
    addEventProperty(mouseEventMove_);
    addEventProperty(mouseEventShift_);

    addInteractionHandler(mwheelCycleHandler_);
    addInteractionHandler(mwheelZoomHandler_);

    //color/arrow setting
    addProperty(mappingOption_);
        mappingOption_.addOption("fromtf", "from Arrow Length", false);
        mappingOption_.addOption("fromdir", "from Arrow Direction", true);

    addProperty(arrowSize_);
    addProperty(magnitudeTF_);

    addProperty(arrowTF_);


    // slice arrangement
    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    sliceAlignment_.onChange(
        MemberFunctionCallback<FlowArrowRenderer2D>(this, &FlowArrowRenderer2D::onSliceAlignmentChange) );
    addProperty(sliceAlignment_);

    addProperty(sliceIndex_);

    addProperty(renderSliceBoundaries_);
    addProperty(boundaryColor_);
    addProperty(mouseXCoord_);
    addProperty(mouseYCoord_);
    addProperty(mouseZCoord_);

    // group slice arrangement properties
    sliceAlignment_.setGroupID("sliceArrangement");
    sliceIndex_.setGroupID("sliceArrangement");
    renderSliceBoundaries_.setGroupID("sliceArrangement");
    boundaryColor_.setGroupID("sliceArrangement");
    mouseXCoord_.setGroupID("sliceArrangement");
    mouseYCoord_.setGroupID("sliceArrangement");
    mouseZCoord_.setGroupID("sliceArrangement");

    mouseXCoord_.setVisibleFlag(false);
    mouseYCoord_.setVisibleFlag(false);
    mouseZCoord_.setVisibleFlag(false);

    setPropertyGroupGuiName("sliceArrangement", "Slice Selection");

    // information overlay
    showCursorInfos_.addOption("never", "Never");
    showCursorInfos_.addOption("onClick", "On Mouse Click");
    showCursorInfos_.addOption("onDrag", "On Mouse Drag");
    showCursorInfos_.addOption("onMove", "On Mouse Move");
    showCursorInfos_.select("onMove");
    addProperty(showCursorInfos_);
    addProperty(showSliceNumber_);
    addProperty(showScaleLegend_);
    addProperty(fontSize_);

    // group information overlay properties
    showCursorInfos_.setGroupID("informationOverlay");
    showSliceNumber_.setGroupID("informationOverlay");
    fontSize_.setGroupID("informationOverlay");
    showScaleLegend_.setGroupID("informationOverlay");
    setPropertyGroupGuiName("informationOverlay", "Information Overlay");

    // zooming
    addProperty(voxelOffset_);
    zoomFactor_.setStepping(0.01f);
    addProperty(zoomFactor_);
    pickingMatrix_.setReadOnlyFlag(true);
    addProperty(pickingMatrix_);

    // group zooming props
    voxelOffset_.setGroupID("zooming");
    zoomFactor_.setGroupID("zooming");
    pickingMatrix_.setGroupID("zooming");
    setPropertyGroupGuiName("zooming", "Zooming");

    // call this method to set the correct permutation for the
    // screen-position-to-voxel-position mapping.
    onSliceAlignmentChange();
}

FlowArrowRenderer2D::~FlowArrowRenderer2D() {
    delete mouseEventPress_;
    delete mouseEventMove_;
    delete mouseEventShift_;
}

Processor* FlowArrowRenderer2D::create() const {
    return new FlowArrowRenderer2D();
}

void FlowArrowRenderer2D::initialize() {
    VolumeRenderer::initialize();

    magnitudeShader_ = ShdrMgr.loadSeparate("slice.vert","","slice.frag", generateHeader(), false);
    LGL_ERROR;
    velocityShader_ = ShdrMgr.loadSeparate("velocity.vert","velocity.geom","velocity.frag", generateHeader(), false);
    LGL_ERROR;

    updatePropertyConfiguration();
}

void FlowArrowRenderer2D::deinitialize() {
    ShdrMgr.dispose(magnitudeShader_);
    magnitudeShader_ = 0;
    ShdrMgr.dispose(velocityShader_);
    velocityShader_ = 0;

    VolumeRenderer::deinitialize();
}

void FlowArrowRenderer2D::beforeProcess() {
    VolumeRenderer::beforeProcess();

    if (magnitudeVolumeInport_.hasChanged() || velocityVolumeInport_.hasChanged()) {
        mousePosition_ = tgt::ivec2(-1);
        lastPickingPosition_ = tgt::ivec3(-1);
        updatePropertyConfiguration();

        magnitudeTF_.setVolume(magnitudeVolumeInport_.getData(), 0);
       //arrowTF_.setVolumeHandle(velocityVolumeInport_.getData(), 0);
    }

    if (invalidationLevel_ >= Processor::INVALID_PROGRAM || magnitudeVolumeInport_.hasChanged() || velocityVolumeInport_.hasChanged())
        rebuildShader();
}

void FlowArrowRenderer2D::updatePropertyConfiguration() {
    if (!magnitudeVolumeInport_.hasData() || !velocityVolumeInport_.hasData())
        return;

    tgt::svec3 volumeDim = velocityVolumeInport_.getData()->getDimensions();

    tgtAssert(sliceAlignment_.getValue() >= 0 && sliceAlignment_.getValue() <= 2, "Invalid alignment value");
    size_t numSlices = volumeDim[sliceAlignment_.getValue()];
    if (numSlices == 0)
        return;

    sliceIndex_.setMaxValue((int)numSlices-1);
    if (sliceIndex_.get() >= static_cast<int>(numSlices))
        sliceIndex_.set(static_cast<int>(numSlices / 2));

    mouseXCoord_.setMaxValue((int)volumeDim.x);
    mouseYCoord_.setMaxValue((int)volumeDim.y);
    mouseZCoord_.setMaxValue((int)volumeDim.z);

    tgt::vec2 halfDim = tgt::vec2((float)volumeDim[voxelPosPermutation_.x], (float)volumeDim[voxelPosPermutation_.y] - 1.f) / 2.f;
    voxelOffset_.setMinValue(-halfDim);
    voxelOffset_.setMaxValue(halfDim);
    voxelOffset_.set(tgt::clamp(voxelOffset_.get(), voxelOffset_.getMinValue(), voxelOffset_.getMaxValue()));
}

void FlowArrowRenderer2D::process() {
    renderOutport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    const VolumeBase* magVolume = magnitudeVolumeInport_.getData();

    const SliceAlignment alignment = sliceAlignment_.getValue();
    // get voxel volume dimensions
    tgt::ivec3 volDim = magVolume->getDimensions();
    int numSlices = volDim[alignment];

    // get the textures dimensions
    tgt::vec3 urf = magnitudeVolumeInport_.getData()->getURB();
    tgt::vec3 llb = magnitudeVolumeInport_.getData()->getLLF();

    // Re-calculate texture dimensions, urf, llb and center of the texture
    // for it might be a NPOT texture and therefore might have been inflated.
    // In that case, the inflating need to be undone and the volume texture
    // needs to be "cropped" to its original measures.
    const tgt::vec3 texDim = urf - llb;
    urf = texDim / 2.f;
    llb = texDim / -2.f;
    tgt::vec3 texCenter = (llb + (texDim * 0.5f));

    // Use OpenGL's ability of multiplying texture coordinate vectors with a matrix
    // on the texture matrix stack to permute the components of the texture
    // coordinates to obtain a correct setting for the current slice alignment.
    textureMatrix_ = tgt::mat4::zero;
    tgt::vec2 texDim2D(0.0f);   // laziness... (no matrix multiplication to determine 2D size of texture)

    // set slice's width and height according to currently slice alignment
    // and set the pointer to the slice rendering method to the matching version
    //
    switch (alignment) {
        case XZ_PLANE:
            texDim2D.x = texDim.x;
            texDim2D.y = texDim.z;
            textureMatrix_.t00 = 1.f;   // setup a permutation matrix, swapping z- and y-
            textureMatrix_.t12 = 1.f;   // components on vectors being multiplied with it
            textureMatrix_.t21 = 1.f;
            textureMatrix_.t33 = 1.f;
            break;

        case YZ_PLANE:
            texDim2D.x = texDim.y;
            texDim2D.y = texDim.z;
            textureMatrix_.t02 = 1.f;    // setup a permutation matrix, swapping x-, y- and z-
            textureMatrix_.t10 = -1.f;   // components on vectors being multiplied with it
            textureMatrix_.t21 = 1.f;
            textureMatrix_.t13 = 1.f;
            textureMatrix_.t33 = 1.f;
            break;

        case XY_PLANE:
            texDim2D.x = texDim.x;
            texDim2D.y = texDim.y;
            textureMatrix_ = tgt::mat4::identity;
            textureMatrix_.t11 = -1.f;  // invert y-axis, since
            textureMatrix_.t13 = 1.f;   // we want to look along positive z-axis
            textureMatrix_.t33 = 1.f;
        default:
            break;
    }   // switch

    float canvasWidth = static_cast<float>(renderOutport_.getSize().x);
    float canvasHeight = static_cast<float>(renderOutport_.getSize().y);
    float scaleWidth = canvasWidth / (texDim2D.x);
    float scaleHeight = canvasHeight / (texDim2D.y);

    // find minimal scaling factor (either scale along canvas' width or
    // canvas' height)
    if (scaleWidth <= scaleHeight) {

        sliceSize_ = texDim2D * scaleWidth * (1.f/zoomFactor_.get());

        sliceLowerLeft_.x = 0;
        sliceLowerLeft_.y = (canvasHeight - sliceSize_.y) / 2.f;

        // adapt for zooming
        tgt::vec2 zoomOffset = -((sliceSize_ - (texDim2D * scaleWidth)) / 2.f);
        tgt::vec2 volDimFloat = tgt::vec2((float)volDim[voxelPosPermutation_.x], (float)volDim[voxelPosPermutation_.y]) - 1.f;
        tgt::vec2 focusOffset = voxelOffset_.get() / volDimFloat;
        focusOffset *= sliceSize_;
        sliceLowerLeft_.x += zoomOffset.x + focusOffset.x;
        sliceLowerLeft_.y += focusOffset.y;
    }
    else {
        sliceSize_ = texDim2D * scaleHeight * (1.f/zoomFactor_.get());

        sliceLowerLeft_.x = (canvasWidth - sliceSize_.x) / 2.f;
        sliceLowerLeft_.y = 0.f;

        // adapt for zooming
        tgt::vec2 zoomOffset = -((sliceSize_ - (texDim2D * scaleHeight)) / 2.f);
        tgt::vec2 volDimFloat = tgt::vec2((float)volDim[voxelPosPermutation_.x], (float)volDim[voxelPosPermutation_.y]) - 1.f;
        tgt::vec2 focusOffset = voxelOffset_.get() / volDimFloat;
        focusOffset *= sliceSize_;
        sliceLowerLeft_.x += focusOffset.x;
        sliceLowerLeft_.y += zoomOffset.y + focusOffset.y;
    }

    // setup matrices
    LGL_ERROR;

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    glOrtho(0.0f, canvasWidth, 0.f, canvasHeight, -1.0f, 1.0f);

    LGL_ERROR;

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    LGL_ERROR;

    // setup shader
    magnitudeShader_->activate();
    LGL_ERROR;
    // bind transfer functions
    TextureUnit texUnit, transferUnitMag, transferUnitVel;
    transferUnitMag.activate();
    magnitudeTF_.get()->getTexture()->bind();
    magnitudeTF_.get()->setUniform(magnitudeShader_, "transFuncParams_", "transFuncTex_", transferUnitMag.getUnitNumber());
    LGL_ERROR;
    magnitudeShader_->deactivate();

    // render slice(s) as textured quads
    float depth = 0.0f;
    const size_t sliceIndex = static_cast<size_t>(sliceIndex_.get());

    int sliceNumber = (static_cast<const int>(sliceIndex));
    // calculate depth in llb/urf space
    depth = ((static_cast<float>(sliceNumber) / static_cast<float>(std::max(numSlices - 1, 1))) - 0.5f)
            * texDim[alignment];

    // map depth to [0, 1]
    depth -= texCenter[alignment];  // center around origin
    depth /= texDim[alignment];     // map to [-0.5, -0.5]
    depth += 0.5f;                  // map to [0, 1]

    MatStack.loadIdentity();
    MatStack.translate(sliceLowerLeft_.x, sliceLowerLeft_.y, 0.0f);
    MatStack.scale(sliceSize_.x, sliceSize_.y, 1.0f);

    if(!magnitudeShader_->isActivated())
        magnitudeShader_->activate();

    // render slice
    // extract 2D slice
    SliceTexture* slice = SliceHelper::getVolumeSlice(magVolume, alignment, tgt::iround(depth*(magVolume->getDimensions()[alignment]-1)));
    tgtAssert(slice,"no slice");
    // bind slice texture
    GLint texFilterMode = magnitudeVolumeInport_.getTextureFilterModeProperty().getValue();
    GLint texClampMode = magnitudeVolumeInport_.getTextureClampModeProperty().getValue();
    tgt::vec4 borderColor = tgt::vec4(magnitudeVolumeInport_.getTextureBorderIntensityProperty().get());
    GLSL::bindSliceTexture(slice, &texUnit, texFilterMode, texClampMode, borderColor);

    // pass slice uniforms to shader
    magnitudeShader_->setIgnoreUniformLocationError(true);
    GLSL::setUniform(magnitudeShader_, "sliceTex_", "sliceTexParams_", slice, &texUnit);
    magnitudeShader_->setIgnoreUniformLocationError(false);

    // render slice
    tgt::vec3 texLowerLeft = (textureMatrix_*tgt::vec4(0.f, 0.f, depth, 1.f)).xyz();
    tgt::vec3 texUpperRight = (textureMatrix_*tgt::vec4(1.f, 1.f, depth, 1.f)).xyz();
    magnitudeShader_->setUniform("textureMatrix_", tgt::mat4::identity);
    glBegin(GL_QUADS);
        glTexCoord2f(texLowerLeft[ (alignment+1) % 3], texLowerLeft[ (alignment+2) % 3]);  glVertex2f(0.0f, 0.0f);
        glTexCoord2f(texUpperRight[(alignment+1) % 3], texLowerLeft[ (alignment+2) % 3]);  glVertex2f(1.0f, 0.0f);
        glTexCoord2f(texUpperRight[(alignment+1) % 3], texUpperRight[(alignment+2) % 3]);  glVertex2f(1.0f, 1.0f);
        glTexCoord2f(texLowerLeft[ (alignment+1) % 3], texUpperRight[(alignment+2) % 3]);  glVertex2f(0.0f, 1.0f);
    glEnd();
    LGL_ERROR;

    delete slice;
    LGL_ERROR;

    tgtAssert(magnitudeShader_, "no slice shader");
    magnitudeShader_->deactivate();

    // render arrows
    //
    //bind tftexture
    const VolumeBase* velVolume = velocityVolumeInport_.getData();
    tgt::Texture* tex = arrowTF_.get()->getTexture();
    transferUnitVel.activate();
    tex->bind();
    LGL_ERROR;
    const VolumeRAM_3xFloat* velValues = dynamic_cast<const VolumeRAM_3xFloat*>(velVolume->getRepresentation<VolumeRAM>());
    velocityShader_->activate();
    velocityShader_->setUniform("maxVelocity", velValues->maxNormalizedMagnitude());
    velocityShader_->setUniform("ColorTexture", transferUnitVel.getUnitNumber());
    velocityShader_->setUniform("colorFromDir_", mappingOption_.getValue());

    float xStep = 1.f/velVolume->getDimensions().x;
    float yStep = 1.f/velVolume->getDimensions().y;
    float zStep = 1.f/velVolume->getDimensions().z;
    glDisable(GL_DEPTH_TEST);
    switch (alignment) {
        case XZ_PLANE:
            velocityShader_->setUniform("sliceAlignement_", tgt::vec3(0.f,1.f,0.f));
            velocityShader_->setUniform("hight_", std::min(xStep,zStep)*arrowSize_.get());
            glBegin(GL_POINTS);
                VRN_FOR_EACH_VOXEL(pos,tgt::svec3(0,sliceIndex,0),tgt::svec3(velVolume->getDimensions().x,sliceIndex+1,velVolume->getDimensions().z)) {
                    glNormal3fv(velValues->voxel(pos).elem);
                    glVertex3f((0.5f+pos.x)*xStep,((0.5f+pos.z)*zStep),-0.5);
                }
            glEnd();
            break;

        case YZ_PLANE:
            velocityShader_->setUniform("sliceAlignement_", tgt::vec3(0.f,1.f,0.f));
            velocityShader_->setUniform("hight_", std::min(yStep,zStep)*arrowSize_.get());
            glBegin(GL_POINTS);
                VRN_FOR_EACH_VOXEL(pos,tgt::svec3(sliceIndex,0,0),tgt::svec3(sliceIndex+1,velVolume->getDimensions().y,velVolume->getDimensions().z)) {
                    glNormal3fv(velValues->voxel(pos).elem);
                    glVertex3f(1-(0.5f+pos.y)*yStep,((0.5f+pos.z)*zStep),0.0);
                }
            glEnd();
            break;

        case XY_PLANE:
            velocityShader_->setUniform("sliceAlignement_", tgt::vec3(0.f,0.f,1.f));
            velocityShader_->setUniform("hight_", std::min(xStep,yStep)*arrowSize_.get());
            glBegin(GL_POINTS);
                VRN_FOR_EACH_VOXEL(pos,tgt::svec3(0,0,sliceIndex),tgt::svec3(velVolume->getDimensions().x,velVolume->getDimensions().y,sliceIndex+1)) {
                    glNormal3fv(velValues->voxel(pos).elem);
                    glVertex3f((0.5f+pos.x)*xStep,1-((0.5f+pos.y)*yStep),0.0); //flip y???
                }
            glEnd();
        default:
            break;
    }   // switch

    glEnable(GL_DEPTH_TEST);

    velocityShader_->deactivate();

    // render a border around each slice's boundaries if desired
    //
    if (renderSliceBoundaries_.get()) {
        MatStack.loadIdentity();
        MatStack.translate(sliceLowerLeft_.x, sliceLowerLeft_.y, 0.0f);
        MatStack.scale(sliceSize_.x, sliceSize_.y, 1.0f);
        glDepthFunc(GL_ALWAYS);
        renderSliceBoundaries();
        glDepthFunc(GL_LESS);
    }

    // If freetype is available render the slice's number and cursor information, otherwise this will do nothing
    renderInfoTexts();
    LGL_ERROR;

    // render legend
    if (showScaleLegend_.get()){
        renderLegend();
        LGL_ERROR;
    }

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;

    glActiveTexture(GL_TEXTURE0);
    renderOutport_.deactivateTarget();
    LGL_ERROR;

    // propagate picking matrix, if in single slice mode and a consumer is connected
    if (singleSliceMode() && !pickingMatrix_.getLinks().empty())
        pickingMatrix_.set(generatePickingMatrix());
}

std::string FlowArrowRenderer2D::generateHeader(const tgt::GpuCapabilities::GlVersion*) {
    std::string header = VolumeRenderer::generateHeader();

    header += magnitudeTF_.get()->getShaderDefines();

    return header;
}

bool FlowArrowRenderer2D::rebuildShader() {
    // do nothing if there is no shader at the moment
    if (!magnitudeShader_ || !velocityShader_)
        return false;

    magnitudeShader_->setHeaders(generateHeader());
    velocityShader_->setHeaders(generateHeader());
    return (velocityShader_->rebuild() & magnitudeShader_->rebuild());
}

void FlowArrowRenderer2D::renderLegend() {
    //TODO: adapt to new font rendering.
    /*
    tgtAssert(velocityVolumeInport_.getData(), "No volume");
    const VolumeBase* volume = velocityVolumeInport_.getData();

    // ESSENTIAL: if you don't use this, your text will become texturized!
    glActiveTexture(GL_TEXTURE0);
    glColor4f(1.f, 1.f, 1.f, 1.f);
    glDisable(GL_DEPTH_TEST);
    MatStack.loadIdentity();
    LGL_ERROR;


    tgt::Font legendFont(VoreenApplication::app()->getFontPath(fontName_));
    // note: the font size may not be smaller than 8
    legendFont.setSize(fontSize_.get());

    float legendOffsetX = 20.f;
    float legendOffsetY = 20.f;
    float lineLength = 80.f;
    float lineWidth = 5.f;

    // determine scale String
    tgt::mat4 matrix = generatePickingMatrix();
    tgt::vec4 first = matrix*tgt::vec4(0.f,0.f,0.f,1.f);
    tgt::vec4 second = matrix*tgt::vec4(1.f,0.f,0.f,1.f);
    float scale = tgt::max(tgt::abs((first.xyz() - second.xyz())*volume->getSpacing()));
    std::string scaleStr = formatSpatialLength(scale*lineLength);

    // determine bounds and render the string
    tgt::Bounds bounds = legendFont.getBounds(tgt::vec3(0.f, 0.f, 0.f), scaleStr);
    float textLength = bounds.getURB().x - bounds.getLLF().x;
    float textHight = bounds.getURB().y - bounds.getLLF().y;
    legendFont.render(tgt::vec3(std::max(legendOffsetX,legendOffsetX+(lineLength-textLength)/2.f), legendOffsetY, 0), scaleStr);

    float lineOffsetY = legendOffsetY + textHight + 5.f;
    float lineOffsetX = std::max(legendOffsetX,legendOffsetX+(textLength-lineLength)/2.f);

    //render line
    glBegin(GL_LINES);
        glVertex2f(lineOffsetX,lineOffsetY);
        glVertex2f(lineOffsetX+lineLength,lineOffsetY);
        glVertex2f(lineOffsetX,lineOffsetY-lineWidth);
        glVertex2f(lineOffsetX,lineOffsetY+lineWidth);
        glVertex2f(lineOffsetX+lineLength,lineOffsetY-lineWidth);
        glVertex2f(lineOffsetX+lineLength,lineOffsetY+lineWidth);
    glEnd();

    //reset all
    glEnable(GL_DEPTH_TEST);
    glColor4f(0.f, 0.f, 0.f, 0.f);
    LGL_ERROR;
    */
}

void FlowArrowRenderer2D::renderSliceBoundaries() const {
    glColor4f(boundaryColor_.get().r, boundaryColor_.get().g, boundaryColor_.get().b, boundaryColor_.get().a);
    glDisable(GL_DEPTH_TEST);
    glBegin(GL_LINE_LOOP);
        glVertex2f(0.0f, 0.0f);
        glVertex2f(1.0f, 0.0f);
        glVertex2f(1.0f, 1.0f);
        glVertex2f(0.0f, 1.0f);
    glEnd();
    glEnable(GL_DEPTH_TEST);
    glColor4f(0.f, 0.f, 0.f, 0.f);

    LGL_ERROR;
}

void FlowArrowRenderer2D::renderInfoTexts() const {

    //TODO: adapt to new font rendering.
    /*
    if (showCursorInfos_.isSelected("never") && !showSliceNumber_.get())
        return;

    tgtAssert(velocityVolumeInport_.getData(), "No volume");
    const VolumeBase* velVolume = velocityVolumeInport_.getData();
    const VolumeBase* magVolume = magnitudeVolumeInport_.getData();
    tgt::ivec3 volDim = velVolume->getDimensions();
    int numSlices = volDim[sliceAlignment_.getValue()];

    glDisable(GL_DEPTH_TEST);

    // ESSENTIAL: if you don't use this, your text will become texturized!
    glActiveTexture(GL_TEXTURE0);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    LGL_ERROR;

    // render voxel position information
    if (!showCursorInfos_.isSelected("never")) {
        tgt::Font fontCursorInfos(VoreenApplication::app()->getFontPath(fontName_));
        // note: the font size may not be smaller than 8
        fontCursorInfos.setSize(fontSize_.get());

        // voxel position
        tgt::ivec3 voxelPos = tgt::iround(screenToVoxelPos(mousePosition_));
        if (mousePosition_.x != -1 && voxelPos.x != -1) { // save cursor position, if it is valid
            lastPickingPosition_ = voxelPos;
        }

        // render cursor information, if a valid picking position is available
        if (lastPickingPosition_.x != -1) {
            lastPickingPosition_ = tgt::clamp(lastPickingPosition_, tgt::ivec3(0), volDim-1);

            // determine renderString
            std::string renderStr;
            std::ostringstream oss;
            oss << "[" << lastPickingPosition_.x << " " << lastPickingPosition_.y << " " << lastPickingPosition_.z << "]: ";
            RealWorldMapping rwm = velocityVolumeInport_.getData()->getRealWorldMapping();
            if (velVolume->hasRepresentation<VolumeRAM>()) { //< TODO: use other representations
                const VolumeRAM* volume = velVolume->getRepresentation<VolumeRAM>();
                oss << volume->getVoxelValueAsString(lastPickingPosition_, &rwm);
            }
            renderStr = oss.str();

            // determine bounds and render the string
            tgt::Bounds bounds = fontCursorInfos.getBounds(tgt::vec3(6.f, 6.f, 0.f), renderStr);
            float textHeight = bounds.getURB().y - bounds.getLLF().y;
            MatStack.loadIdentity();
            MatStack.translate(0.f, renderOutport_.getSize().y - textHeight - 25.f, 0.0f);
            fontCursorInfos.render(tgt::vec3(6.f, 6.f, 0), renderStr);

            // render spatial position
            std::string spatialPosStr = formatSpatialLength(tgt::vec3(lastPickingPosition_)*velVolume->getSpacing());
            MatStack.loadIdentity();
            MatStack.translate(0.f, renderOutport_.getSize().y - 2.f*textHeight - 28.f, 0.0f);
            fontCursorInfos.render(tgt::vec3(6.f, 6.f, 0), spatialPosStr);

            MatStack.loadIdentity();
        }

        LGL_ERROR;
    }

    if (showSliceNumber_.get()) {
        // note: the font size may not be smaller than 8
        tgt::Font fontSliceNumber(VoreenApplication::app()->getFontPath(fontName_));
        fontSliceNumber.setSize(fontSize_.get());

        // Therefore calculate the rendered text's bounding box by using a dummy
        // string.
        //
        std::string prefix10;
        std::string prefix100;
        std::string prefix1000;
        std::string dummy;
        if (numSlices < 10) {
            dummy = "8/8";
            prefix10 = "";
        }
        else if (numSlices < 100) {
            dummy = "88/88";
            prefix10 = "0";
        }
        else if (numSlices < 1000) {
            dummy = "888/888";
            prefix10 = "00";
            prefix100 = "0";
        }
        else {
            dummy = "8888/8888";
            prefix10 = "000";
            prefix100 = "00";
            prefix1000 = "0";
        }

        tgt::Bounds bounds = fontSliceNumber.getBounds(tgt::vec3(6.f, 6.f, 0.f), dummy);
        float textWidth = ceilf(bounds.getURB().x - bounds.getLLF().x);

        // Do not render the slice numbers if the slice width becomes too small to
        // prevent FTGL from creating OpenGL state errors.
        //
        if ((floorf(sliceSize_.x) > textWidth)) {
            // render the slice number information
            int sliceNumber = sliceIndex_.get();
            std::string prefix;

            MatStack.loadIdentity();
            MatStack.translate(sliceLowerLeft_.x + sliceSize_.x - textWidth - 12.f, sliceLowerLeft_.y , 0.0f);

            std::ostringstream oss;
            if ((sliceNumber) < 10)
                oss << prefix10;
            else if ((sliceNumber) < 100)
                oss << prefix100;
            else if ((sliceNumber) < 1000)
                oss << prefix1000;

            oss << (sliceNumber) << "/" << numSlices - 1;

            fontSliceNumber.render(tgt::vec3(6,6,0), oss.str());
            LGL_ERROR;
        }
    }

    MatStack.loadIdentity();
    glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
    glEnable(GL_DEPTH_TEST);
    LGL_ERROR;
    */
}

tgt::vec3 FlowArrowRenderer2D::screenToVoxelPos(tgt::ivec2 screenPos) const {

    if (!velocityVolumeInport_.getData() || !renderOutport_.getRenderTarget())
        return tgt::vec3(-1.f);

    tgt::vec3 volumeDim(velocityVolumeInport_.getData()->getDimensions());
    tgt::ivec2 screenDim = renderOutport_.getSize();

    tgt::ivec2 p(0, 0);
    p.x = screenPos.x - static_cast<int>(tgt::round(sliceLowerLeft_.x));
    p.y = (screenDim.y - screenPos.y) - static_cast<int>(tgt::round(sliceLowerLeft_.y));

    // if coordinates are negative, no slice could be hit
    if (tgt::hor(tgt::lessThan(p, tgt::ivec2(0))))
        return tgt::vec3(-1.f);

    const tgt::ivec2 sliceSizeInt = static_cast<tgt::ivec2>(sliceSize_);

    // if coordinates are greater than the number of slices per direction
    // times their extension in that direction, no slice could be hit either
    if ((p.x >= (sliceSizeInt.x)) || (p.y >= (sliceSizeInt.y)))
        return tgt::vec3(-1.f);

    // determine the picked slice
    const int sliceColID = p.x / sliceSizeInt.x;
    const int sliceRowID = -(p.y / sliceSizeInt.y);
    const int slice = sliceColID + sliceRowID + sliceIndex_.get();

    // calculate the normalized position within the picked slice
    tgt::vec2 posWithinSlice(
        static_cast<float>(p.x % sliceSizeInt.x),
        static_cast<float>(p.y % sliceSizeInt.y));
    posWithinSlice /= sliceSize_;

    // calculate the normalized depth of the picked slice (texture z coordinate)
    float depth = slice / std::max(volumeDim[voxelPosPermutation_.z] - 1.f, 1.f);

    // now we have the assigned texture coordinates of the picked fragment
    tgt::vec4 texCoords(posWithinSlice, depth, 1.f);
    texCoords = tgt::clamp(texCoords, tgt::vec4(0.f), tgt::vec4(1.f));

    // apply current texture matrix to assigned tex coords
    tgt::vec3 texCoordsTransformed = (textureMatrix_ * texCoords).xyz();

    // transform final tex coords into volume coordinates
    tgt::vec3 voxPos = texCoordsTransformed * (volumeDim-1.f);
    voxPos = tgt::clamp(voxPos, tgt::vec3(0.f), tgt::vec3(volumeDim-1.f));

    return voxPos;
}

tgt::mat4 FlowArrowRenderer2D::generatePickingMatrix() const {

    if (!velocityVolumeInport_.hasData())
        return tgt::mat4::createIdentity();

    tgt::vec3 volumeDim(velocityVolumeInport_.getData()->getDimensions());

    // 1. translate slice to origin
    tgt::mat4 originTranslation = tgt::mat4::createTranslation(tgt::vec3(-sliceLowerLeft_.x, -sliceLowerLeft_.y, 0.f));

    // 2. normalize screen coords with regard to the slice
    tgt::mat4 sliceScale = tgt::mat4::createScale(tgt::vec3(1.f / sliceSize_.x, 1.f / sliceSize_.y, 1.f / (volumeDim[voxelPosPermutation_.z] - 1.f)));

    // 3. apply current texture matrix
    //tgt::mat4 textureMatrix = textureMatrix_; // not used

    // 4. scale normalized coordinates to volume dimensions
    tgt::mat4 volumeScale = tgt::mat4::createScale(volumeDim - 1.f);

    // compose transformation matrix
    tgt::mat4 result = volumeScale * textureMatrix_ * sliceScale * originTranslation;

    return result;
}

void FlowArrowRenderer2D::onSliceAlignmentChange() {
    sliceAlignment_.getValue();
    switch (sliceAlignment_.getValue()) {
        case XY_PLANE:
            voxelPosPermutation_ = tgt::ivec3(0, 1, 2);
            break;
        case XZ_PLANE:
            voxelPosPermutation_ = tgt::ivec3(0, 2, 1);
            break;
        case YZ_PLANE:
            voxelPosPermutation_ = tgt::ivec3(2, 1, 0);
            break;
        default:
            break;
    }
    updatePropertyConfiguration();
}

void FlowArrowRenderer2D::mouseLocalization(tgt::MouseEvent* e) {

    if(e->getEventType() == tgt::MouseEvent::MOUSEPRESSEVENT)
        mouseIsPressed_ = true;
    else if(e->getEventType() == tgt::MouseEvent::MOUSERELEASEEVENT)
        mouseIsPressed_ = false;

    if (showCursorInfos_.isSelected("onMove") ||
        (showCursorInfos_.isSelected("onClick") && (e->getEventType() == tgt::MouseEvent::MOUSEPRESSEVENT || e->getEventType() == tgt::MouseEvent::WHEELEVENT)) ||
        (showCursorInfos_.isSelected("onDrag") && (e->getEventType() == tgt::MouseEvent::MOUSEMOVEEVENT && mouseIsPressed_))) {

            e->accept();
            if (mousePosition_ != e->coord()) {
                mousePosition_ = e->coord();

                tgt::ivec3 voxelPos = tgt::iround(screenToVoxelPos(mousePosition_));
                if (mousePosition_.x != -1 && voxelPos.x != -1) {
                    mouseXCoord_.set(voxelPos.x);
                    mouseYCoord_.set(voxelPos.y);
                    mouseZCoord_.set(voxelPos.z);
                }

                invalidate();
            }
    }
    // Don't save picking information in dragging mode when mouse is not pressed
    else if(e->getEventType() == tgt::MouseEvent::MOUSERELEASEEVENT && showCursorInfos_.isSelected("onDrag")) {
        e->accept();
        mousePosition_ = tgt::vec2(-1);
        lastPickingPosition_ = tgt::vec3(-1);
        invalidate();
    }
    else
        e->ignore();
}

void FlowArrowRenderer2D::shiftEvent(tgt::MouseEvent* e) {

    e->ignore();
    if (!velocityVolumeInport_.isReady() || !renderOutport_.isReady())
        return;

    if (e->action() == tgt::MouseEvent::PRESSED) {
        mousePosition_ = e->coord();
        return;
    }

    tgt::vec3 volDim = tgt::vec3(velocityVolumeInport_.getData()->getDimensions()) - 1.f;
    tgt::vec2 mouseCoords((float)e->coord().x, (float)e->coord().y);

    tgt::vec2 mouseOffset = mouseCoords - tgt::vec2(mousePosition_);
    mouseOffset.y *= -1.f;
    tgt::vec2 voxelOffset = voxelOffset_.get() +
        (mouseOffset / sliceSize_) * tgt::vec2(volDim[voxelPosPermutation_.x], volDim[voxelPosPermutation_.y]);
    voxelOffset = tgt::clamp(voxelOffset, voxelOffset_.getMinValue(), voxelOffset_.getMaxValue());
    voxelOffset_.set(voxelOffset);

    mousePosition_ = e->coord();
    e->accept();
}

bool FlowArrowRenderer2D::singleSliceMode() const {
    return true;
}

} // namespace voreen
