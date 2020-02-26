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

#include "unalignedsliceviewer.h"

#include "tgt/tgt_math.h"
#include <sstream>

#include "tgt/gpucapabilities.h"
#include "tgt/glmath.h"
#include "tgt/font.h"
#include "tgt/tgt_gl.h"
#include "tgt/textureunit.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "voreen/core/datastructures/octree/volumeoctreebase.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumegl.h"

using tgt::TextureUnit;

namespace voreen {

const std::string UnalignedSliceViewer::loggerCat_("voreen.base.SliceViewer");
const std::string UnalignedSliceViewer::fontName_("Vera.ttf");

UnalignedSliceViewer::UnalignedSliceViewer()
    : VolumeRenderer()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , transferFunc1_("transferFunction", "Transfer Function")
    , transferFunc2_("transferFunction2", "Transfer Function 2")
    , transferFunc3_("transferFunction3", "Transfer Function 3")
    , transferFunc4_("transferFunction4", "Transfer Function 4")
    , sliceAlignment_("sliceAlignmentProp", "Slice Alignment")
    , planeNormal_("planeNormal", "Plane Normal", tgt::vec3(0, 1, 0), -tgt::vec3::one, tgt::vec3::one)
    , planeDistance_("planeDistance", "Plane Distance", 0, -1000.0f, 1000.0f)
    , sliceIndex_("sliceIndex", "Slice Number ", 0, 0, 10000)
    , numGridRows_("numSlicesPerRow", "Num Rows", 1, 1, 5)
    , numGridCols_("numSlicesPerCol", "Num Columns", 1, 1, 5)
    , selectCenterSliceOnInputChange_("selectCenterSliceOnInputChange", "Auto-select Center Slice", true)
    , mouseCoord_("sliceMousePos", "Mouse Position", -tgt::ivec3::one, -tgt::ivec3::one, tgt::ivec3(10000))
    , renderSliceBoundaries_("renderSliceBoundaries", "Render Slice Boundaries", true)
    , boundaryColor_("boundaryColor", "Boundary Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , boundaryWidth_("boundaryWidth", "Boundary Width", 1, 1, 10)
    , texMode_("textureMode", "Texture Mode", Processor::INVALID_PROGRAM)
    , sliceLevelOfDetail_("sliceLevelOfDetail", "Slice Level of Detail", 0, 0, 3)
    , interactionLevelOfDetail_("interactionLevelOfDetail", "Interaction Level of Detail", 1, 0, 3, Processor::VALID)
    , sliceExtractionTimeLimit_("sliceExtracionTimeLimit", "Slice Creation Time Limit (ms)", 200, 0, 1000, Processor::VALID)
    , sliceCacheSize_("sliceCacheSize", "Slice Cache Size", 10, 0, 100, Processor::VALID)
    , showCursorInfos_("showCursorInformation", "Show Cursor Info")
    , showSliceNumber_("showSliceNumber", "Show Slice Number", true)
    , infoAlignment_("legendAlignment", "Legend Alignment")
    , showScaleLegend_("renderLegend", "Show Scale Legend", false)
    , legendLineLength_("legendLineLength", "Scale Legend Size (Pixels)", 80.f, 20.f, 250.f)
    , fontSize_("fontSize", "Font Size", 14, 8, 48)
    , voxelOffset_("voxelOffset", "Voxel Offset", tgt::vec2(0.f), tgt::vec2(std::numeric_limits<float>::lowest()), tgt::vec2(std::numeric_limits<float>::max()))
    , zoomFactor_("zoomFactor", "Zoom Factor", 1.f, 0.01f, 1.f)
    , resetViewButton_("resetViewButton", "Reset Shift and Zoom")
    , applyChannelShift_("applyChannelShift", "Apply Channel Shift", false, Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , channelShift1_("channelShift0", "Channel Shift 1", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift2_("channelShift1", "Channel Shift 2", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift3_("channelShift2", "Channel Shift 3", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift4_("channelShift3", "Channel Shift 4", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , resetChannelShift_("resetchannelshift", "Reset Channel Shift")
    , pickingMatrix_("pickingMatrix", "Picking Matrix", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f), Processor::VALID)
    , mwheelCycleHandler_("mouseWheelHandler", "Slice Cycling", &sliceIndex_)
    , mwheelZoomHandler_("zoomHandler", "Slice Zoom", &zoomFactor_, tgt::MouseEvent::CTRL)
    , sliceShader_(0)
    , sliceCache_(this, 10)
    , voxelPosPermutation_(0, 1, 2)
    , sliceLowerLeft_(1.f)
    , sliceSize_(0.f)
    , mousePosition_(-1, -1)
    , mouseIsPressed_(false)
    , lastPickingPosition_(-1, -1, -1)
    , sliceComplete_(true)
    , planeNeedsUpdate_(true)
{
    // texture mode (2D/3D)
    texMode_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode_.selectByKey("3d-texture");
    addProperty(texMode_);
    texMode_.onChange(MemberFunctionCallback<UnalignedSliceViewer>(this, &UnalignedSliceViewer::updatePropertyConfiguration));
    //texMode_.onChange(MemberFunctionCallback<UnalignedSliceViewer>(this, &UnalignedSliceViewer::rebuildShader()));  // does not work because of return value
    texMode_.setGroupID(inport_.getID());
    sliceLevelOfDetail_.setGroupID(inport_.getID());
    addProperty(sliceLevelOfDetail_);
    interactionLevelOfDetail_.setGroupID(inport_.getID());
    addProperty(interactionLevelOfDetail_);
    sliceExtractionTimeLimit_.setGroupID(inport_.getID());
    addProperty(sliceExtractionTimeLimit_);
    sliceCacheSize_.setGroupID(inport_.getID());
    sliceCacheSize_.onChange(MemberFunctionCallback<UnalignedSliceViewer>(this, &UnalignedSliceViewer::updatePropertyConfiguration));
    addProperty(sliceCacheSize_);

    inport_.addCondition(new PortConditionVolumeTypeGL());
    inport_.showTextureAccessProperties(true);
    addPort(inport_);
    addPort(outport_);

    setPropertyGroupGuiName(inport_.getID(), "Slice Technical Properties");

    // interaction
    mouseEventShift_ = new EventProperty<UnalignedSliceViewer>("mouseEvent.Shift", "Slice Shift",
        this, &UnalignedSliceViewer::shiftEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::MOTION, tgt::Event::CTRL);

    mouseEventMove_ = new EventProperty<UnalignedSliceViewer>("mouseEvent.cursorPositionMove", "Cursor Position Move",
        this, &UnalignedSliceViewer::mouseLocalization,
        tgt::MouseEvent::MOUSE_BUTTON_NONE,
        tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE,
        true);

    mouseEventPress_ = new EventProperty<UnalignedSliceViewer>("mouseEvent.cursorPositionPress", "Cursor Position Press",
        this, &UnalignedSliceViewer::mouseLocalization,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_MIDDLE),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED | tgt::MouseEvent::WHEEL | tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE,
        true);

    addEventProperty(mouseEventPress_);
    addEventProperty(mouseEventMove_);
    addEventProperty(mouseEventShift_);

    addInteractionHandler(mwheelCycleHandler_);
    addInteractionHandler(mwheelZoomHandler_);

    addProperty(transferFunc1_);
    addProperty(transferFunc2_);
    addProperty(transferFunc3_);
    addProperty(transferFunc4_);

    // slice arrangement
    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    sliceAlignment_.addOption("unaligned-plane", "Unaligned plane", UNALIGNED_PLANE);
    sliceAlignment_.onChange(
        MemberFunctionCallback<UnalignedSliceViewer>(this, &UnalignedSliceViewer::onSliceAlignmentChange) );
    addProperty(sliceAlignment_);

    addProperty(planeNormal_);
    ON_CHANGE_LAMBDA(planeNormal_, [this]{ planeNeedsUpdate_ = true; })
    addProperty(planeDistance_);
    ON_CHANGE_LAMBDA(planeDistance_, [this]{ planeNeedsUpdate_ = true; })

    addProperty(sliceIndex_);
    addProperty(numGridRows_);
    addProperty(numGridCols_);
    addProperty(selectCenterSliceOnInputChange_);

    addProperty(renderSliceBoundaries_);
    addProperty(boundaryColor_);
    addProperty(boundaryWidth_);
    addProperty(mouseCoord_);

    // group slice arrangement properties
    sliceAlignment_.setGroupID("sliceArrangement");
    planeNormal_.setGroupID("sliceArrangement");
    planeDistance_.setGroupID("sliceArrangement");
    sliceIndex_.setGroupID("sliceArrangement");
    numGridRows_.setGroupID("sliceArrangement");
    numGridCols_.setGroupID("sliceArrangement");
    selectCenterSliceOnInputChange_.setGroupID("sliceArrangement");
    renderSliceBoundaries_.setGroupID("sliceArrangement");
    boundaryColor_.setGroupID("sliceArrangement");
    boundaryWidth_.setGroupID("sliceArrangement");
    mouseCoord_.setGroupID("sliceArrangement");
    mouseCoord_.setVisibleFlag(false);
    mouseCoord_.setReadOnlyFlag(true);

    setPropertyGroupGuiName("sliceArrangement", "Slice Selection");

    // information overlay
    showCursorInfos_.addOption("never", "Never");
    showCursorInfos_.addOption("onClick", "On Mouse Click");
    showCursorInfos_.addOption("onDrag", "On Mouse Drag");
    showCursorInfos_.addOption("onMove", "On Mouse Move");
    showCursorInfos_.select("onMove");
    addProperty(showCursorInfos_);
    infoAlignment_.addOption("NW", "North-West", ALIGNMENT_NW);
    infoAlignment_.addOption("N", "North", ALIGNMENT_N);
    infoAlignment_.addOption("NE", "North-East", ALIGNMENT_NE);
    infoAlignment_.addOption("E", "East", ALIGNMENT_E);
    infoAlignment_.addOption("SE", "South-East", ALIGNMENT_SE);
    infoAlignment_.addOption("S", "South", ALIGNMENT_S);
    infoAlignment_.addOption("SW", "South-West", ALIGNMENT_SW);
    infoAlignment_.addOption("W", "West", ALIGNMENT_W);
    infoAlignment_.select("NW");
    addProperty(infoAlignment_);
    addProperty(showSliceNumber_);
    addProperty(showScaleLegend_);
    legendLineLength_.setStepping(1.f);
    addProperty(legendLineLength_);
    addProperty(fontSize_);

    // group information overlay properties
    showCursorInfos_.setGroupID("informationOverlay");
    infoAlignment_.setGroupID("informationOverlay");
    showSliceNumber_.setGroupID("informationOverlay");
    fontSize_.setGroupID("informationOverlay");
    showScaleLegend_.setGroupID("informationOverlay");
    legendLineLength_.setGroupID("informationOverlay");
    setPropertyGroupGuiName("informationOverlay", "Information Overlay");

    // zooming
    addProperty(voxelOffset_);
    zoomFactor_.setStepping(0.01f);
    addProperty(zoomFactor_);
    addProperty(resetViewButton_);
    pickingMatrix_.setReadOnlyFlag(true);
    addProperty(pickingMatrix_);

    // group zooming props
    voxelOffset_.setGroupID("zooming");
    zoomFactor_.setGroupID("zooming");
    resetViewButton_.setGroupID("zooming");
    pickingMatrix_.setGroupID("zooming");
    setPropertyGroupGuiName("zooming", "Zooming");

    resetViewButton_.onChange(MemberFunctionCallback<UnalignedSliceViewer>(this, &UnalignedSliceViewer::resetView));

    // channel shift
    addProperty(applyChannelShift_);
    addProperty(channelShift1_);
    addProperty(channelShift2_);
    addProperty(channelShift3_);
    addProperty(channelShift4_);
    addProperty(resetChannelShift_);
    resetChannelShift_.onClick(MemberFunctionCallback<UnalignedSliceViewer>(this, &UnalignedSliceViewer::resetChannelShift));
    applyChannelShift_.onChange(MemberFunctionCallback<UnalignedSliceViewer>(this, &UnalignedSliceViewer::updatePropertyConfiguration));
    applyChannelShift_.setGroupID("channel-shift");
    channelShift1_.setGroupID("channel-shift");
    channelShift2_.setGroupID("channel-shift");
    channelShift3_.setGroupID("channel-shift");
    channelShift4_.setGroupID("channel-shift");
    resetChannelShift_.setGroupID("channel-shift");
    setPropertyGroupGuiName("channel-shift", "Channel Shift Correction (in voxels)");

    // call this method to set the correct permutation for the
    // screen-position-to-voxel-position mapping.
    onSliceAlignmentChange();
}

UnalignedSliceViewer::~UnalignedSliceViewer() {
    delete mouseEventPress_;
    delete mouseEventMove_;
    delete mouseEventShift_;
}

Processor* UnalignedSliceViewer::create() const {
    return new UnalignedSliceViewer();
}

void UnalignedSliceViewer::initialize() {
    VolumeRenderer::initialize();

    sliceShader_ = ShdrMgr.load("sl_base", generateHeader(), false);
    LGL_ERROR;

    QualityMode.addObserver(this);

    updatePropertyConfiguration();
}

void UnalignedSliceViewer::deinitialize() {
    ShdrMgr.dispose(sliceShader_);
    sliceShader_ = 0;

    sliceCache_.clear();

    QualityMode.removeObserver(this);

    VolumeRenderer::deinitialize();
}

void UnalignedSliceViewer::adjustPropertiesToInput() {
    const VolumeBase* inputVolume = inport_.getData();

    if (selectCenterSliceOnInputChange_.get() && !firstProcessAfterDeserialization()/* && sliceIndex_.get() == 0*/) {
        if (inputVolume && sliceAlignment_.getValue() != UNALIGNED_PLANE) {
            int alignmentIndex = sliceAlignment_.getValue();
            tgtAssert(alignmentIndex >= 0 && alignmentIndex <= 2, "invalid alignment index");
            int centerSlice = (int)inputVolume->getDimensions()[alignmentIndex] / 2;
            sliceIndex_.set(centerSlice);
        }
    }
}

void UnalignedSliceViewer::beforeProcess() {
    VolumeRenderer::beforeProcess();

    if (inport_.hasChanged()) {
        mousePosition_ = tgt::ivec2(-1);
        lastPickingPosition_ = tgt::ivec3(-1);
        updatePropertyConfiguration();

        const VolumeBase* inputVolume = inport_.getData();
        transferFunc1_.setVolume(inputVolume, 0);
        if (inputVolume->getNumChannels() > 1)
            transferFunc2_.setVolume(inputVolume, 1);
        if (inputVolume->getNumChannels() > 2)
            transferFunc3_.setVolume(inputVolume, 2);
        if (inputVolume->getNumChannels() > 3)
            transferFunc4_.setVolume(inputVolume, 3);
    }

    if (invalidationLevel_ >= Processor::INVALID_PROGRAM || inport_.hasChanged())
        rebuildShader();
}

void UnalignedSliceViewer::afterProcess() {
    VolumeRenderer::afterProcess();

    if(QualityMode.isInteractionMode())
        processedInInteraction_ = true;

    if (!sliceComplete_)
        invalidate();
}

void UnalignedSliceViewer::qualityModeChanged() {
    if (!QualityMode.isInteractionMode() && processedInInteraction_) {
        processedInInteraction_ = false;
        invalidate();
    }
}

void UnalignedSliceViewer::updatePropertyConfiguration() {

    // properties not depending on the input volume
    bool sliceMode2D = texMode_.isSelected("2d-texture");
    sliceLevelOfDetail_.setReadOnlyFlag(!sliceMode2D);
    interactionLevelOfDetail_.setReadOnlyFlag(!sliceMode2D);
    sliceExtractionTimeLimit_.setReadOnlyFlag(!sliceMode2D);
    sliceCacheSize_.setReadOnlyFlag(!sliceMode2D);

    if (static_cast<size_t>(sliceCacheSize_.get()) != sliceCache_.getCacheSize())
        sliceCache_.setCacheSize(static_cast<size_t>(sliceCacheSize_.get()));

    // input-dependent properties
    if (!inport_.hasData()) {
        transferFunc2_.setVisibleFlag(false);
        transferFunc3_.setVisibleFlag(false);
        transferFunc4_.setVisibleFlag(false);
        return;
    }

    tgt::svec3 volumeDim = inport_.getData()->getDimensions();
    size_t numChannels = inport_.getData()->getNumChannels();

    bool unalignedPlane = sliceAlignment_.getValue() == UNALIGNED_PLANE;
    planeNormal_.setVisibleFlag(unalignedPlane);
    planeDistance_.setVisibleFlag(unalignedPlane);
    sliceIndex_.setVisibleFlag(!unalignedPlane);
    numGridCols_.setVisibleFlag(!unalignedPlane);
    numGridRows_.setVisibleFlag(!unalignedPlane);
    selectCenterSliceOnInputChange_.setVisibleFlag(!unalignedPlane);

    tgtAssert(sliceAlignment_.getValue() >= 0 && sliceAlignment_.getValue() <= 3, "Invalid alignment value");
    size_t numSlices = unalignedPlane ? 1 : volumeDim[sliceAlignment_.getValue()];
    if (numSlices == 0)
        return;

    sliceIndex_.setMaxValue((int)numSlices-1);
    if (sliceIndex_.get() >= static_cast<int>(numSlices))
        sliceIndex_.set(static_cast<int>(numSlices / 2));

    mouseCoord_.setMaxValue(volumeDim);

    numGridCols_.setMaxValue((int)numSlices);
    numGridRows_.setMaxValue((int)numSlices);

    if (numGridRows_.get() >= static_cast<int>(numSlices))
        numGridRows_.set((int)numSlices);

    if (numGridCols_.get() >= static_cast<int>(numSlices))
        numGridCols_.set((int)numSlices);

    tgt::vec2 halfDim = tgt::vec2((float)volumeDim[voxelPosPermutation_.x], (float)volumeDim[voxelPosPermutation_.y] - 1.f) / 2.f;
    voxelOffset_.setMinValue(-halfDim);
    voxelOffset_.setMaxValue(halfDim);
    voxelOffset_.set(tgt::clamp(voxelOffset_.get(), voxelOffset_.getMinValue(), voxelOffset_.getMaxValue()));

    // transfer functions
    transferFunc2_.setVisibleFlag(numChannels > 1);
    transferFunc3_.setVisibleFlag(numChannels > 2);
    transferFunc4_.setVisibleFlag(numChannels > 3);

    // channel shift
    channelShift2_.setVisibleFlag(numChannels > 1);
    channelShift3_.setVisibleFlag(numChannels > 2);
    channelShift4_.setVisibleFlag(numChannels > 3);

    bool channelShift = applyChannelShift_.get();
    channelShift1_.setReadOnlyFlag(!channelShift);
    channelShift2_.setReadOnlyFlag(!channelShift);
    channelShift3_.setReadOnlyFlag(!channelShift);
    channelShift4_.setReadOnlyFlag(!channelShift);
    resetChannelShift_.setReadOnlyFlag(!channelShift);
}

void UnalignedSliceViewer::resetChannelShift() {
    channelShift1_.set(tgt::vec3::zero);
    channelShift2_.set(tgt::vec3::zero);
    channelShift3_.set(tgt::vec3::zero);
    channelShift4_.set(tgt::vec3::zero);
}

void UnalignedSliceViewer::process() {
    const VolumeBase* volume = inport_.getData();

    // make sure VolumeGL is available in 3D texture mode
    if (texMode_.isSelected("3d-texture")) {
        if (sliceAlignment_.getValue() == UNALIGNED_PLANE || !volume->getRepresentation<VolumeGL>()) {
            LERROR("3D texture could not be created. Falling back to 2D texture mode.");
            texMode_.select("2d-texture");
            rebuildShader();
            return;
        }
    }

    sliceComplete_ = true;

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    const SliceAlignment alignment = sliceAlignment_.getValue();
    // get voxel volume dimensions
    tgt::ivec3 volDim = volume->getDimensions();
    int numSlices = alignment == UNALIGNED_PLANE ? 1 : volDim[alignment];

    // get the textures dimensions
    tgt::vec3 urf = volume->getURB();
    tgt::vec3 llb = volume->getLLF();

    // Re-calculate texture dimensions, urf, llb and center of the texture
    // for it might be a NPOT texture and therefore might have been inflated.
    // In that case, the inflating need to be undone and the volume texture
    // needs to be "cropped" to its original measures.
    const tgt::vec3 texDim = urf - llb;

    // Use OpenGL's ability of multiplying texture coordinate vectors with a matrix
    // on the texture matrix stack to permute the components of the texture
    // coordinates to obtain a correct setting for the current slice alignment.
    textureMatrix_ = tgt::mat4::zero;
    tgt::vec2 texDim2D(0.0f);   // laziness... (no matrix multiplication to determine 2D size of texture)

    tgt::mat4 toSliceCoordMatrix = tgt::mat4::zero;

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

            toSliceCoordMatrix.t00 = 1.f;
            toSliceCoordMatrix.t12 = 1.f;
            break;

        case YZ_PLANE:
            texDim2D.x = texDim.y;
            texDim2D.y = texDim.z;
            textureMatrix_.t02 = 1.f;    // setup a permutation matrix, swapping x-, y- and z-
            textureMatrix_.t10 = -1.f;   // components on vectors being multiplied with it
            textureMatrix_.t21 = 1.f;
            textureMatrix_.t13 = 1.f;
            textureMatrix_.t33 = 1.f;

            toSliceCoordMatrix.t01 = 1.f;
            toSliceCoordMatrix.t12 = 1.f;
            break;

        case XY_PLANE:
            texDim2D.x = texDim.x;
            texDim2D.y = texDim.y;
            textureMatrix_ = tgt::mat4::identity;
            textureMatrix_.t11 = -1.f;  // invert y-axis, since
            textureMatrix_.t13 = 1.f;   // we want to look along positive z-axis
            textureMatrix_.t33 = 1.f;

            toSliceCoordMatrix.t00 = 1.f;
            toSliceCoordMatrix.t11 = 1.f;
            break;

        case UNALIGNED_PLANE:
            texDim2D.x = resolution_.x;
            texDim2D.y = resolution_.y;
            textureMatrix_ = tgt::mat4::identity;

            toSliceCoordMatrix.t00 = 1.f;
            toSliceCoordMatrix.t11 = 1.f;
            break;
        default:
            break;
    }   // switch

    float canvasWidth = static_cast<float>(outport_.getSize().x);
    float canvasHeight = static_cast<float>(outport_.getSize().y);
    int numSlicesCol = numGridCols_.get();
    int numSlicesRow = numGridRows_.get();
    float scaleWidth = canvasWidth / (texDim2D.x * numSlicesCol);
    float scaleHeight = canvasHeight / (texDim2D.y * numSlicesRow);

    // find minimal scaling factor (either scale along canvas' width or
    // canvas' height)
    if (scaleWidth <= scaleHeight) {

        sliceSize_ = texDim2D * scaleWidth * (1.f/zoomFactor_.get());

        sliceLowerLeft_.x = 0;
        sliceLowerLeft_.y = (canvasHeight - (numSlicesRow * sliceSize_.y)) / 2.f;

        // adapt for zooming
        tgt::vec2 zoomOffset = -((sliceSize_ - (texDim2D * scaleWidth)) / 2.f) * (float)numSlicesCol;
        tgt::vec2 volDimFloat = tgt::vec2((float)volDim[voxelPosPermutation_.x], (float)volDim[voxelPosPermutation_.y]) - 1.f;
        tgt::vec2 focusOffset = voxelOffset_.get() / volDimFloat;
        focusOffset *= sliceSize_;
        sliceLowerLeft_.x += zoomOffset.x + focusOffset.x;
        sliceLowerLeft_.y += focusOffset.y;
    }
    else {
        sliceSize_ = texDim2D * scaleHeight * (1.f/zoomFactor_.get());

        sliceLowerLeft_.x = (canvasWidth - (numSlicesCol * sliceSize_.x)) / 2.f;
        sliceLowerLeft_.y = 0.f;

        // adapt for zooming
        tgt::vec2 zoomOffset = -((sliceSize_ - (texDim2D * scaleHeight)) / 2.f) * (float)numSlicesRow;
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
    MatStack.multMatrix(tgt::mat4::createOrtho(0.0f, canvasWidth, 0.f, canvasHeight, -1.0f, 1.0f));

    LGL_ERROR;

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    LGL_ERROR;

    // setup shader
    sliceShader_->activate();
    LGL_ERROR;

    // bind volume/slice texture
    bool setupSuccessful = true;
    TextureUnit texUnit;
    if (texMode_.isSelected("3d-texture")) { // 3D texture
        // bind volume
        std::vector<VolumeStruct> volumeTextures;
        volumeTextures.push_back(VolumeStruct(
            volume,
            &texUnit,
            "volume_","volumeParams_",
            inport_.getTextureClampModeProperty().getValue(),
            tgt::vec4(inport_.getTextureBorderIntensityProperty().get()),
            inport_.getTextureFilterModeProperty().getValue())
        );
        setupSuccessful = bindVolumes(sliceShader_, volumeTextures, 0, lightPosition_.get());
        LGL_ERROR;
    }
    else if (!texMode_.isSelected("2d-texture")){
        LERROR("unknown texture mode: " << texMode_.get());
        setupSuccessful = false;
    }
    if (!setupSuccessful) {
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.popMatrix();
        outport_.deactivateTarget();
        return;
    }

    // bind transfer functions
    TextureUnit transferUnit1, transferUnit2, transferUnit3, transferUnit4;
    transferUnit1.activate();
    transferFunc1_.get()->getTexture()->bind();
    transferFunc1_.get()->setUniform(sliceShader_, "transFuncParams_", "transFuncTex_", transferUnit1.getUnitNumber());
    LGL_ERROR;
    if (volume->getNumChannels() > 1) {
        transferUnit2.activate();
        transferFunc2_.get()->getTexture()->bind();
        transferFunc2_.get()->setUniform(sliceShader_, "transFuncParams2_", "transFuncTex2_", transferUnit2.getUnitNumber());
        LGL_ERROR;
    }
    if (volume->getNumChannels() > 2) {
        transferUnit3.activate();
        transferFunc3_.get()->getTexture()->bind();
        transferFunc3_.get()->setUniform(sliceShader_, "transFuncParams3_", "transFuncTex3_", transferUnit3.getUnitNumber());
        LGL_ERROR;
    }
    if (volume->getNumChannels() > 3) {
        transferUnit4.activate();
        transferFunc4_.get()->getTexture()->bind();
        transferFunc4_.get()->setUniform(sliceShader_, "transFuncParams4_", "transFuncTex4_", transferUnit4.getUnitNumber());
        LGL_ERROR;
    }

    // channel shift
    if (applyChannelShift_.get()) {

        sliceShader_->setUniform("channelShift_", channelShift1_.get() / tgt::vec3(volDim));
        if (volume->getNumChannels() > 1) {
            sliceShader_->setUniform("channelShift2_", channelShift2_.get() / tgt::vec3(volDim));
        }
        if (volume->getNumChannels() > 2) {
            sliceShader_->setUniform("channelShift3_", channelShift3_.get() / tgt::vec3(volDim));
        }
        if (volume->getNumChannels() > 3) {
            sliceShader_->setUniform("channelShift4_", channelShift4_.get() / tgt::vec3(volDim));
        }
    }

    // render slice(s) as textured quads
    float depth = 0.0f;
    const size_t sliceIndex = static_cast<size_t>(sliceIndex_.get());
    for (int pos = 0, x = 0, y = 0; pos < (numSlicesCol * numSlicesRow);
        ++pos, x = pos % numSlicesCol, y = pos / numSlicesCol)
    {
        if(alignment != UNALIGNED_PLANE) {
            int sliceNumber = (pos + static_cast<const int>(sliceIndex));
            if (sliceNumber >= numSlices)
                break;

            // compute depth in texture coordinates and check if it is not below the first or above the last slice
            depth = (static_cast<float>(sliceNumber) + 0.5f) / static_cast<float>(volume->getDimensions()[alignment]);
            float minDepth = 0.5f / static_cast<float>(volume->getDimensions()[alignment]);
            float maxDepth = (static_cast<float>(volume->getDimensions()[alignment]) - 0.5f) / static_cast<float>(volume->getDimensions()[alignment]);
            if (depth < minDepth || depth > maxDepth)
                continue;
        }

        MatStack.loadIdentity();
        MatStack.translate(sliceLowerLeft_.x + (x * sliceSize_.x),
            sliceLowerLeft_.y + ((numSlicesRow - (y + 1)) * sliceSize_.y), 0.0f);
        MatStack.scale(sliceSize_.x, sliceSize_.y, 1.0f);

        if(!sliceShader_->isActivated())
            sliceShader_->activate();

        setGlobalShaderParameters(sliceShader_);

        // extract 2D slice
        SliceTexture* slice = nullptr;
        if (texMode_.isSelected("2d-texture")) {

            bool singleSliceComplete = true;

            if(alignment != UNALIGNED_PLANE) {

                const size_t sliceID = tgt::iround(depth * volume->getDimensions()[alignment] - 0.5f);
                int* shiftArray = 0;
                if (applyChannelShift_.get()) {
                    //create shift array
                    shiftArray = new int[volume->getNumChannels()];
                    switch (volume->getNumChannels()) {
                        case 4:
                            shiftArray[3] = channelShift4_.get()[alignment];
                            //no break;
                        case 3:
                            shiftArray[2] = channelShift3_.get()[alignment];
                            //no break;
                        case 2:
                            shiftArray[1] = channelShift2_.get()[alignment];
                            //no break;
                        case 1:
                            shiftArray[0] = channelShift1_.get()[alignment];
                            break;
                        default:
                            tgtAssert(false, "unsupported channel count!");
                    }
                }
                switch (QualityMode.getQuality()) {
                    case VoreenQualityMode::RQ_INTERACTIVE:
                        slice = sliceCache_.getVolumeSlice(volume, alignment, sliceID, shiftArray,
                                                           interactionLevelOfDetail_.get(),
                                                           static_cast<clock_t>(sliceExtractionTimeLimit_.get()),
                                                           &singleSliceComplete, false);
                        break;
                    case VoreenQualityMode::RQ_DEFAULT:
                        slice = sliceCache_.getVolumeSlice(volume, alignment, sliceID, shiftArray,
                                                           sliceLevelOfDetail_.get(),
                                                           static_cast<clock_t>(sliceExtractionTimeLimit_.get()),
                                                           &singleSliceComplete, false);
                        break;
                    case VoreenQualityMode::RQ_HIGH:
                        //no time limit, octree level 0
                        slice = sliceCache_.getVolumeSlice(volume, alignment, sliceID, shiftArray, 0,
                                                           static_cast<clock_t>(0), &singleSliceComplete, false);
                        break;
                    default:
                        tgtAssert(false, "unknown rendering quality");
                }

                delete[] shiftArray;
            }
            else {

                // Update plane, if necessary.
                //if(planeNeedsUpdate_) //TODO: find out why this needs to be done each call of process().
                {
                    updatePlane();
                    planeNeedsUpdate_ = false;
                }

                switch (QualityMode.getQuality()) {
                    case VoreenQualityMode::RQ_INTERACTIVE:
                        samplingRate_ = 1.0f / (1 << interactionLevelOfDetail_.get());
                        break;
                    case VoreenQualityMode::RQ_DEFAULT:
                        samplingRate_ = 1.0f / (1 << sliceLevelOfDetail_.get());
                        break;
                    case VoreenQualityMode::RQ_HIGH:
                        samplingRate_ = 2.0f;
                        break;
                    default:
                        tgtAssert(false, "unknown rendering quality");
                }
                slice = sliceCache_.getVolumeSlice(volume, plane_, samplingRate_);
                singleSliceComplete = slice != nullptr;
            }

            sliceComplete_ &= singleSliceComplete;

            if (!slice)
                continue;

            // bind slice texture
            GLint texFilterMode = inport_.getTextureFilterModeProperty().getValue();
            GLint texClampMode = inport_.getTextureClampModeProperty().getValue();
            tgt::vec4 borderColor = tgt::vec4(inport_.getTextureBorderIntensityProperty().get());
            if (!GLSL::bindSliceTexture(slice, &texUnit, texFilterMode, texClampMode, borderColor))
                continue;

            // pass slice uniforms to shader
            sliceShader_->setIgnoreUniformLocationError(true);
            GLSL::setUniform(sliceShader_, "sliceTex_", "sliceTexParams_", slice, &texUnit);
            sliceShader_->setIgnoreUniformLocationError(false);


            sliceShader_->setUniform("toSliceCoordMatrix", toSliceCoordMatrix);

            LGL_ERROR;
        }

        sliceShader_->setUniform("textureMatrix_", textureMatrix_);
        tgt::vec2 texLowerLeft = tgt::vec2(0.f);
        tgt::vec2 texUpperRight = tgt::vec2(1.f);
        renderSliceGeometry(tgt::vec4(texLowerLeft.x, texLowerLeft.y, depth, 1.f),
            tgt::vec4(texUpperRight.x, texLowerLeft.y, depth, 1.f),
            tgt::vec4(texUpperRight.x, texUpperRight.y, depth, 1.f),
            tgt::vec4(texLowerLeft.x, texUpperRight.y, depth, 1.f));

    }   // textured slices

    tgtAssert(sliceShader_, "no slice shader");
    sliceShader_->deactivate();

    // render a border around each slice's boundaries if desired
    //
    if (renderSliceBoundaries_.get()) {
        MatStack.loadIdentity();
        MatStack.translate(sliceLowerLeft_.x, sliceLowerLeft_.y, 0.0f);
        MatStack.scale(sliceSize_.x * numSlicesCol, sliceSize_.y * numSlicesRow, 1.0f);
        glDepthFunc(GL_ALWAYS);
        renderSliceBoundaries();
        glDepthFunc(GL_LESS);
    }

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
    outport_.deactivateTarget();
    LGL_ERROR;

    // propagate picking matrix, if in single slice mode and a consumer is connected
    if (singleSliceMode() && !pickingMatrix_.getLinks().empty())
        pickingMatrix_.set(generatePickingMatrix());

    // propagate voxel position, if in show-cursor-info-on-move mode
    if (showCursorInfos_.isSelected("onMove")) {
        tgt::ivec3 voxelPos = tgt::iround(screenToVoxelPos(mousePosition_));
        mouseCoord_.set(voxelPos);
    }

}

void UnalignedSliceViewer::renderSliceGeometry(const tgt::vec4& t0, const tgt::vec4& t1, const tgt::vec4& t2, const tgt::vec4& t3) const {
    GLuint vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    std::vector<VertexHelper> buffer;
    buffer.push_back(VertexHelper(tgt::vec2(0.f, 0.f), t0));
    buffer.push_back(VertexHelper(tgt::vec2(1.f, 0.f), t1));
    buffer.push_back(VertexHelper(tgt::vec2(1.f, 1.f), t2));
    buffer.push_back(VertexHelper(tgt::vec2(0.f, 1.f), t3));

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, 4 * sizeof(VertexHelper), &(*buffer.begin()), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(VertexHelper), nullptr);
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(VertexHelper), (void*)sizeof(tgt::vec2));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(3);

    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

    glDeleteBuffers(1, &vbo);
    glBindVertexArray(0);
    glDeleteVertexArrays(1, &vao);
}

std::string UnalignedSliceViewer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = VolumeRenderer::generateHeader();

    if (texMode_.isSelected("2d-texture"))
        header += "#define SLICE_TEXTURE_MODE_2D \n";
    else if (texMode_.isSelected("3d-texture"))
        header += "#define SLICE_TEXTURE_MODE_3D \n";
    else {
        LWARNING("Unknown texture mode: " << texMode_.get());
    }

    header += "#define NUM_CHANNELS " + (inport_.hasData() ? itos(inport_.getData()->getNumChannels()) : "1") + " \n";

    header += transferFunc1_.get()->getShaderDefines();

    if (applyChannelShift_.get())
        header += "#define APPLY_CHANNEL_SHIFT\n";

    return header;
}

bool UnalignedSliceViewer::rebuildShader() {
    // do nothing if there is no shader at the moment
    if (!isInitialized() || !sliceShader_)
        return false;

    sliceShader_->setHeaders(generateHeader());
    return sliceShader_->rebuild();
}

void UnalignedSliceViewer::renderLegend() {

    tgtAssert(inport_.getData(), "No volume");
    const VolumeBase* volume = inport_.getData();

    // ESSENTIAL: if you don't use this, your text will become texturized!
    glActiveTexture(GL_TEXTURE0);
    glColor4f(1.f, 1.f, 1.f, 1.f);
    glDisable(GL_DEPTH_TEST);
    MatStack.loadIdentity();
    LGL_ERROR;


    tgt::Font legendFont(VoreenApplication::app()->getFontPath(fontName_));
    // note: the font size may not be smaller than 8
    legendFont.setFontSize(fontSize_.get());

    float legendOffsetX = 20.f;
    float legendOffsetY = 20.f;
    float lineLength = /*80.f;*/ legendLineLength_.get();
    float lineWidth = 5.f;
    tgt::ivec2 screensize = outport_.getSize();

    // determine scale String
    tgt::mat4 matrix = generatePickingMatrix();
    tgt::vec4 first = matrix*tgt::vec4(0.f,0.f,0.f,1.f);
    tgt::vec4 second = matrix*tgt::vec4(1.f,0.f,0.f,1.f);
    float scale = tgt::max(tgt::abs((first.xyz() - second.xyz())*volume->getSpacing()));
    std::string scaleStr = formatSpatialLength(scale*lineLength);

    // determine bounds and render the string
    tgt::vec2 size = legendFont.getSize(tgt::vec3(0.0f), scaleStr, screensize);
    float textLength = size.x;
    float textHeight = size.y;
    legendFont.render(tgt::vec3(std::max(legendOffsetX,legendOffsetX+(lineLength-textLength)/2.f), legendOffsetY, 0),
                      scaleStr, screensize);

    float lineOffsetY = legendOffsetY + textHeight + 5.f;
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
    glColor4f(1.f, 1.f, 1.f, 1.f);
    LGL_ERROR;

}

void UnalignedSliceViewer::renderSliceBoundaries() const {

    int numSlicesRow = numGridRows_.get();
    int numSlicesCol = numGridCols_.get();
    tgtAssert(numSlicesRow > 0 && numSlicesCol > 0, "Invalid slice counts");

    glLineWidth(static_cast<GLfloat>(boundaryWidth_.get()));
    glColor4f(boundaryColor_.get().r, boundaryColor_.get().g, boundaryColor_.get().b, boundaryColor_.get().a);
    glDisable(GL_DEPTH_TEST);
    glBegin(GL_LINE_LOOP);
        glVertex2f(0.0f, 0.0f);
        glVertex2f(1.0f, 0.0f);
        glVertex2f(1.0f, 1.0f);
        glVertex2f(0.0f, 1.0f);
    glEnd();

    float delta_x = 1.f / numSlicesCol;
    float delta_y = 1.f / numSlicesRow;
    glBegin(GL_LINES);
    for (int x = 1; x < numSlicesCol; ++x) {
        glVertex2f(delta_x * x, 0.f);
        glVertex2f(delta_x * x, 1.f);
    }

    for (int y = 1; y < numSlicesRow; ++y) {
        glVertex3f(0.0f, delta_y * y, 0.f);
        glVertex3f(1.0f, delta_y * y, 0.f);
    }
    glEnd();
    glEnable(GL_DEPTH_TEST);
    glColor4f(1.f, 1.f, 1.f, 1.f);
    glLineWidth(1.f);
    LGL_ERROR;
}

void UnalignedSliceViewer::renderInfoTexts() const {

    if (showCursorInfos_.isSelected("never") && !showSliceNumber_.get())
        return;

    int numSlicesRow = numGridRows_.get();
    int numSlicesCol = numGridCols_.get();
    tgtAssert(numSlicesRow > 0 && numSlicesCol > 0, "Invalid slice counts");

    tgtAssert(inport_.getData(), "No volume");
    const VolumeBase* volume = inport_.getData();
    tgt::ivec3 volDim = volume->getDimensions();
    int numSlices = sliceAlignment_.getValue() == UNALIGNED_PLANE ? 1 : volDim[sliceAlignment_.getValue()];

    glDisable(GL_DEPTH_TEST);

    // ESSENTIAL: if you don't use this, your text will become texturized!
    glActiveTexture(GL_TEXTURE0);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    LGL_ERROR;

    // render voxel position information
    if (!showCursorInfos_.isSelected("never")) {
        tgt::Font fontCursorInfos(VoreenApplication::app()->getFontPath(fontName_));
        // note: the font size may not be smaller than 8
        fontCursorInfos.setFontSize(fontSize_.get());

        // voxel position
        tgt::ivec3 voxelPos = tgt::iround(screenToVoxelPos(mousePosition_));
        if (mousePosition_.x != -1 && voxelPos.x != -1) { // save cursor position, if it is valid
            lastPickingPosition_ = voxelPos;
        }

        // render cursor information, if a valid picking position is available
        if (lastPickingPosition_.x != -1) {
            lastPickingPosition_ = tgt::clamp(lastPickingPosition_, tgt::ivec3(0), volDim-1);

            std::ostringstream oss;

            // voxel position string
            oss <<  "Voxel Pos:   " << "[" << lastPickingPosition_.x << " " << lastPickingPosition_.y << " " << lastPickingPosition_.z << "]";
            std::string voxelPosStr = oss.str();

            // spatial position
            oss.clear();
            oss.str("");
            oss << "Spatial Pos: " << formatSpatialLength(tgt::vec3(lastPickingPosition_)*volume->getSpacing());
            std::string spatialPosStr = oss.str();

            // data value
            oss.clear();
            oss.str("");
            oss << "Data Value: ";

            RealWorldMapping rwm = inport_.getData()->getRealWorldMapping();
            //get voxel without channel shift
            if(!applyChannelShift_.get()) {
                if (volume->hasRepresentation<VolumeRAM>()) {
                    const VolumeRAM* volume = inport_.getData()->getRepresentation<VolumeRAM>();
                    oss << volume->getVoxelValueAsString(lastPickingPosition_, &rwm);
                }
                else if (volume->hasRepresentation<VolumeOctreeBase>()) {
                    const VolumeOctreeBase* octree = volume->getRepresentation<VolumeOctreeBase>();
                    tgtAssert(octree, "no octree returned");
                    // apply real-world mapping to the normalized voxel value
                    size_t numChannels = octree->getNumChannels();
                    if (numChannels > 1)
                        oss << "[";
                    for (size_t i=0; i<numChannels; i++) {
                        oss << rwm.normalizedToRealWorld(static_cast<float>(octree->getVoxel(lastPickingPosition_, i)) / 65535.f);
                        if (i < numChannels-1)
                            oss << " ";
                    }
                    if (numChannels > 1)
                        oss << "]";
                    oss << " " << rwm.getUnit();
                }
                else if (volume->hasRepresentation<VolumeDisk>()) {
                    VolumeRAM* volume = 0;
                    try {
                        volume = inport_.getData()->getRepresentation<VolumeDisk>()->loadBrick(lastPickingPosition_,tgt::svec3::one);
                        oss << volume->getVoxelValueAsString(tgt::svec3::zero, &rwm);
                    }
                    catch (tgt::Exception&) {
                        // VolumeDisk does not support to load bricks -> try to load a slice instead
                        volume = inport_.getData()->getRepresentation<VolumeDisk>()->loadSlices(lastPickingPosition_.z, lastPickingPosition_.z);
                        oss << volume->getVoxelValueAsString(tgt::svec3(lastPickingPosition_.xy(), 0), &rwm);
                    }
                    delete volume;
                }
                else if (volume->hasRepresentation<VolumeGL>()) {
                    oss << "volume gl not supported yet";
                }
                else
                    oss << "unknown volume representation";
            } else {
                //voxel with channel shift
                size_t numChannels = volume->getNumChannels();
                tgt::ivec3 tmp;
                if (numChannels > 1)
                    oss << "[";
                for (size_t i=0; i<numChannels; i++) {
                        switch(i) {
                        case 0:
                            tmp = lastPickingPosition_ + tgt::ivec3(channelShift1_.get());
                            break;
                        case 1:
                            tmp = lastPickingPosition_ + tgt::ivec3(channelShift2_.get());
                            break;
                        case 2:
                            tmp = lastPickingPosition_ + tgt::ivec3(channelShift3_.get());
                            break;
                        case 3:
                            tmp = lastPickingPosition_ + tgt::ivec3(channelShift4_.get());
                            break;
                        default:
                            tgtAssert(false,"Should not get here");
                        }

                        if(tgt::max(tgt::lessThan(tmp,tgt::ivec3::zero) + tgt::lessThan(tgt::ivec3(volume->getDimensions())-tgt::ivec3::one,tmp)))
                            oss << rwm.normalizedToRealWorld(0.f);
                        else {
                            if (volume->hasRepresentation<VolumeRAM>()) {
                                const VolumeRAM* volume = inport_.getData()->getRepresentation<VolumeRAM>();
                                oss << volume->getVoxelValueAsString(tmp, &rwm,i);
                            }
                            else if (volume->hasRepresentation<VolumeOctreeBase>()) {
                                const VolumeOctreeBase* octree = volume->getRepresentation<VolumeOctreeBase>();
                                oss << rwm.normalizedToRealWorld(octree->getVoxel(tmp, i) / 65535.f);
                            }
                            else if (volume->hasRepresentation<VolumeDisk>()) {
                                VolumeRAM* volume = inport_.getData()->getRepresentation<VolumeDisk>()->loadBrick(tmp,tgt::svec3::one);
                                oss << volume->getVoxelValueAsString(tgt::svec3::zero, &rwm,i);
                                delete volume;
                            }
                            else if (volume->hasRepresentation<VolumeGL>()) {
                                oss << "volume gl not supported yet";
                            }
                            else
                                oss << "unknown volume representation";
                        }

                        //add empty space
                        if (i < numChannels-1)
                            oss << " ";
                    } // end for
                    //close line
                    if (numChannels > 1)
                        oss << "]";
                    oss << " " << rwm.getUnit();
            }

            std::string dataValueStr = oss.str();
            tgt::ivec2 screensize = outport_.getSize();

            std::string text = voxelPosStr+"\n"+spatialPosStr+"\n"+dataValueStr;

            fontCursorInfos.setTextAlignment(tgt::Font::BottomLeft);
            tgt::vec2 padding(3, 3);
            fontCursorInfos.setPadding(padding);
            tgt::vec2 boundingBox = fontCursorInfos.getSize(tgt::vec3::zero, text, screensize);

            // Calculate minimal width (using dummmy text) such that text not being aligned to the left won't change position everytime value changes.
            float minWidth = fontCursorInfos.getSize(tgt::vec3::zero, "Spatial Pos: [888.8, 888.8, 888.8] mm", screensize).x;
            if (boundingBox.x < minWidth)
                boundingBox.x = minWidth;

            tgt::vec3 offset = tgt::vec3::zero;
            switch (infoAlignment_.getValue()) {
            case ALIGNMENT_N:
                offset.x = (screensize.x - boundingBox.x - 2*padding.x) / 2;
                offset.y = 0;
                break;
            case ALIGNMENT_NE:
                offset.x = screensize.x - boundingBox.x - 2*padding.x;
                offset.y = 0;
                break;
            case ALIGNMENT_E:
                offset.x = screensize.x - boundingBox.x - 2*padding.x;
                offset.y = -(screensize.y - boundingBox.y - 2 * padding.y) / 2;
                break;
            case ALIGNMENT_SE:
                offset.x = screensize.x - boundingBox.x - 2*padding.x;
                offset.y = -(screensize.y - boundingBox.y - 2 * padding.y);
                break;
            case ALIGNMENT_S:
                offset.x = (screensize.x - boundingBox.x - 2*padding.x) / 2;
                offset.y = -(screensize.y - boundingBox.y - 2 * padding.y);
                break;
            case ALIGNMENT_SW:
                offset.x = 0;
                offset.y = -(screensize.y - boundingBox.y - 2 * padding.y);
                break;
            case ALIGNMENT_W:
                offset.x = 0;
                offset.y = -(screensize.y - boundingBox.y - 2 * padding.y) / 2;
                break;
            case ALIGNMENT_NW: // equal to zero
            default:
                break;
            }

            fontCursorInfos.render(offset+tgt::vec3(0, screensize.y, 0), text, screensize);

            MatStack.loadIdentity();
        }

        LGL_ERROR;
    }

    if (showSliceNumber_.get() && sliceAlignment_.getValue() != UNALIGNED_PLANE) {
        // note: the font size may not be smaller than 8
        tgt::Font fontSliceNumber(VoreenApplication::app()->getFontPath(fontName_));
        fontSliceNumber.setFontSize(fontSize_.get());

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


        tgt::ivec2 screensize = outport_.getSize();
        float textWidth = fontSliceNumber.getSize(tgt::vec3(6.f, 6.f, 0.f), dummy, screensize).x;

        // Do not render the slice numbers if the slice width becomes too small.
        if ((floorf(sliceSize_.x) > textWidth)) {
            // render the slice number information
            int sliceNumber = sliceIndex_.get();
            std::string prefix;
            for (int pos = 0, x = 0, y = 0; pos < numSlicesCol * numSlicesRow;
                ++pos, x = pos % numSlicesCol, y = pos / numSlicesCol)
            {
                MatStack.loadIdentity();
                if ((pos + sliceNumber) >= numSlices)
                    break;

                tgt::vec3 offset = tgt::vec3(sliceLowerLeft_.x + (x * sliceSize_.x) + sliceSize_.x,
                    sliceLowerLeft_.y + ((numSlicesRow - (y + 1)) * sliceSize_.y), 0.0f);

                std::ostringstream oss;
                if ((sliceNumber + pos) < 10)
                    oss << prefix10;
                else if ((sliceNumber + pos) < 100)
                    oss << prefix100;
                else if ((sliceNumber + pos) < 1000)
                    oss << prefix1000;

                oss << (sliceNumber + pos) << "/" << numSlices - 1;

                fontSliceNumber.setPadding(3);
                fontSliceNumber.setTextAlignment(tgt::Font::TopRight);;
                fontSliceNumber.render(offset, oss.str(), screensize);
            }
            LGL_ERROR;
        }

    }

    MatStack.loadIdentity();
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    LGL_ERROR;
}

tgt::vec3 UnalignedSliceViewer::screenToVoxelPos(tgt::ivec2 screenPos) const {

    const VolumeBase* volume = inport_.getData();

    if (!volume || !outport_.getRenderTarget())
        return tgt::vec3(-1.f);

    tgt::vec3 volumeDim(volume->getDimensions());
    tgt::ivec2 screenDim = outport_.getSize();

    tgt::ivec2 p(0, 0);
    p.x = screenPos.x - static_cast<int>(tgt::round(sliceLowerLeft_.x));
    p.y = (screenDim.y - screenPos.y) - static_cast<int>(tgt::round(sliceLowerLeft_.y));

    // if coordinates are negative, no slice could be hit
    if (tgt::hor(tgt::lessThan(p, tgt::ivec2(0))))
        return tgt::vec3(-1.f);

    const int numSlicesRow = numGridRows_.get();
    const int numSlicesCol = numGridCols_.get();

    const tgt::ivec2 sliceSizeInt = static_cast<tgt::ivec2>(sliceSize_);

    // if coordinates are greater than the number of slices per direction
    // times their extension in that direction, no slice could be hit either
    if ((p.x >= (sliceSizeInt.x * numSlicesCol)) || (p.y >= (sliceSizeInt.y * numSlicesRow)))
        return tgt::vec3(-1.0f);

    // calculate the normalized position within the picked slice
    tgt::vec2 posWithinSlice(
            static_cast<float>(p.x % sliceSizeInt.x),
            static_cast<float>(p.y % sliceSizeInt.y));
    posWithinSlice /= sliceSize_;

    tgt::vec3 voxPos;
    if (sliceAlignment_.getValue() == UNALIGNED_PLANE) {
        // TODO: use (picking) matrix for calculation!
        //voxPos = generatePickingMatrix() * (tgt::vec3(screenPos.x, screenPos.y, 0.0f));
        tgt::vec2 pos = posWithinSlice * tgt::vec2(resolution_);
        tgt::vec2 sp(tgt::min(volume->getSpacing()) / samplingRate_);
        tgt::vec3 fetchX = normalize(xVec_) * sp.x;
        tgt::vec3 fetchY = normalize(yVec_) * sp.y;
        tgt::vec3 fetchOrigin = origin_ + (0.5f * fetchX) + (0.5f * fetchY);
        tgt::mat4 wToV = volume->getWorldToVoxelMatrix();
        voxPos = fetchOrigin + (pos.x * fetchX) + (pos.y * fetchY);
        voxPos = wToV * voxPos;
    }
    else {

        // determine the picked slice
        const int sliceColID = p.x / sliceSizeInt.x;
        const int sliceRowID = (numSlicesRow - 1) - (p.y / sliceSizeInt.y);
        const int slice = sliceColID + (sliceRowID * numSlicesCol) + sliceIndex_.get();

        // calculate the normalized depth of the picked slice (texture z coordinate)
        float depth = (static_cast<float>(slice) + 0.5f) / std::max(volumeDim[voxelPosPermutation_.z], 1.f);

        // now we have the assigned texture coordinates of the picked fragment
        tgt::vec4 texCoords(posWithinSlice, depth, 1.f);
        texCoords = tgt::clamp(texCoords, tgt::vec4(0.f), tgt::vec4(1.f));

        // apply current texture matrix to assigned tex coords
        tgt::vec3 texCoordsTransformed = (textureMatrix_ * texCoords).xyz();

        // transform final tex coords into volume coordinates
        voxPos = texCoordsTransformed * (volumeDim) - tgt::vec3(0.5f);
    }

    if(tgt::clamp(voxPos, tgt::vec3::zero, volumeDim - tgt::vec3::one) != voxPos) {
        return -tgt::vec3::one;
    }
    return voxPos;
}

tgt::mat4 UnalignedSliceViewer::generatePickingMatrix() const {

    if (!inport_.hasData())
        return tgt::mat4::createIdentity();

    tgt::vec3 volumeDim(inport_.getData()->getDimensions());

    if(sliceAlignment_.getValue() != UNALIGNED_PLANE) {

        // 1. translate slice to origin (also add 0.5 in depth direction, since we subtract this afterwards and the slice number should already be correct)
        tgt::mat4 originTranslation = tgt::mat4::createTranslation(
                tgt::vec3(-sliceLowerLeft_.x, -sliceLowerLeft_.y, 0.5f));

        // 2. normalize screen coords with regard to the slice
        tgt::mat4 sliceScale = tgt::mat4::createScale(
                tgt::vec3(1.f / sliceSize_.x, 1.f / sliceSize_.y, 1.f / (volumeDim[voxelPosPermutation_.z])));

        // 3. apply current texture matrix
        //tgt::mat4 textureMatrix = textureMatrix_; // used directly (see below)

        // 4. transform normalized coordinates to volume dimensions
        tgt::mat4 volumeScale = tgt::mat4::createTranslation(tgt::vec3(-0.5f)) * tgt::mat4::createScale(volumeDim);

        // compose transformation matrix
        tgt::mat4 result = volumeScale * textureMatrix_ * sliceScale * originTranslation;

        return result;
    }
    else {

        // TODO: fix  the code below!
        const VolumeBase* volume = inport_.getData();

        tgt::plane plane = tgt::plane(planeNormal_.get(), planeDistance_.get());

        tgt::vec3 urb = volume->getURB();
        tgt::vec3 llf = volume->getLLF();
        tgt::vec3 center = (urb + llf) * 0.5f;

        tgt::vec3 xMax = center;
        xMax.x = urb.x;
        tgt::vec3 yMax = center;
        yMax.y = urb.y;
        tgt::vec3 zMax = center;
        zMax.z = urb.z;

        // transform to world coordinates:
        tgt::mat4 pToW = volume->getPhysicalToWorldMatrix();
        center = pToW * center;
        xMax = pToW * xMax;
        yMax = pToW * yMax;
        zMax = pToW * zMax;

        // project to plane:
        float d = plane.distance(center);
        center = center - (plane.n * d);
        d = plane.distance(xMax);
        xMax = xMax - (plane.n * d);
        d = plane.distance(yMax);
        yMax = yMax - (plane.n * d);
        d = plane.distance(zMax);
        zMax = zMax - (plane.n * d);

        // find max axis in plane:
        tgt::vec3 maxVec = xMax - center;
        if(distance(yMax, center) > length(maxVec))
            maxVec = yMax - center;
        if(distance(zMax, center) > length(maxVec))
            maxVec = zMax - center;

        maxVec = normalize(maxVec);
        tgt::vec3 temp = normalize(cross(maxVec, plane.n));

        // construct transformation to temporary system:
        tgt::mat4 m(maxVec.x, temp.x, plane.n.x, center.x,
                    maxVec.y, temp.y, plane.n.y, center.y,
                    maxVec.z, temp.z, plane.n.z, center.z,
                    0.0f,     0.0f,   0.0f,   1.0f);
        tgt::mat4 mInv = tgt::mat4::identity;
        m.invert(mInv);
        m = mInv * pToW;

        // 1. translate slice to origin (also add 0.5 in depth direction, since we subtract this afterwards and the slice number should already be correct)
        tgt::mat4 originTranslation = tgt::mat4::createTranslation(
                tgt::vec3(-sliceLowerLeft_.x, -sliceLowerLeft_.y, 0.5f));

        // 2. normalize screen coords with regard to the slice
        tgt::mat4 sliceScale = tgt::mat4::createScale(
                tgt::vec3(1.f / sliceSize_.x, 1.f / sliceSize_.y, 1.f / (volumeDim[voxelPosPermutation_.z])));

        // 3. apply current texture matrix
        //tgt::mat4 textureMatrix = textureMatrix_; // used directly (see below)

        // 4. transform normalized coordinates to volume dimensions
        tgt::mat4 volumeScale = tgt::mat4::createTranslation(tgt::vec3(-0.5f)) * tgt::mat4::createScale(volumeDim);

        // compose transformation matrix
        tgt::mat4 result = m * volumeScale * textureMatrix_ * sliceScale * originTranslation;

        return result;
    }
}

void UnalignedSliceViewer::updatePlane() {

    const VolumeBase* volume = inport_.getData();

    plane_ = tgt::plane(planeNormal_.get(), planeDistance_.get());

    tgt::vec3 urb = volume->getURB();
    tgt::vec3 llf = volume->getLLF();
    tgt::vec3 center = (urb + llf) * 0.5f;

    tgt::vec3 xMax = center;
    xMax.x = urb.x;
    tgt::vec3 yMax = center;
    yMax.y = urb.y;
    tgt::vec3 zMax = center;
    zMax.z = urb.z;

    // transform to world coordinates:
    tgt::mat4 pToW = volume->getPhysicalToWorldMatrix();
    center = pToW * center;
    xMax = pToW * xMax;
    yMax = pToW * yMax;
    zMax = pToW * zMax;

    // project to plane:
    float d = plane_.distance(center);
    center = center - (plane_.n * d);
    d = plane_.distance(xMax);
    xMax = xMax - (plane_.n * d);
    d = plane_.distance(yMax);
    yMax = yMax - (plane_.n * d);
    d = plane_.distance(zMax);
    zMax = zMax - (plane_.n * d);

    // find max axis in plane:
    tgt::vec3 maxVec = xMax - center;
    if(distance(yMax, center) > length(maxVec))
        maxVec = yMax - center;
    if(distance(zMax, center) > length(maxVec))
        maxVec = zMax - center;

    maxVec = normalize(maxVec);
    tgt::vec3 temp = normalize(cross(maxVec, plane_.n));

    // construct transformation to temporary system:
    tgt::mat4 m(maxVec.x, temp.x, plane_.n.x, center.x,
                maxVec.y, temp.y, plane_.n.y, center.y,
                maxVec.z, temp.z, plane_.n.z, center.z,
                0.0f,     0.0f,   0.0f,   1.0f);
    tgt::mat4 mInv = tgt::mat4::identity;
    m.invert(mInv);

    // transform bounds to temp system in order to construct new coordinate frame
    tgt::Bounds b(volume->getLLF(), volume->getURB());
    b = b.transform(mInv*pToW);

    // construct new coordinate frame:
    origin_ = center;
    origin_ += b.getLLF().x * maxVec;
    origin_ += b.getLLF().y * temp;

    tgt::vec2 sp(tgt::min(volume->getSpacing()) / samplingRate_);
    resolution_ = tgt::ivec2(tgt::iceil(b.diagonal().x / sp.x), tgt::iceil(b.diagonal().y / sp.y));

    xVec_ = maxVec * (sp.x * resolution_.x);
    yVec_ = temp * (sp.y * resolution_.y);
}

void UnalignedSliceViewer::onSliceAlignmentChange() {
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
        case UNALIGNED_PLANE:
            voxelPosPermutation_ = tgt::ivec3(0, 1, 2);
            break;
        default:
            break;
    }
    updatePropertyConfiguration();
}

void UnalignedSliceViewer::mouseLocalization(tgt::MouseEvent* e) {

    if(e->getEventType() == tgt::MouseEvent::MOUSEPRESSEVENT)
        mouseIsPressed_ = true;
    else if(e->getEventType() == tgt::MouseEvent::MOUSERELEASEEVENT)
        mouseIsPressed_ = false;

    //FIXME: HACK: to work properbly with point viewer
    if(!mouseCoord_.getLinks().empty()) {
        tgt::ivec3 tmpVoxelPos = tgt::iround(screenToVoxelPos(e->coord()));
        if (e->coord().x != -1 && tmpVoxelPos.x != -1) {
            mouseCoord_.set(tmpVoxelPos);
        }
    }


    if (showCursorInfos_.isSelected("onMove") ||
        (showCursorInfos_.isSelected("onClick") && (e->getEventType() == tgt::MouseEvent::MOUSEPRESSEVENT || e->getEventType() == tgt::MouseEvent::WHEELEVENT)) ||
        (showCursorInfos_.isSelected("onDrag") && (e->getEventType() == tgt::MouseEvent::MOUSEMOVEEVENT && mouseIsPressed_))) {

            e->accept();
            if (mousePosition_ != e->coord()) {
                mousePosition_ = e->coord();

                tgt::ivec3 voxelPos = tgt::iround(screenToVoxelPos(mousePosition_));
                if (mousePosition_.x != -1 && voxelPos.x != -1) {
                    mouseCoord_.set(voxelPos);
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

void UnalignedSliceViewer::shiftEvent(tgt::MouseEvent* e) {

    e->ignore();
    if (!inport_.isReady() || !outport_.isReady())
        return;

    if (e->action() == tgt::MouseEvent::PRESSED) {
        mousePosition_ = e->coord();
        return;
    }

    tgt::vec3 volDim = tgt::vec3(inport_.getData()->getDimensions()) - 1.f;
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

bool UnalignedSliceViewer::singleSliceMode() const {
    return (numGridRows_.get() == 1 && numGridCols_.get() == 1);
}

void UnalignedSliceViewer::resetView()
{
    zoomFactor_.reset();
    voxelOffset_.reset();
}

} // namespace voreen
