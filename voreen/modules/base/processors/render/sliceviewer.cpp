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

#include "sliceviewer.h"

#include "tgt/tgt_math.h"
#include <sstream>
#include <chrono>

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
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/datastructures/octree/volumeoctreebase.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumegl.h"

using tgt::TextureUnit;

namespace voreen {

OctreeSliceTexture::OctreeSliceTexture()
    : buf_()
    , texture_(tgt::svec3(2,2,1), GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, tgt::Texture::LINEAR, tgt::Texture::CLAMP_TO_EDGE, nullptr, false)
{
}

void OctreeSliceTexture::updateDimensions(tgt::svec2 dim) {
    if(dim != tgt::svec2(texture_.getDimensions().xy())) {
        texture_.setCpuTextureData(nullptr, false);
        texture_.updateDimensions(tgt::ivec3(dim.x, dim.y, 1), false);

        buf_.resize(tgt::hmul(dim), OctreeSliceTexture::Pixel(0));
        texture_.setCpuTextureData(reinterpret_cast<GLubyte*>(buf_.data()), false);
        LGL_ERROR;
    }
}

void OctreeSliceTexture::uploadTexture() {
    texture_.uploadTexture();
    LGL_ERROR;
}

void OctreeSliceTexture::bindTexture() {
    texture_.bind();
    LGL_ERROR;
}
OctreeSliceTexture::Pixel* OctreeSliceTexture::buf() {
    return buf_.data();
}
tgt::ivec2 OctreeSliceTexture::dimensions() {
    return texture_.getDimensions().xy();
}
void OctreeSliceTexture::clear() {
    std::fill(buf_.begin(), buf_.end(), OctreeSliceTexture::Pixel(0));
}
OctreeSliceViewProgress::OctreeSliceViewProgress()
    : nextSlice_(0)
    , nextTileX_(0)
    , nextTileY_(0)
{}

const std::string SliceViewer::loggerCat_("voreen.base.SliceViewer");
const std::string SliceViewer::fontName_("Vera.ttf");

SliceViewer::SliceViewer()
    : VolumeRenderer()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , transferFunc1_("transferFunction", "Transfer Function")
    , transferFunc2_("transferFunction2", "Transfer Function 2")
    , transferFunc3_("transferFunction3", "Transfer Function 3")
    , transferFunc4_("transferFunction4", "Transfer Function 4")
    , sliceAlignment_("sliceAlignmentProp", "Slice Alignment")
    , sliceIndex_("sliceIndex", "Slice Number ", 0, 0, 10000)
    , numGridRows_("numSlicesPerRow", "Num Rows", 1, 1, 5)
    , numGridCols_("numSlicesPerCol", "Num Columns", 1, 1, 5)
    , selectCenterSliceOnInputChange_("selectCenterSliceOnInputChange", "Auto-select Center Slice", true)
    , mouseCoord_("sliceMousePos", "Mouse Position", -tgt::ivec3::one, -tgt::ivec3::one, tgt::ivec3(100000))
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
    , copyImageShader_(0)
    , sliceCache_(this, 10)
    , voxelPosPermutation_(0, 1, 2)
    , sliceLowerLeft_(1.f)
    , sliceSize_(0.f)
    , mousePosition_(-1, -1)
    , mouseIsPressed_(false)
    , lastPickingPosition_(-1, -1, -1)
    , sliceComplete_(true)
    , octreeTexture_(nullptr)
    , octreeRenderProgress_()
{
    // texture mode (2D/3D)
    texMode_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode_.addOption("octree", "Octree", OCTREE);
    texMode_.selectByKey("3d-texture");
    addProperty(texMode_);
    ON_CHANGE(texMode_, SliceViewer, invalidateOctreeTexture);
    texMode_.onChange(MemberFunctionCallback<SliceViewer>(this, &SliceViewer::updatePropertyConfiguration));
    //texMode_.onChange(MemberFunctionCallback<SliceViewer>(this, &SliceViewer::rebuildShader()));  // does not work because of return value
    texMode_.setGroupID(inport_.getID());
    sliceLevelOfDetail_.setGroupID(inport_.getID());
    ON_CHANGE(sliceLevelOfDetail_, SliceViewer, invalidateOctreeTexture);
    addProperty(sliceLevelOfDetail_);
    ON_CHANGE(interactionLevelOfDetail_, SliceViewer, invalidateOctreeTexture);
    interactionLevelOfDetail_.setGroupID(inport_.getID());
    addProperty(interactionLevelOfDetail_);
    sliceExtractionTimeLimit_.setGroupID(inport_.getID());
    addProperty(sliceExtractionTimeLimit_);
    sliceCacheSize_.setGroupID(inport_.getID());
    sliceCacheSize_.onChange(MemberFunctionCallback<SliceViewer>(this, &SliceViewer::updatePropertyConfiguration));
    addProperty(sliceCacheSize_);

    inport_.addCondition(new PortConditionVolumeTypeGL());
    inport_.showTextureAccessProperties(true);
    addPort(inport_);
    outport_.onSizeReceiveChange<SliceViewer>(this, &SliceViewer::invalidateOctreeTexture);
    addPort(outport_);

    setPropertyGroupGuiName(inport_.getID(), "Slice Technical Properties");

    // interaction
    mouseEventShift_ = new EventProperty<SliceViewer>("mouseEvent.Shift", "Slice Shift",
        this, &SliceViewer::shiftEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::MOTION, tgt::Event::CTRL);

    mouseEventMove_ = new EventProperty<SliceViewer>("mouseEvent.cursorPositionMove", "Cursor Position Move",
        this, &SliceViewer::mouseLocalization,
        tgt::MouseEvent::MOUSE_BUTTON_NONE,
        tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE,
        true);

    mouseEventPress_ = new EventProperty<SliceViewer>("mouseEvent.cursorPositionPress", "Cursor Position Press",
        this, &SliceViewer::mouseLocalization,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_MIDDLE),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED | tgt::MouseEvent::WHEEL | tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE,
        true);

    addEventProperty(mouseEventPress_);
    addEventProperty(mouseEventMove_);
    addEventProperty(mouseEventShift_);

    addInteractionHandler(mwheelCycleHandler_);
    addInteractionHandler(mwheelZoomHandler_);

    addProperty(transferFunc1_);
    ON_CHANGE(transferFunc1_, SliceViewer, invalidateOctreeTexture);
    addProperty(transferFunc2_);
    ON_CHANGE(transferFunc2_, SliceViewer, invalidateOctreeTexture);
    addProperty(transferFunc3_);
    ON_CHANGE(transferFunc3_, SliceViewer, invalidateOctreeTexture);
    addProperty(transferFunc4_);
    ON_CHANGE(transferFunc4_, SliceViewer, invalidateOctreeTexture);

    // slice arrangement
    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    sliceAlignment_.onChange(
        MemberFunctionCallback<SliceViewer>(this, &SliceViewer::onSliceAlignmentChange) );
    ON_CHANGE(sliceAlignment_, SliceViewer, invalidateOctreeTexture);
    addProperty(sliceAlignment_);

    addProperty(sliceIndex_);
    ON_CHANGE(sliceIndex_, SliceViewer, invalidateOctreeTexture);
    addProperty(numGridRows_);
    ON_CHANGE(numGridRows_, SliceViewer, invalidateOctreeTexture);
    addProperty(numGridCols_);
    ON_CHANGE(numGridCols_, SliceViewer, invalidateOctreeTexture);
    addProperty(selectCenterSliceOnInputChange_);

    addProperty(renderSliceBoundaries_);
    addProperty(boundaryColor_);
    addProperty(boundaryWidth_);
    addProperty(mouseCoord_);

    // group slice arrangement properties
    sliceAlignment_.setGroupID("sliceArrangement");
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
    ON_CHANGE(voxelOffset_, SliceViewer, invalidateOctreeTexture);
    zoomFactor_.setStepping(0.01f);
    addProperty(zoomFactor_);
    ON_CHANGE(zoomFactor_, SliceViewer, invalidateOctreeTexture);
    addProperty(resetViewButton_);
    pickingMatrix_.setReadOnlyFlag(true);
    addProperty(pickingMatrix_);

    // group zooming props
    voxelOffset_.setGroupID("zooming");
    zoomFactor_.setGroupID("zooming");
    resetViewButton_.setGroupID("zooming");
    pickingMatrix_.setGroupID("zooming");
    setPropertyGroupGuiName("zooming", "Zooming");

    resetViewButton_.onChange(MemberFunctionCallback<SliceViewer>(this, &SliceViewer::resetView));

    // channel shift
    addProperty(applyChannelShift_);
    ON_CHANGE(applyChannelShift_, SliceViewer, invalidateOctreeTexture);
    addProperty(channelShift1_);
    ON_CHANGE(channelShift1_, SliceViewer, invalidateOctreeTexture);
    addProperty(channelShift2_);
    ON_CHANGE(channelShift2_, SliceViewer, invalidateOctreeTexture);
    addProperty(channelShift3_);
    ON_CHANGE(channelShift3_, SliceViewer, invalidateOctreeTexture);
    addProperty(channelShift4_);
    ON_CHANGE(channelShift4_, SliceViewer, invalidateOctreeTexture);
    addProperty(resetChannelShift_);
    resetChannelShift_.onClick(MemberFunctionCallback<SliceViewer>(this, &SliceViewer::resetChannelShift));
    applyChannelShift_.onChange(MemberFunctionCallback<SliceViewer>(this, &SliceViewer::updatePropertyConfiguration));
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

SliceViewer::~SliceViewer() {
    delete mouseEventPress_;
    delete mouseEventMove_;
    delete mouseEventShift_;
}

Processor* SliceViewer::create() const {
    return new SliceViewer();
}

void SliceViewer::initialize() {
    VolumeRenderer::initialize();

    copyImageShader_ = ShdrMgr.load("passthrough", RenderProcessor::generateHeader(), false);
    LGL_ERROR;
    sliceShader_ = ShdrMgr.load("sl_base", generateSliceShaderHeader(), false);
    LGL_ERROR;

    octreeTexture_.reset(new OctreeSliceTexture());

    QualityMode.addObserver(this);

    updatePropertyConfiguration();
}

void SliceViewer::deinitialize() {
    ShdrMgr.dispose(sliceShader_);
    sliceShader_ = 0;
    ShdrMgr.dispose(copyImageShader_);
    copyImageShader_ = 0;

    octreeTexture_.reset();

    sliceCache_.clear();

    QualityMode.removeObserver(this);

    VolumeRenderer::deinitialize();
}

void SliceViewer::adjustPropertiesToInput() {
    const VolumeBase* inputVolume = inport_.getData();

    if (selectCenterSliceOnInputChange_.get() && !firstProcessAfterDeserialization()/* && sliceIndex_.get() == 0*/) {
        if (inputVolume) {
            int alignmentIndex = sliceAlignment_.getValue();
            tgtAssert(alignmentIndex >= 0 && alignmentIndex <= 2, "invalid alignment index");
            int centerSlice = (int)inputVolume->getDimensions()[alignmentIndex] / 2;
            sliceIndex_.set(centerSlice);
        }
    }
}

void SliceViewer::invalidateOctreeTexture() {
    if(texMode_.getValue() == OCTREE && octreeTexture_) {
        octreeTexture_->clear();
        octreeRenderProgress_ = OctreeSliceViewProgress();
        invalidate();
    }
}

void SliceViewer::beforeProcess() {
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

void SliceViewer::afterProcess() {
    VolumeRenderer::afterProcess();

    if(QualityMode.isInteractionMode())
        processedInInteraction_ = true;

    if (!sliceComplete_) {
        invalidate();
    }
}

void SliceViewer::qualityModeChanged() {
    if (!QualityMode.isInteractionMode() && processedInInteraction_) {
        processedInInteraction_ = false;
        invalidate();
    }
}

void SliceViewer::updatePropertyConfiguration() {

    // properties not depending on the input volume
    bool mode2d = texMode_.isSelected("2d-texture");
    bool mode3d = texMode_.isSelected("3d-texture");
    bool modeoctree = texMode_.isSelected("octree");
    sliceLevelOfDetail_.setReadOnlyFlag(!(mode2d || modeoctree));
    interactionLevelOfDetail_.setReadOnlyFlag(!(mode2d || modeoctree));
    sliceExtractionTimeLimit_.setReadOnlyFlag(!(mode2d || modeoctree));
    sliceCacheSize_.setReadOnlyFlag(!mode2d);

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

    tgtAssert(sliceAlignment_.getValue() >= 0 && sliceAlignment_.getValue() <= 2, "Invalid alignment value");
    size_t numSlices = volumeDim[sliceAlignment_.getValue()];
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

void SliceViewer::resetChannelShift() {
    channelShift1_.set(tgt::vec3::zero);
    channelShift2_.set(tgt::vec3::zero);
    channelShift3_.set(tgt::vec3::zero);
    channelShift4_.set(tgt::vec3::zero);
}

static std::vector<tgt::Vector4<uint8_t>> extractLUT(const TransFunc1DKeys& tf) {
    tgt::Texture* texture = tf.getTexture();

    texture->downloadTexture();
    const int mapDim = texture->getDimensions().x;

    std::vector<tgt::Vector4<uint8_t>> valueMap;
    valueMap.reserve(mapDim);
    for(size_t i=0; i<mapDim; ++i) {
        valueMap.push_back(texture->texel<tgt::Vector4<uint8_t>>(i));
    }
    return valueMap;
}

struct TransferFuncCPU {
    std::vector<tgt::Vector4<uint8_t>> lut_;
    RealWorldMapping inputToLut_;

    TransferFuncCPU(TransFunc1DKeys* tf, RealWorldMapping volumeRwm)
        : lut_()
        , inputToLut_()
    {
        std::vector<uint8_t> valueMap;
        RealWorldMapping tfRwm;

        if(tf) {
            lut_ = extractLUT(*tf);
            tfRwm = RealWorldMapping(tf->getDomain(), "");
        } else {
            lut_ = extractLUT(TransFunc1DKeys());
            tfRwm = RealWorldMapping(tgt::vec2(0.0, 1.0), "");
        }

        inputToLut_ = RealWorldMapping::combine(volumeRwm, tfRwm.getInverseMapping());
    }

    tgt::Vector4<uint8_t> lookUp(float val) const {
        float lutVal = inputToLut_.normalizedToRealWorld(val);
        size_t size = lut_.size();
        size_t index = tgt::clamp<size_t>(tgt::round(lutVal * (size-1)), 0, size-1);
        return lut_[index];
    }
};

typedef std::chrono::steady_clock Clock;
typedef std::chrono::time_point<Clock> TimePoint;

enum class DeadlineResult {
    Succeeded,
    TimedOut,
};

static DeadlineResult renderOctreeSlice(OctreeSliceTexture& texture, const VolumeOctree& octree, tgt::ivec2 pixBegin, tgt::ivec2 pixEnd, tgt::mat4 pixelToVoxelMats[4], size_t level, const TransferFuncCPU tfs[4], OctreeSliceViewProgress& progress, TimePoint deadline) {
    int rowSize = texture.dimensions().x;
    auto* pixels = texture.buf();

    tgt::mat4 voxelToPixel;
    bool success = pixelToVoxelMats[0].invert(voxelToPixel);
    tgtAssert(success, "Failed to invert pixelToVoxel matrix");

    tgt::ivec3 volDim = octree.getDimensions();
    LocatedVolumeOctreeNodeConst root = octree.getLocatedRootNode();
    const OctreeBrickPoolManagerBase& brickPoolManager = *octree.getBrickPoolManager();
    const tgt::svec3 brickDataSize = octree.getBrickDim();
    const size_t numChannels = octree.getNumChannels();

    // Try to hit the same brick as many times as possible by choosing a tile
    // size that corresponds to one (or slightly less than one) brick...
    const tgt::svec3 brickSizeInVoxels = brickDataSize * (1UL << level);
    const tgt::ivec2 tileSizeRaw = tgt::abs(voxelToPixel*tgt::vec3::zero-voxelToPixel*tgt::vec3(brickSizeInVoxels)).xy();
    // ... but a tile size of one or smaller doesn't make sense either.
    const tgt::ivec2 tileSize = tgt::max(tgt::ivec2::two,tileSizeRaw);

    const tgt::ivec2 begin = tgt::max(pixBegin, tgt::ivec2::zero);
    const tgt::ivec2 end = tgt::min(pixEnd, texture.dimensions());

    struct CacheData {
        tgt::mat4 voxelToBrick_;
        uint64_t addr_;
        union {
            const uint16_t* data_;
            uint16_t mean_;
        };
    };

    const size_t cacheSize = 8;
    std::array<tgt::ivec3, cacheSize> cacheLlf { tgt::ivec3( 0) };
    std::array<tgt::ivec3, cacheSize> cacheUrb { tgt::ivec3(-1) };
    std::array<CacheData,  cacheSize> cacheData { CacheData { tgt::mat4::identity, 0, nullptr } };
    int nextOut = 0;

    //TODO:
    // linear sampling
    // multithreading?
    // screen reduction?

    for(int ty=begin.y + progress.nextTileY_*tileSize.y; ty<end.y; ty+=tileSize.y) {
        for(int tx=begin.x + progress.nextTileX_*tileSize.x; tx<end.x; tx+=tileSize.x) {
            if(Clock::now() > deadline) {
                return DeadlineResult::TimedOut;
            }

            const tgt::ivec2 tileBegin(tx,ty);
            const tgt::ivec2 tileEnd = tgt::min(tileBegin + tileSize, end);

            if(numChannels != 1) {
                for(int py=tileBegin.y; py<tileEnd.y; ++py) {
                    for(int px=tileBegin.x; px<tileEnd.x; ++px) {
                        int index = px + py*rowSize;
                        pixels[index] = OctreeSliceTexture::Pixel(0);
                    }
                }
            }
            for(size_t channel = 0; channel < numChannels; ++channel) {
                const tgt::mat4 pixelToVoxel = pixelToVoxelMats[channel];
                const TransferFuncCPU& tf = tfs[channel];
                for(int py=tileBegin.y; py<tileEnd.y; ++py) {
                    for(int px=tileBegin.x; px<tileEnd.x; ++px) {
                        tgt::vec4 pixelPos = tgt::vec4(px, py, 0.0, 1.0);
                        tgt::vec3 pos = (pixelToVoxel*pixelPos).xyz();

                        tgt::ivec3 posi = tgt::round(pos);

                        int index = px + py*rowSize;
                        OctreeSliceTexture::Pixel& p = pixels[index];

                        if(tgt::hor(tgt::lessThan(posi, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(posi, volDim))) {
                            // This can happen due to channel shift

                            if(numChannels == 1) {
                                p = OctreeSliceTexture::Pixel(0);
                            }
                            continue;
                        }

                        int cacheEntry = -1;
                        for(int i = 0; i<cacheSize; ++i) {
                            if(tgt::hand(tgt::greaterThanEqual(posi, cacheLlf[i])) && tgt::hand(tgt::lessThan(posi, cacheUrb[i]))) {
                                cacheEntry = i;
                            }
                        }

                        if(cacheEntry == -1) {
                            cacheEntry = nextOut;
                            nextOut = (nextOut+1)%cacheSize;
                            auto& d = cacheData[cacheEntry];
                            if(d.addr_ != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
                                brickPoolManager.releaseBrick(d.addr_, OctreeBrickPoolManagerBase::READ);
                            }

                            LocatedVolumeOctreeNodeConst node = root.findChildNode(posi, brickDataSize, level);
                            auto& location = node.location();
                            d.voxelToBrick_ = location.voxelToBrick();
                            d.addr_ = node.node().getBrickAddress();
                            if(d.addr_ == OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
                                d.mean_ = node.node().getAvgValues()[channel];
                            } else {
                                d.data_ = brickPoolManager.getBrick(d.addr_);
                            }
                            cacheLlf[cacheEntry] = location.voxelLLF();
                            cacheUrb[cacheEntry] = location.voxelURB();
                        }
                        auto& d = cacheData[cacheEntry];

                        tgt::svec3 brickPos = tgt::round(d.voxelToBrick_ * pos);

                        uint16_t rawVal;
                        if(d.addr_ != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
                            size_t brickIndex = channel+numChannels*cubicCoordToLinear(brickPos, brickDataSize);

                            rawVal = d.data_[brickIndex];
                        } else {
                            rawVal = d.mean_;
                        }
                        const float factor = 1.0f/0xffff;
                        float voxelValue = static_cast<float>(rawVal) * factor;

                        OctreeSliceTexture::Pixel np = tf.lookUp(voxelValue);

                        // Compositing:
                        // Looking at sl_base frag: For some reason the compositing
                        // premultiplies alpha for multi channel, but not for
                        // single? I guess we mirror that here...
                        if(numChannels == 1) {
                            p = np;
                        } else {
                            tgt::Vector3<uint8_t> multiplied = (tgt::Vector3<uint16_t>(np.xyz()) * static_cast<uint16_t>(np.a)) / static_cast<uint16_t>(256);
                            p.xyz() += multiplied;
                            p.a = std::max(np.a, p.a);
                        }
                    }
                }
            }
            progress.nextTileX_ += 1;
        }
        progress.nextTileX_ = 0;
        progress.nextTileY_ += 1;
    }
    for(auto d : cacheData) {
        if(d.addr_ != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
            brickPoolManager.releaseBrick(d.addr_, OctreeBrickPoolManagerBase::READ);
        }
    }

    return DeadlineResult::Succeeded;
}

void SliceViewer::renderFromOctree() {
    octreeTexture_->updateDimensions(outport_.getSize());


    // First update texture buffer (on cpu)
    const VolumeBase* volume = inport_.getData();
    tgtAssert(volume, "No volume");
    tgtAssert(volume->hasRepresentation<VolumeOctree>(), "No octree");
    const VolumeOctree& octree = *volume->getRepresentation<VolumeOctree>();
    tgt::ivec3 volDim = volume->getDimensions();

    const SliceAlignment alignment = sliceAlignment_.getValue();
    int numSlices = volDim[alignment];
    const size_t sliceIndex = static_cast<size_t>(sliceIndex_.get());
    int numSlicesCol = numGridCols_.get();
    int numSlicesRow = numGridRows_.get();

    TransferFuncCPU tfs[4] = {
        {transferFunc1_.get(), volume->getRealWorldMapping()},
        {transferFunc2_.get(), volume->getRealWorldMapping()},
        {transferFunc3_.get(), volume->getRealWorldMapping()},
        {transferFunc4_.get(), volume->getRealWorldMapping()},
    };

    TimePoint deadline = Clock::now();
    size_t renderTimeMs = sliceExtractionTimeLimit_.get();
    if(QualityMode.getQuality() == VoreenQualityMode::RQ_HIGH || renderTimeMs == 0) {
        deadline += std::chrono::hours(999);
    } else {
        deadline += std::chrono::milliseconds(renderTimeMs);
    }

    tgt::ivec2 ll = tgt::round(sliceLowerLeft_);
    tgt::ivec2 sliceSize = tgt::round(sliceSize_); //TODO check rounding here is correct

    sliceComplete_ = true; // Will be reset if we encounter a timeout

    for (int pos = octreeRenderProgress_.nextSlice_; pos < (numSlicesCol * numSlicesRow); ++pos) {
        int x = pos % numSlicesCol;
        int y = pos / numSlicesCol;

        int sliceNumber = (pos + static_cast<const int>(sliceIndex));
        if (sliceNumber >= numSlices)
            break;

        tgt::ivec2 begin(ll.x + x*sliceSize.x, ll.y + ((numSlicesRow - (y + 1)) * sliceSize.y));
        tgt::ivec2 end = begin+sliceSize;
        tgt::vec3 dimInv = tgt::vec3(1.0)/tgt::max(tgt::vec3(1.0f),tgt::vec3(end-begin, numSlices));
        tgt::vec3 offset = tgt::vec3(-begin, sliceNumber) + tgt::vec3(0.5);
        tgt::mat4 pixelToScreenNorm = tgt::mat4::createScale(dimInv) * tgt::mat4::createTranslation(offset);
        tgt::mat4 screenNormToTexture = textureMatrix_;
        tgt::mat4 textureToVoxel = volume->getTextureToVoxelMatrix();

        tgt::mat4 pixelToVoxel = textureToVoxel * screenNormToTexture * pixelToScreenNorm;

        tgt::mat4 pixelToVoxelMats[4] {
            applyChannelShift_.get() ? tgt::mat4::createTranslation(channelShift1_.get()) * pixelToVoxel : pixelToVoxel,
            applyChannelShift_.get() ? tgt::mat4::createTranslation(channelShift2_.get()) * pixelToVoxel : pixelToVoxel,
            applyChannelShift_.get() ? tgt::mat4::createTranslation(channelShift3_.get()) * pixelToVoxel : pixelToVoxel,
            applyChannelShift_.get() ? tgt::mat4::createTranslation(channelShift4_.get()) * pixelToVoxel : pixelToVoxel,
        };

        float pixelDistX = tgt::distance(pixelToVoxel*tgt::vec3::zero, pixelToVoxel*tgt::vec3(1,0,0));
        float pixelDistY = tgt::distance(pixelToVoxel*tgt::vec3::zero, pixelToVoxel*tgt::vec3(0,1,0));
        int baseLevel = std::floor(std::log2(std::min(pixelDistX, pixelDistY)));

        int rawLevel = baseLevel;
        switch(QualityMode.getQuality()) {
            case VoreenQualityMode::RQ_INTERACTIVE:
                rawLevel += interactionLevelOfDetail_.get();
                break;
            case VoreenQualityMode::RQ_DEFAULT:
                rawLevel += sliceLevelOfDetail_.get();
                break;
            case VoreenQualityMode::RQ_HIGH:
                //no time limit, octree level 0
                rawLevel = 0;
                break;
            default:
                tgtAssert(false,"unknown rendering quality");
        }

        size_t level = tgt::clamp(rawLevel, 0, static_cast<int>(octree.getNumLevels()-1));

        if(renderOctreeSlice(*octreeTexture_, octree, begin, end, pixelToVoxelMats, level, tfs, octreeRenderProgress_, deadline) == DeadlineResult::Succeeded) {
            octreeRenderProgress_.nextSlice_ += 1;
            octreeRenderProgress_.nextTileY_ = 0;
            octreeRenderProgress_.nextTileX_ = 0;
        } else {
            sliceComplete_ = false;
            break;
        }
    }


    // Update and render texture
    octreeTexture_->uploadTexture();

    tgt::TextureUnit unit;
    unit.activate();
    octreeTexture_->bindTexture();

    copyImageShader_->activate();
    copyImageShader_->setUniform("colorTex_", unit.getUnitNumber());
    renderQuad();
    copyImageShader_->deactivate();
    LGL_ERROR;
}

void SliceViewer::renderFromVolumeTexture(tgt::mat4 toSliceCoordMatrix) {
    const VolumeBase* volume = inport_.getData();
    tgt::ivec3 volDim = volume->getDimensions();
    const SliceAlignment alignment = sliceAlignment_.getValue();
    int numSlices = volDim[alignment];
    float canvasWidth = static_cast<float>(outport_.getSize().x);
    float canvasHeight = static_cast<float>(outport_.getSize().y);
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
    int numSlicesCol = numGridCols_.get();
    int numSlicesRow = numGridRows_.get();

    // render slice(s) as textured quads
    float depth = 0.0f;
    const size_t sliceIndex = static_cast<size_t>(sliceIndex_.get());
    for (int pos = 0, x = 0, y = 0; pos < (numSlicesCol * numSlicesRow);
        ++pos, x = pos % numSlicesCol, y = pos / numSlicesCol)
    {
        int sliceNumber = (pos + static_cast<const int>(sliceIndex));
        if (sliceNumber >= numSlices)
            break;

        // compute depth in texture coordinates and check if it is not below the first or above the last slice
        depth = (static_cast<float>(sliceNumber) + 0.5f) / static_cast<float>(volume->getDimensions()[alignment]);
        float minDepth = 0.5f / static_cast<float>(volume->getDimensions()[alignment]);
        float maxDepth = (static_cast<float>(volume->getDimensions()[alignment]) - 0.5f) / static_cast<float>(volume->getDimensions()[alignment]);
        if ((depth < minDepth) || depth > maxDepth)
            continue;

        MatStack.loadIdentity();
        MatStack.translate(sliceLowerLeft_.x + (x * sliceSize_.x),
            sliceLowerLeft_.y + ((numSlicesRow - (y + 1)) * sliceSize_.y), 0.0f);
        MatStack.scale(sliceSize_.x, sliceSize_.y, 1.0f);

        if(!sliceShader_->isActivated())
            sliceShader_->activate();

        setGlobalShaderParameters(sliceShader_);

        // extract 2D slice
        if (texMode_.isSelected("2d-texture")) {

            const size_t sliceID = tgt::iround(depth * volume->getDimensions()[alignment] - 0.5f);

            bool singleSliceComplete = true;
            SliceTexture* slice = 0;
            int* shiftArray = 0;
            if(applyChannelShift_.get()) {
                //create shift array
                shiftArray = new int[volume->getNumChannels()];
                switch(volume->getNumChannels()) {
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
                    tgtAssert(false,"unsupported channel count!");
                }
            }
            switch(QualityMode.getQuality()) {
            case VoreenQualityMode::RQ_INTERACTIVE:
                slice = sliceCache_.getVolumeSlice(volume, alignment, sliceID, shiftArray, interactionLevelOfDetail_.get(),
                    static_cast<clock_t>(sliceExtractionTimeLimit_.get()), &singleSliceComplete, false);
            break;
            case VoreenQualityMode::RQ_DEFAULT:
                slice = sliceCache_.getVolumeSlice(volume, alignment, sliceID, shiftArray, sliceLevelOfDetail_.get(),
                    static_cast<clock_t>(sliceExtractionTimeLimit_.get()), &singleSliceComplete, false);
            break;
            case VoreenQualityMode::RQ_HIGH:
                //no time limit, octree level 0
                slice = sliceCache_.getVolumeSlice(volume, alignment, sliceID, shiftArray, 0,
                    static_cast<clock_t>(0), &singleSliceComplete, false);
            break;
            default:
                tgtAssert(false,"unknown rendering quality");
            }

            delete[] shiftArray;

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
        LGL_ERROR;
        tgt::vec2 texLowerLeft = tgt::vec2(0.f);
        tgt::vec2 texUpperRight = tgt::vec2(1.f);
        renderSliceGeometry(tgt::vec4(texLowerLeft.x, texLowerLeft.y, depth, 1.f),
            tgt::vec4(texUpperRight.x, texLowerLeft.y, depth, 1.f),
            tgt::vec4(texUpperRight.x, texUpperRight.y, depth, 1.f),
            tgt::vec4(texLowerLeft.x, texUpperRight.y, depth, 1.f));

    }   // textured slices

    tgtAssert(sliceShader_, "no slice shader");
    sliceShader_->deactivate();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;
}

void SliceViewer::process() {
    const VolumeBase* volume = inport_.getData();

    // make sure VolumeGL is available in 3D texture mode
    if (texMode_.isSelected("3d-texture")) {
        if (!volume->getRepresentation<VolumeGL>()) {
            LERROR("3D texture could not be created. Falling back to 2D texture mode.");
            texMode_.select("2d-texture");
            rebuildShader();
            return;
        }
    }
    // ... and octree is available for octree mode
    if (texMode_.isSelected("octree")) {
        if (!volume->hasRepresentation<VolumeOctree>()) {
            LERROR("Octree representation is not available. Falling back to 2D texture mode.");
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
    int numSlices = volDim[alignment];

    // get the textures dimensions
    tgt::vec3 urf = inport_.getData()->getURB();
    tgt::vec3 llb = inport_.getData()->getLLF();

    // Re-calculate texture dimensions, urf, llb and center of the texture
    // for it might be a NPOT texture and therefore might have been inflated.
    // In that case, the inflating need to be undone and the volume texture
    // needs to be "cropped" to its original measures.
    const tgt::vec3 texDim = urf - llb;
    urf = texDim / 2.f;
    llb = texDim / -2.f;

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

    if (texMode_.isSelected("3d-texture") || texMode_.isSelected("2d-texture")) {
        renderFromVolumeTexture(toSliceCoordMatrix);
    } else if (texMode_.isSelected("octree")) {
        renderFromOctree();
    } else {
        LWARNING("Unknown texture mode: " << texMode_.get());
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

void SliceViewer::renderSliceGeometry(const tgt::vec4& t0, const tgt::vec4& t1, const tgt::vec4& t2, const tgt::vec4& t3) const {
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
    glVertexAttribPointer(0, 2, GL_FLOAT, false, sizeof(VertexHelper), 0);
    glVertexAttribPointer(3, 4, GL_FLOAT, false, sizeof(VertexHelper), (void*)sizeof(tgt::vec2));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(3);

    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

    glDeleteBuffers(1, &vbo);
    glBindVertexArray(0);
    glDeleteVertexArrays(1, &vao);
}

std::string SliceViewer::generateSliceShaderHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = VolumeRenderer::generateHeader();

    if (texMode_.isSelected("2d-texture"))
        header += "#define SLICE_TEXTURE_MODE_2D \n";
    else if (texMode_.isSelected("3d-texture") || texMode_.isSelected("octree"))
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

bool SliceViewer::rebuildShader() {
    // do nothing if there is no shader at the moment
    if (!isInitialized() || !sliceShader_)
        return false;

    if(texMode_.isSelected("3d-texture") || texMode_.isSelected("2d-texture")) {
        sliceShader_->setHeaders(generateSliceShaderHeader());
        return sliceShader_->rebuild();
    } else {
        return true;
    }
}

void SliceViewer::renderLegend() {

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

void SliceViewer::renderSliceBoundaries() const {

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

void SliceViewer::renderInfoTexts() const {

    if (showCursorInfos_.isSelected("never") && !showSliceNumber_.get())
        return;

    int numSlicesRow = numGridRows_.get();
    int numSlicesCol = numGridCols_.get();
    tgtAssert(numSlicesRow > 0 && numSlicesCol > 0, "Invalid slice counts");

    tgtAssert(inport_.getData(), "No volume");
    const VolumeBase* volume = inport_.getData();
    tgt::ivec3 volDim = volume->getDimensions();
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
                            break;
                        case 1:
                            tmp = lastPickingPosition_ + tgt::ivec3(channelShift2_.get());
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
                            break;
                        case 2:
                            tmp = lastPickingPosition_ + tgt::ivec3(channelShift3_.get());
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
                            break;
                        case 3:
                            tmp = lastPickingPosition_ + tgt::ivec3(channelShift4_.get());
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
                            break;
                        default:
                            tgtAssert(false,"Should not get here");
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

    if (showSliceNumber_.get()) {
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

tgt::vec3 SliceViewer::screenToVoxelPos(tgt::ivec2 screenPos) const {

    if (!inport_.getData() || !outport_.getRenderTarget())
        return tgt::vec3(-1.f);

    tgt::vec3 volumeDim(inport_.getData()->getDimensions());
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
        return tgt::vec3(-1.f);

    // determine the picked slice
    const int sliceColID = p.x / sliceSizeInt.x;
    const int sliceRowID = (numSlicesRow-1) - (p.y / sliceSizeInt.y);
    const int slice = sliceColID + (sliceRowID * numSlicesCol) + sliceIndex_.get();

    // calculate the normalized position within the picked slice
    tgt::vec2 posWithinSlice(
        static_cast<float>(p.x % sliceSizeInt.x),
        static_cast<float>(p.y % sliceSizeInt.y));
    posWithinSlice /= sliceSize_;

    // calculate the normalized depth of the picked slice (texture z coordinate)
    float depth = (static_cast<float>(slice) + 0.5f) / std::max(volumeDim[voxelPosPermutation_.z], 1.f);

    // now we have the assigned texture coordinates of the picked fragment
    tgt::vec4 texCoords(posWithinSlice, depth, 1.f);
    texCoords = tgt::clamp(texCoords, tgt::vec4(0.f), tgt::vec4(1.f));

    // apply current texture matrix to assigned tex coords
    tgt::vec3 texCoordsTransformed = (textureMatrix_ * texCoords).xyz();

    // transform final tex coords into volume coordinates
    tgt::vec3 voxPos = texCoordsTransformed * (volumeDim) - tgt::vec3(0.5f);
    voxPos = tgt::clamp(voxPos, tgt::vec3(0.f), tgt::vec3(volumeDim-1.f));

    return voxPos;
}

tgt::mat4 SliceViewer::generatePickingMatrix() const {

    if (!inport_.hasData())
        return tgt::mat4::createIdentity();

    tgt::vec3 volumeDim(inport_.getData()->getDimensions());

    // 1. translate slice to origin (also add 0.5 in depth direction, since we subtract this afterwards and the slice number should already be correct)
    tgt::mat4 originTranslation = tgt::mat4::createTranslation(tgt::vec3(-sliceLowerLeft_.x, -sliceLowerLeft_.y, 0.5f));

    // 2. normalize screen coords with regard to the slice
    tgt::mat4 sliceScale = tgt::mat4::createScale(tgt::vec3(1.f / sliceSize_.x, 1.f / sliceSize_.y, 1.f / (volumeDim[voxelPosPermutation_.z])));

    // 3. apply current texture matrix
    //tgt::mat4 textureMatrix = textureMatrix_; // used directly (see below)

    // 4. transform normalized coordinates to volume dimensions
    tgt::mat4 volumeScale = tgt::mat4::createTranslation(tgt::vec3(-0.5f)) * tgt::mat4::createScale(volumeDim);

    // compose transformation matrix
    tgt::mat4 result = volumeScale * textureMatrix_ * sliceScale * originTranslation;

    return result;
}

void SliceViewer::onSliceAlignmentChange() {
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

void SliceViewer::mouseLocalization(tgt::MouseEvent* e) {

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

void SliceViewer::shiftEvent(tgt::MouseEvent* e) {

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

bool SliceViewer::singleSliceMode() const {
    return (numGridRows_.get() == 1 && numGridCols_.get() == 1);
}

void SliceViewer::resetView()
{
    zoomFactor_.reset();
    voxelOffset_.reset();
}

} // namespace voreen
