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

#include "floworthogonalslicerenderer.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/numericproperty.h"

#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"

using tgt::TextureUnit;

namespace voreen {

FlowOrthogonalSliceRenderer::FlowOrthogonalSliceRenderer()
    : VolumeRenderer(),
    outport_(Port::OUTPORT, "image.outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
    inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input"),
    transferFunc_("transferFunction", "Transfer Function"),
    texMode_("textureMode", "Texture Mode", Processor::INVALID_PROGRAM),
    sliceShader_(0),
    volumeDimensions_(1),
    slicePositions_(-1, -1, -1),
    useXYSliceProp_("useXYSliceProp", "render xy-slice:", true),
    useXZSliceProp_("useZXSliceProp", "render zx-slice:", false),
    useZYSliceProp_("useZYSliceProp", "render zy-slice:", false),
    sliceNoXYProp_("sliceNoXYProp", "xy-slice number:", 1, 1, 100),
    sliceNoXZProp_("sliceNoZXProp", "zx-slice number:", 1, 1, 100),
    sliceNoZYProp_("sliceNoZYProp", "zy-slice number:", 1, 1, 100),
    camProp_("camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f))),
    cameraHandler_(0),
    cpInport_(Port::INPORT, "coprocessor.slicepositionsInput", "coprocessor.slicepositionsInput", false),
    cpOutport_(Port::OUTPORT, "coprocessor.slicepositions")
{
    //inport_.addCondition(new PortConditionVolumeTypeGL());
    inport_.showTextureAccessProperties(true);
    addPort(inport_);
    addPort(outport_);

    addProperty(transferFunc_);

    // texture mode (2D/3D)
    texMode_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode_.selectByKey("3d-texture");
    addProperty(texMode_);
    texMode_.setGroupID(inport_.getID() + ".textureAccess");

    MemberFunctionCallback<FlowOrthogonalSliceRenderer> sliceChangeAction(this,
        &FlowOrthogonalSliceRenderer::onSlicePositionChange);
    useXYSliceProp_.onChange(sliceChangeAction);
    useXZSliceProp_.onChange(sliceChangeAction);
    useZYSliceProp_.onChange(sliceChangeAction);
    sliceNoXYProp_.onChange(sliceChangeAction);
    sliceNoXZProp_.onChange(sliceChangeAction);
    sliceNoZYProp_.onChange(sliceChangeAction);

    addProperty(useXYSliceProp_);
    addProperty(sliceNoXYProp_);
    addProperty(useXZSliceProp_);
    addProperty(sliceNoXZProp_);
    addProperty(useZYSliceProp_);
    addProperty(sliceNoZYProp_);
    addProperty(camProp_);
    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &camProp_);
    addInteractionHandler(cameraHandler_);

    addPort(cpInport_);
    addPort(cpOutport_);

    onSlicePositionChange();
}

FlowOrthogonalSliceRenderer::~FlowOrthogonalSliceRenderer() {
    delete cameraHandler_;
}

Processor* FlowOrthogonalSliceRenderer::create() const {
    return new FlowOrthogonalSliceRenderer();
}

bool FlowOrthogonalSliceRenderer::isReady() const {
    const std::vector<Port*>& inports = getInports();
    for (size_t i = 0; i < inports.size(); ++i)
        if(!inports[i]->isReady())
            return false;

    const std::vector<Port*>& outports = getOutports();
    for (size_t i = 0; i < outports.size(); ++i)
        if(!outports[i]->isReady())
            return false;

    // The co-processor is optional, so the processor is ready, even if the
    // processor's co-processor inport is not connected()!

    return true;
}

void FlowOrthogonalSliceRenderer::process() {

    if (!inport_.isReady())
        return;

    if (inport_.hasChanged()) {
        volumeDimensions_ = inport_.getData()->getRepresentation<VolumeRAM>()->getDimensions();
        updateNumSlices();  // validate the currently set values and adjust them if necessary
        transferFunc_.setVolume(inport_.getData());
        rebuildShader(); //set num channels
    }

    outport_.activateTarget("FlowOrthogonalSliceRenderer::process()");

    const VolumeGL* volumeGL = inport_.getData()->getRepresentation<VolumeGL>();
    const VolumeTexture* const tex = volumeGL->getTexture();
    if (!tex) {
        LERROR("setVolumeHandle(): VolumeTexture in VolumGL is NULL!");
        return;
    }

    // activate the shader and set the needed uniforms
    TextureUnit volumeUnit, transferUnit;
    tgt::Camera cam = camProp_.get();
    sliceShader_->activate();
    LGL_ERROR;
    if(inport_.getData()->getNumChannels() >= 1)
        transferFunc_.get()->setUniform(sliceShader_, "transFuncParams_", "transFuncTex_", transferUnit.getUnitNumber());
    if(inport_.getData()->getNumChannels() >= 2)
        transferFunc_.get()->setUniform(sliceShader_, "transFuncParams2_", "transFuncTex2_", transferUnit.getUnitNumber());
    if(inport_.getData()->getNumChannels() >= 3)
        transferFunc_.get()->setUniform(sliceShader_, "transFuncParams3_", "transFuncTex3_", transferUnit.getUnitNumber());
    if(inport_.getData()->getNumChannels() >= 4)
        transferFunc_.get()->setUniform(sliceShader_, "transFuncParams4_", "transFuncTex4_", transferUnit.getUnitNumber());

    //sliceShader_->setUniform("numChannels_", static_cast<GLint>(inport_.getData()->getNumChannels())); done in header
    LGL_ERROR;

    // bind volume
    std::vector<VolumeStruct> volumeTextures;
    volumeTextures.push_back(VolumeStruct(
        inport_.getData(),
        &volumeUnit,
        "volume_","volumeParams_",
        inport_.getTextureClampModeProperty().getValue(),
        tgt::vec4(inport_.getTextureBorderIntensityProperty().get()),
        inport_.getTextureFilterModeProperty().getValue())
        );
    bindVolumes(sliceShader_, volumeTextures, &cam, lightPosition_.get());
    LGL_ERROR;

    FlowOrthogonalSliceRenderer* other = cpInport_.getConnectedProcessor();
    if (other != 0) {
        tgt::ivec3 positions = tgt::clamp(other->getSlicePositions(), tgt::ivec3(-1),
            (volumeDimensions_ - tgt::ivec3(1)));
        setupSliceProperties(positions);
    } else {
        onSlicePositionChange();
        setSlicePropertiesVisible(true);
    }

    transferUnit.activate();
    transferFunc_.get()->getTexture()->bind();

    sliceShader_->setUniform("textureMatrix_", tgt::mat4::identity);

    // important: save current camera state before using the processor's camera or
    // successive processors will use those settings!
    //
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();

    camProp_.look(outport_.getSize());

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if (useXYSliceProp_.get() == true) {
        float sliceNo = static_cast<float>(sliceNoXYProp_.get() - 1);
        renderSlice(SLICE_XY, sliceNo);
    }

    if (useXZSliceProp_.get() == true) {
        float sliceNo = static_cast<float>(sliceNoXZProp_.get() - 1);
        renderSlice(SLICE_XZ, sliceNo);
    }

    if (useZYSliceProp_.get() == true) {
        float sliceNo = static_cast<float>(sliceNoZYProp_.get() - 1);
        renderSlice(SLICE_ZY, sliceNo);
    }

    sliceShader_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
}

const tgt::ivec3& FlowOrthogonalSliceRenderer::getSlicePositions() const {
    return slicePositions_;
}

// protected methods
//

void FlowOrthogonalSliceRenderer::onSlicePositionChange() {
    slicePositions_.x = sliceNoZYProp_.get();
    slicePositions_.y = sliceNoXZProp_.get();
    slicePositions_.z = sliceNoXYProp_.get();

    if (useXYSliceProp_.get() == false)
        slicePositions_.z *= -1;
    if (useXZSliceProp_.get() == false)
        slicePositions_.y *= -1;
    if (useZYSliceProp_.get() == false)
        slicePositions_.x *= -1;
}

void FlowOrthogonalSliceRenderer::renderSlice(const SliceAlignment& sliceAlign, const float sliceNo)
{
    tgt::vec3 dim = static_cast<tgt::vec3>(volumeDimensions_ - tgt::ivec3(1));
    tgt::vec3 llb = (dim / -tgt::max(dim));
    tgt::vec3 textureSize = (llb * -2.0f);

    tgt::ivec3 permutation(0, 1, 2);
    switch (sliceAlign) {
        default:
        case SLICE_XY:
            break;
        case SLICE_XZ:
            permutation = tgt::ivec3(0, 2, 1);
            break;
        case SLICE_ZY:
            permutation = tgt::ivec3(2, 1, 0);
            break;
    }
    float s = sliceNo / (volumeDimensions_[permutation.z] - 1);

    tgt::vec3 ll = permuteComponents(tgt::vec3(0.0f, 0.0f, s), permutation);
    tgt::vec3 lr = permuteComponents(tgt::vec3(1.0f, 0.0f, s), permutation);
    tgt::vec3 ur = permuteComponents(tgt::vec3(1.0f, 1.0f, s), permutation);
    tgt::vec3 ul = permuteComponents(tgt::vec3(0.0f, 1.0f, s), permutation);

    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.texcoord(ll);
        IMode.vertex( ((ll * textureSize) + llb) );

        IMode.texcoord(lr);
        IMode.vertex( ((lr * textureSize) + llb) );

        IMode.texcoord(ur);
        IMode.vertex( ((ur * textureSize) + llb) );

        IMode.texcoord(ul);
        IMode.vertex( ((ul * textureSize) + llb) );
    IMode.end();
}

void FlowOrthogonalSliceRenderer::initialize() {
    VolumeRenderer::initialize();

    sliceShader_ = ShdrMgr.load("sl_base", buildShaderHeader(), false);
    LGL_ERROR;
}

void FlowOrthogonalSliceRenderer::deinitialize() {
    ShdrMgr.dispose(sliceShader_);
    sliceShader_ = 0;
    VolumeRenderer::deinitialize();
}

// protected methods
//

void FlowOrthogonalSliceRenderer::beforeProcess() {
    VolumeRenderer::beforeProcess();

    if (invalidationLevel_ >= Processor::INVALID_PROGRAM)
        rebuildShader();
}

std::string FlowOrthogonalSliceRenderer::buildShaderHeader() {
    std::string header = VolumeRenderer::generateHeader();

    if (texMode_.isSelected("2d-texture"))
        header += "#define SLICE_TEXTURE_MODE_2D \n";
    else if (texMode_.isSelected("3d-texture"))
        header += "#define SLICE_TEXTURE_MODE_3D \n";
    else {
        LWARNING("Unknown texture mode: " << texMode_.get());
    }

    if(inport_.getData())
        header += "#define NUM_CHANNELS " + itos(inport_.getData()->getNumChannels()) + " \n";

    header += transferFunc_.get()->getShaderDefines();

    return header;
}

bool FlowOrthogonalSliceRenderer::rebuildShader() {
    // do nothing if there is no shader at the moment
    if (!sliceShader_)
        return false;

    sliceShader_->setHeaders(buildShaderHeader());
    return sliceShader_->rebuild();
}

void FlowOrthogonalSliceRenderer::setSlicePropertiesVisible(const bool visible) {
    sliceNoXYProp_.setVisibleFlag(visible);
    useXYSliceProp_.setVisibleFlag(visible);
    sliceNoXZProp_.setVisibleFlag(visible);
    useXZSliceProp_.setVisibleFlag(visible);
    sliceNoZYProp_.setVisibleFlag(visible);
    useZYSliceProp_.setVisibleFlag(visible);
}

void FlowOrthogonalSliceRenderer::setupSliceProperties(const tgt::ivec3& positions) {
    if (slicePositions_ == positions)
        return;

    slicePositions_ = positions;
    bool useForeignSlices = false;
    if (slicePositions_.x > 0) {
        sliceNoZYProp_.set(slicePositions_.x);
        useZYSliceProp_.set(true);
        useForeignSlices = true;
    } else
        useZYSliceProp_.set(false);

    if (slicePositions_.y > 0) {
        sliceNoXZProp_.set(slicePositions_.y);
        useXZSliceProp_.set(true);
        useForeignSlices = true;
    } else
        useXZSliceProp_.set(false);

    if (slicePositions_.z > 0) {
        sliceNoXYProp_.set(slicePositions_.z);
        useXYSliceProp_.set(true);
        useForeignSlices = true;
    } else
        useXYSliceProp_.set(false);

    if (useForeignSlices == true)
        setSlicePropertiesVisible(false);
    else
        setSlicePropertiesVisible(true);

    updateNumSlices();
}

void FlowOrthogonalSliceRenderer::updateNumSlices() {
    // set number of slice for xy-plane (along z-axis)
    //
    int numSlicesXY = volumeDimensions_.z;
    sliceNoXYProp_.setMaxValue(numSlicesXY);
    if (sliceNoXYProp_.get() > numSlicesXY)
        sliceNoXYProp_.set(numSlicesXY / 2);
    sliceNoXYProp_.updateWidgets();

    // set number of slices for zx-plane (along y-axis)
    //
    int numSlicesZX = volumeDimensions_.y;
    sliceNoXZProp_.setMaxValue(numSlicesZX);
    if (sliceNoXZProp_.get() > numSlicesZX)
        sliceNoXZProp_.set(numSlicesZX / 2);
    sliceNoXZProp_.updateWidgets();

    // set number of slices for sagittal plane (along x-axis)
    //
    int numSlicesZY = volumeDimensions_.x;
    sliceNoZYProp_.setMaxValue(numSlicesZY);
    if (sliceNoZYProp_.get() > numSlicesZY)
        sliceNoZYProp_.set(numSlicesZY / 2);
    sliceNoZYProp_.updateWidgets();
}

inline tgt::vec3 FlowOrthogonalSliceRenderer::permuteComponents(const tgt::vec3& input,
                                                            const tgt::ivec3& permutation)
{
    return tgt::vec3(input[permutation.x], input[permutation.y], input[permutation.z]);
}

}   // namespace
