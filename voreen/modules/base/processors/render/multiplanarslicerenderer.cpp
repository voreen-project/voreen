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

#include "multiplanarslicerenderer.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/numericproperty.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "tgt/textureunit.h"
#include "tgt//immediatemode/immediatemode.h"

using tgt::TextureUnit;

inline tgt::vec3 permuteComponents(const tgt::vec3& input,
                                   const tgt::ivec3& permutation) {
    return tgt::vec3(input[permutation.x], input[permutation.y], input[permutation.z]);
}

namespace voreen {

const std::string MultiplanarSliceRenderer::loggerCat_("voreen.base.MultiplanarSliceRenderer");

MultiplanarSliceRenderer::MultiplanarSliceRenderer()
    : VolumeRenderer()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , transferFunc1_("transferFunction", "Transfer Function")
    , transferFunc2_("transferFunction2", "Transfer Function 2")
    , transferFunc3_("transferFunction3", "Transfer Function 3")
    , transferFunc4_("transferFunction4", "Transfer Function 4")
    , texMode_("textureMode", "Texture Mode", Processor::INVALID_PROGRAM)
    , renderXYSlice_("renderSlice.XY", "Render XY Slice", true)
    , renderXZSlice_("renderSlice.XZ", "Render XZ Slice", true)
    , renderYZSlice_("renderSlice.YZ", "Render YZ Slice", true)
    , sliceNumberXY_("sliceNumber.XY", "XY Slice Number", 0, 0, 10000)
    , sliceNumberXZ_("sliceNumber.XZ", "XZ Slice Number", 0, 0, 10000)
    , sliceNumberYZ_("sliceNumber.YZ", "YZ Slice Number", 0, 0, 10000)
    , camProp_("camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f)), true)
    , cameraHandler_(0)
    , sliceShader_(0)
{
    // texture mode (2D/3D)
    texMode_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode_.selectByKey("3d-texture");
    ON_CHANGE(texMode_, MultiplanarSliceRenderer, rebuildShader);
    addProperty(texMode_);
    texMode_.setGroupID(inport_.getID() /*+ ".textureAccess"*/);

    //inport_.addCondition(new PortConditionVolumeTypeGL());
    //inport_.showTextureAccessProperties(true);
    addPort(inport_);
    addPort(outport_);

    addProperty(transferFunc1_);
    addProperty(transferFunc2_);
    transferFunc2_.setVisibleFlag(false);
    addProperty(transferFunc3_);
    transferFunc3_.setVisibleFlag(false);
    addProperty(transferFunc4_);
    transferFunc4_.setVisibleFlag(false);

    addProperty(renderXYSlice_);
    addProperty(sliceNumberXY_);
    addProperty(renderXZSlice_);
    addProperty(sliceNumberXZ_);
    addProperty(renderYZSlice_);
    addProperty(sliceNumberYZ_);
    addProperty(camProp_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &camProp_);
    addInteractionHandler(cameraHandler_);
}

MultiplanarSliceRenderer::~MultiplanarSliceRenderer() {
    delete cameraHandler_;
}

Processor* MultiplanarSliceRenderer::create() const {
    return new MultiplanarSliceRenderer();
}

void MultiplanarSliceRenderer::initialize() {
    VolumeRenderer::initialize();

    sliceShader_ = ShdrMgr.load("sl_base", generateHeader(), false);
    LGL_ERROR;
}

void MultiplanarSliceRenderer::deinitialize() {
    tgtAssert(sliceShader_, "no shader");
    ShdrMgr.dispose(sliceShader_);
    sliceShader_ = 0;
    VolumeRenderer::deinitialize();
}

void MultiplanarSliceRenderer::process() {
    tgtAssert(sliceShader_, "no shader");
    tgtAssert(inport_.isReady(), "Inport not ready");
    const VolumeBase* inputVolume = inport_.getData();

    // make sure VolumeGL is available in 3D texture mode
    if (texMode_.isSelected("3d-texture")) {
        if (!inputVolume->getRepresentation<VolumeGL>()) {
            LERROR("3D texture could not be created. Falling back to 2D texture mode.");
            texMode_.select("2d-texture");
            rebuildShader();
            return;
        }
    }

    outport_.activateTarget("OrthogonalSliceRenderer::process()");
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // get GL resources and setup shader
    sliceShader_->activate();
    LGL_ERROR;

    // bind texture
    bool setupSuccessful = false;
    TextureUnit texUnit;
    if (texMode_.isSelected("2d-texture")) {      // 2D texture
        setupSuccessful = true;
    }
    else if (texMode_.isSelected("3d-texture")) { // 3D texture
        // bind volume
        std::vector<VolumeStruct> volumeTextures;
        volumeTextures.push_back(VolumeStruct(
            inport_.getData(),
            &texUnit,
            "volume_","volumeParams_",
            inport_.getTextureClampModeProperty().getValue(),
            tgt::vec4(inport_.getTextureBorderIntensityProperty().get()),
            inport_.getTextureFilterModeProperty().getValue())
        );
        setupSuccessful = bindVolumes(sliceShader_, volumeTextures, 0, lightPosition_.get());
        LGL_ERROR;

        if(!setupSuccessful)
            LERROR("3D texture could not been created. Try 2D texture mode.");
    }
    else {
        LERROR("unknown texture mode: " << texMode_.get());
        setupSuccessful = false;
    }
    if (!setupSuccessful) {
        outport_.deactivateTarget();
        return;
    }

    // bind transfer functions
    TextureUnit transferUnit1, transferUnit2, transferUnit3, transferUnit4;
    transferUnit1.activate();
    transferFunc1_.get()->getTexture()->bind();
    transferFunc1_.get()->setUniform(sliceShader_, "transFuncParams_", "transFuncTex_", transferUnit1.getUnitNumber());
    LGL_ERROR;
    if (inputVolume->getNumChannels() > 1) {
        transferUnit2.activate();
        transferFunc2_.get()->getTexture()->bind();
        transferFunc2_.get()->setUniform(sliceShader_, "transFuncParams2_", "transFuncTex2_", transferUnit2.getUnitNumber());
        LGL_ERROR;
    }
    if (inputVolume->getNumChannels() > 2) {
        transferUnit3.activate();
        transferFunc3_.get()->getTexture()->bind();
        transferFunc3_.get()->setUniform(sliceShader_, "transFuncParams3_", "transFuncTex3_", transferUnit3.getUnitNumber());
        LGL_ERROR;
    }
    if (inputVolume->getNumChannels() > 3) {
        transferUnit4.activate();
        transferFunc4_.get()->getTexture()->bind();
        transferFunc4_.get()->setUniform(sliceShader_, "transFuncParams4_", "transFuncTex4_", transferUnit4.getUnitNumber());
        LGL_ERROR;
    }

    sliceShader_->setUniform("textureMatrix_", tgt::mat4::identity);

    // important: save current camera state before using the processor's camera or
    // successive processors will use those settings!
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();

    camProp_.look(outport_.getSize());

    // render in texture coordinates and use the volume's textureToWorld matrix
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(inport_.getData()->getTextureToWorldMatrix());

    if (renderXYSlice_.get()) {
        renderSlice(XY_PLANE, sliceNumberXY_.get(), texUnit);
    }
    if (renderXZSlice_.get()) {
        renderSlice(XZ_PLANE, sliceNumberXZ_.get(), texUnit);
    }
    if (renderYZSlice_.get()) {
        renderSlice(YZ_PLANE, sliceNumberYZ_.get(), texUnit);
    }

    // restore matrix stack
    MatStack.popMatrix();

    tgtAssert(sliceShader_, "no shader");
    sliceShader_->deactivate();

    glActiveTexture(GL_TEXTURE0);
    outport_.deactivateTarget();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
}

std::string MultiplanarSliceRenderer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = VolumeRenderer::generateHeader();

    if (texMode_.isSelected("2d-texture"))
        header += "#define SLICE_TEXTURE_MODE_2D \n";
    else if (texMode_.isSelected("3d-texture"))
        header += "#define SLICE_TEXTURE_MODE_3D \n";
    else {
        LWARNING("Unknown texture mode: " << texMode_.get());
    }

    header += "#define NUM_CHANNELS " + (inport_.hasData() ? itos(inport_.getData()->getNumChannels()) : "1") + "\n";

    header += transferFunc1_.get()->getShaderDefines();

    return header;
}

void MultiplanarSliceRenderer::rebuildShader() {
    // do nothing if there is no shader at the moment
    if (!sliceShader_)
        return;

    sliceShader_->setHeaders(generateHeader());
    sliceShader_->rebuild();
    LGL_ERROR;
}

// protected methods
//

void MultiplanarSliceRenderer::renderSlice(SliceAlignment sliceAlign, int sliceNo, tgt::TextureUnit& texUnit) {
    tgtAssert(sliceShader_, "no shader");
    tgtAssert(inport_.hasData(), "No volume");
    const VolumeBase* volume =  inport_.getData();
    tgt::ivec3 volDim = volume->getDimensions();

    tgt::ivec3 permutation(0, 1, 2);
    tgt::mat4 toSliceCoordMatrix = tgt::mat4::zero;
    switch (sliceAlign) {
        default:
    case XY_PLANE:
            toSliceCoordMatrix.t00 = 1.f;
            toSliceCoordMatrix.t11 = 1.f;
            break;
        case XZ_PLANE:
            permutation = tgt::ivec3(0, 2, 1);
            toSliceCoordMatrix.t00 = 1.f;
            toSliceCoordMatrix.t12 = 1.f;
            break;
        case YZ_PLANE:
            permutation = tgt::ivec3(2, 1, 0);
            toSliceCoordMatrix.t01 = 1.f;
            toSliceCoordMatrix.t12 = 1.f;
            break;
    }
    // compute texture coordinates from slice number and clamp to first and last slice
    float s = (static_cast<float>(sliceNo) + 0.5f) / static_cast<float>(volDim[permutation.z]);

    tgt::vec3 ll = permuteComponents(tgt::vec3(0.f, 0.f, s), permutation);
    tgt::vec3 lr = permuteComponents(tgt::vec3(1.f, 0.f, s), permutation);
    tgt::vec3 ur = permuteComponents(tgt::vec3(1.f, 1.f, s), permutation);
    tgt::vec3 ul = permuteComponents(tgt::vec3(0.f, 1.f, s), permutation);

    if (texMode_.isSelected("2d-texture")) {
        // extract 2D slice
        SliceTexture* slice = SliceHelper::getVolumeSlice(volume, sliceAlign, sliceNo);
        if (!slice) {
            LERROR("Failed to extract slice");
            return;
        }

        // bind slice texture
        GLint texFilterMode = inport_.getTextureFilterModeProperty().getValue();
        GLint texClampMode = inport_.getTextureClampModeProperty().getValue();
        tgt::vec4 borderColor = tgt::vec4(inport_.getTextureBorderIntensityProperty().get());
        if (!GLSL::bindSliceTexture(slice, &texUnit, texFilterMode, texClampMode, borderColor)) {
            LERROR("Failed to bind slice texture");
            return;
        }

        // pass slice uniforms to shader
        sliceShader_->setIgnoreUniformLocationError(true);
        GLSL::setUniform(sliceShader_, "sliceTex_", "sliceTexParams_", slice, &texUnit);
        sliceShader_->setUniform("toSliceCoordMatrix", toSliceCoordMatrix);
        sliceShader_->setIgnoreUniformLocationError(false);
        sliceShader_->setUniform("textureMatrix_", tgt::mat4::identity);
    }

    // render slice (vertices are rendered using texture coordinates)
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.texcoord(ll);
        IMode.vertex(ll);

        IMode.texcoord(lr);
        IMode.vertex(lr);

        IMode.texcoord(ur);
        IMode.vertex(ur);

        IMode.texcoord(ul);
        IMode.vertex(ul);
    IMode.end();
}

void MultiplanarSliceRenderer::adjustPropertiesToInput() {
    // get volume
    const VolumeBase* inputVolume = inport_.getData();

    if (!inputVolume)
        return;

    // adapt camera if necessary
    camProp_.adaptInteractionToScene(inputVolume->getBoundingBox().getBoundingBox());

    tgt::ivec3 volDim = static_cast<tgt::ivec3>(inputVolume->getDimensions());
    int numChannels = static_cast<int>(inputVolume->getNumChannels());

    // set number of slice for xy-plane (along z-axis)
    int tmpNo = sliceNumberXY_.get();
    sliceNumberXY_.setMaxValue(volDim.z-1);
    if (tmpNo >= volDim.z-1)
        sliceNumberXY_.set(volDim.z / 2);
    sliceNumberXY_.updateWidgets();

    // set number of slices for zx-plane (along y-axis)
    tmpNo = sliceNumberXZ_.get();
    sliceNumberXZ_.setMaxValue(volDim.y-1);
    if (tmpNo >= volDim.y-1)
        sliceNumberXZ_.set(volDim.y / 2);
    sliceNumberXZ_.updateWidgets();

    // set number of slices for sagittal plane (along x-axis)
    tmpNo = sliceNumberYZ_.get();
    sliceNumberYZ_.setMaxValue(volDim.x-1);
    if (tmpNo >= volDim.x-1)
        sliceNumberYZ_.set(volDim.x / 2);
    sliceNumberYZ_.updateWidgets();

    // transfer functions
    transferFunc2_.setVisibleFlag(numChannels > 1);
    transferFunc3_.setVisibleFlag(numChannels > 2);
    transferFunc4_.setVisibleFlag(numChannels > 3);

    transferFunc1_.setVolume(inputVolume, 0);
    if (inputVolume->getNumChannels() > 1)
        transferFunc2_.setVolume(inputVolume, 1);
    if (inputVolume->getNumChannels() > 2)
        transferFunc3_.setVolume(inputVolume, 2);
    if (inputVolume->getNumChannels() > 3)
        transferFunc4_.setVolume(inputVolume, 3);

    rebuildShader();
}

}   // namespace
