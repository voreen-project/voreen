/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/volume/slice/slicehelper.h"

namespace voreen {

const std::string UnalignedSliceViewer::loggerCat_("voreen.flowreen.UnalignedSliceViewer");

UnalignedSliceViewer::UnalignedSliceViewer()
    : VolumeRenderer()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , normal_("normal", "Normal", tgt::vec3(0, 1, 0), -tgt::vec3::one, tgt::vec3::one)
    , distance_("distance", "Distance", 0, -1000.0f, 1000.0f)
    , samplingRate_("samplingRate", "Sampling Rate", 2.0f, 0.5f, 8.0f)
    , transferFunc_("transferFunction", "Transfer Function")
    , sliceShader_(nullptr)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(normal_);
    addProperty(distance_);
    addProperty(samplingRate_);
    addProperty(transferFunc_);
}

Processor* UnalignedSliceViewer::create() const {
    return new UnalignedSliceViewer();
}

void UnalignedSliceViewer::initialize() {
    VolumeRenderer::initialize();

    sliceShader_ = ShdrMgr.load("sl_base", generateHeader(), false);
    LGL_ERROR;
}

void UnalignedSliceViewer::deinitialize() {
    ShdrMgr.dispose(sliceShader_);
    sliceShader_ = nullptr;

    VolumeRenderer::deinitialize();
}

void UnalignedSliceViewer::beforeProcess() {
    VolumeRenderer::beforeProcess();

    if (inport_.hasChanged()) {
        const VolumeBase* inputVolume = inport_.getData();
        transferFunc_.setVolume(inputVolume, 0);
    }

    if (invalidationLevel_ >= Processor::INVALID_PROGRAM || inport_.hasChanged())
        rebuildShader();
}

void UnalignedSliceViewer::afterProcess() {
    VolumeRenderer::afterProcess();
}

void UnalignedSliceViewer::process() {
    const VolumeBase* volume = inport_.getData();

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;
/*
    float canvasWidth = static_cast<float>(outport_.getSize().x);
    float canvasHeight = static_cast<float>(outport_.getSize().y);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    MatStack.multMatrix(tgt::mat4::createOrtho(0.0f, canvasWidth, 0.f, canvasHeight, -1.0f, 1.0f));
    LGL_ERROR;

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    LGL_ERROR;
*/
    sliceShader_->activate();
    tgt::TextureUnit texUnit;
    setGlobalShaderParameters(sliceShader_);

    tgt::TextureUnit transferUnit;
    transferUnit.activate();
    transferFunc_.get()->getTexture()->bind();
    transferFunc_.get()->setUniform(sliceShader_, "transFuncParams_", "transFuncTex_", transferUnit.getUnitNumber());
    LGL_ERROR;

    tgt::plane plane(normal_.get(), distance_.get());
    texUnit.activate();
    std::unique_ptr<SliceTexture> slice(SliceHelper::getVolumeSlice(volume, plane, samplingRate_.get()));
    LGL_ERROR;

    // bind slice texture
    GLint texFilterMode = inport_.getTextureFilterModeProperty().getValue();
    GLint texClampMode = inport_.getTextureClampModeProperty().getValue();
    tgt::vec4 borderColor = tgt::vec4(inport_.getTextureBorderIntensityProperty().get());
    GLSL::bindSliceTexture(slice.get(), &texUnit, texFilterMode, texClampMode, borderColor);

    // pass slice uniforms to shader
    sliceShader_->setIgnoreUniformLocationError(true);
    GLSL::setUniform(sliceShader_, "sliceTex_", "sliceTexParams_", slice.get(), &texUnit);
    sliceShader_->setIgnoreUniformLocationError(false);

    sliceShader_->setUniform("toSliceCoordMatrix", tgt::mat4::identity);
    sliceShader_->setUniform("textureMatrix_", tgt::mat4::identity);

    renderQuad();

    sliceShader_->deactivate();
    slice->disable();
/*
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;
*/
    outport_.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void UnalignedSliceViewer::adjustPropertiesToInput() {
    const VolumeBase* inputVolume = inport_.getData();

    if(inputVolume) {
        float diagonal = tgt::length(inputVolume->getBoundingBox(false).getBoundingBox(false).diagonal());
        distance_.set(0.0f);
        distance_.setMinValue(-diagonal);
        distance_.setMaxValue(diagonal);
    }
}

std::string UnalignedSliceViewer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = VolumeRenderer::generateHeader();

    header += "#define SLICE_TEXTURE_MODE_2D \n";
    header += "#define NUM_CHANNELS 1 \n";
    header += transferFunc_.get()->getShaderDefines();

    return header;
}

bool UnalignedSliceViewer::rebuildShader() {
    // do nothing if there is no shader at the moment
    if (!isInitialized() || !sliceShader_)
        return false;

    sliceShader_->setHeaders(generateHeader());
    return sliceShader_->rebuild();
}

}
