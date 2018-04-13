/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "singlevolumeslicer.h"
#include "tgt/glmath.h"
#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"

using tgt::TextureUnit;

namespace voreen {

SingleVolumeSlicer::SingleVolumeSlicer()
    : VolumeSlicer()
    , shaderProp_("slicing.shader", "Slicing shader", "sl_singlevolume.frag", "sl_singlevolume.vert")
    , outport_(Port::OUTPORT, "image.output", "Image Output", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER, GL_RGBA16F)
    , cameraHandler_("cameraHandler", "Camera Handler", &camera_)
{

    addInteractionHandler(cameraHandler_);

    addProperty(shaderProp_);

    addPort(outport_);
}

Processor* SingleVolumeSlicer::create() const {
    return new SingleVolumeSlicer();
}

void SingleVolumeSlicer::initialize() {
    VolumeSlicer::initialize();
    compile();
}

void SingleVolumeSlicer::compile() {
    shaderProp_.setHeader(generateHeader());
    shaderProp_.rebuild();
}

bool SingleVolumeSlicer::isReady() const {
    //check if all inports are connected:
    if(!volumeInport_.isReady())
        return false;

    //check if outport is connected:
    if (!outport_.isReady())
        return false;

    if(!shaderProp_.hasValidShader())
        return false;

    return true;
}

void SingleVolumeSlicer::process() {

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    // bind transfer function
    TextureUnit transferUnit;
    transferUnit.activate();
    if (transferFunc_.get())
        transferFunc_.get()->getTexture()->bind();

    transferFunc_.setVolume(volumeInport_.getData());

    // vector containing the volumes to bind; is passed to bindVolumes()
    std::vector<VolumeStruct> volumeTextures;

    // add main volume
    TextureUnit volUnit;
    volumeTextures.push_back(VolumeStruct(
        volumeInport_.getData(),
        &volUnit,
        "volume_","volumeStruct_")
    );

    // initialize slicing shader
    tgt::Shader* slicingPrg = shaderProp_.getShader();
    slicingPrg->activate();
    // fragment shader uniforms
    transferFunc_.get()->setUniform(slicingPrg, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());

    setupUniforms(slicingPrg);

    glDisable(GL_DEPTH_TEST);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(outport_.getSize()));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadMatrix(camera_.get().getViewMatrix());

    // set common uniforms used by all shaders
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(slicingPrg, &cam);
    // bind the volumes and pass the necessary information to the shader
    bindVolumes(slicingPrg, volumeTextures, &cam, lightPosition_.get());

    unsigned int numSlices = static_cast<unsigned int>(maxLength_ / sliceDistance_);

    slicingPrg->activate();
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_ONE);

    for (unsigned int curSlice=0; curSlice<numSlices; curSlice++) {
        IMode.begin(tgt::ImmediateMode::TRIANGLE_FAN);
            for (unsigned int curPoint=0; curPoint<6; curPoint++)
                IMode.vertex(tgt::ivec2(curPoint, curSlice));
        IMode.end();
    }

    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    slicingPrg->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();

    LGL_ERROR;
}

std::string SingleVolumeSlicer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string headerSource = VolumeRenderer::generateHeader();
    headerSource += transferFunc_.get()->getShaderDefines();
    return headerSource;
}



} // namespace
