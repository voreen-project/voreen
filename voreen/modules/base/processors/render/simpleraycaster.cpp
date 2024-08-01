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

#include "simpleraycaster.h"

#include "voreen/core/utils/voreenqualitymode.h"
#include "voreen/core/properties/cameraproperty.h"
#include "tgt/textureunit.h"

using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

SimpleRaycaster::SimpleRaycaster()
    : VolumeRaycaster()
    , shader_("shader", "Shader", "rc_simple.frag", "passthrough.vert")
    , transferFunc_("transferFunction", "Transfer Function", Processor::INVALID_RESULT)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)), true)
{
    volumeInport_.showTextureAccessProperties(true);

    addProperty(shader_);
    addProperty(transferFunc_);

    // camera is required for depth value calculation
    addProperty(camera_);
}

Processor* SimpleRaycaster::create() const {
    return new SimpleRaycaster();
}

void SimpleRaycaster::initialize() {
    VolumeRaycaster::initialize();
    rebuildShader();
}

void SimpleRaycaster::deinitialize() {
    VolumeRaycaster::deinitialize();
}

void SimpleRaycaster::beforeProcess() {
    VolumeRaycaster::beforeProcess();

    // rebuild shader, if changed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        rebuildShader();

    // assign volume to transfer function
    transferFunc_.setVolume(volumeInport_.getData());

    if(volumeInport_.hasChanged() && volumeInport_.getData())
        camera_.adaptInteractionToScene(volumeInport_.getData()->getBoundingBox().getBoundingBox());
    LGL_ERROR;
}

void SimpleRaycaster::process() {
    // determine, activate and clear render target
    const bool renderCoarse = QualityMode.isInteractionMode() && interactionCoarseness_.get() > 1;
    tgt::svec2 renderSize;
    if (renderCoarse) {
        renderSize = outport1_.getSize() / interactionCoarseness_.get();
        internalRenderPort1_.resize(renderSize);
    }
    else {
        renderSize = outport1_.getSize();
    }
    RenderPort& renderDestination = (renderCoarse ? internalRenderPort1_ : outport1_);
    renderDestination.activateTarget();
    renderDestination.clearTarget();

    // retrieve shader from shader property
    tgt::Shader* shader = shader_.getShader();
    if (!shader || !shader->isLinked()) {
        outport1_.deactivateTarget();
        return;
    }

    // activate shader and set common uniforms
    shader->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(shader, &cam, renderSize);
    LGL_ERROR;

    // bind entry and exit params and pass texture units to the shader
    TextureUnit entryUnit, entryDepthUnit, exitUnit, exitDepthUnit;
    entryPort_.bindTextures(entryUnit, entryDepthUnit, GL_NEAREST);
    shader->setUniform("entryPoints_", entryUnit.getUnitNumber());
    shader->setUniform("entryPointsDepth_", entryDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(shader, "entryParameters_");

    exitPort_.bindTextures(exitUnit, exitDepthUnit, GL_NEAREST);
    shader->setUniform("exitPoints_", exitUnit.getUnitNumber());
    shader->setUniform("exitPointsDepth_", exitDepthUnit.getUnitNumber());
    exitPort_.setTextureParameters(shader, "exitParameters_");
    LGL_ERROR;

    // bind volume texture and pass it to the shader
    std::vector<VolumeStruct> volumeTextures;
    TextureUnit volUnit;
     volumeTextures.push_back(VolumeStruct(
        volumeInport_.getData(),
        &volUnit,
        "volume_",
        "volumeStruct_",
        volumeInport_.getTextureClampModeProperty().getValue(),
        tgt::vec4(volumeInport_.getTextureBorderIntensityProperty().get()),
        volumeInport_.getTextureFilterModeProperty().getValue())
    );
    bindVolumes(shader, volumeTextures, &cam, lightPosition_.get());
    LGL_ERROR;

    // bind transfer function and pass it to the shader
    TextureUnit transferUnit;
    if (transferFunc_.get()) {
        transferUnit.activate();
        transferFunc_.get()->getTexture()->bind();
        transferFunc_.get()->setUniform(shader, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());
    }

    // render screen aligned quad
    renderQuad();

    // clean up
    shader->deactivate();
    renderDestination.deactivateTarget();

    // copy over and rescale image from internal render port to outport in coarseness mode
    if (renderCoarse) {
        rescaleRendering(renderDestination, outport1_);
    }

    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void SimpleRaycaster::rebuildShader() {
    shader_.setHeader(generateHeader());
    shader_.rebuild();
}

std::string SimpleRaycaster::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = VolumeRaycaster::generateHeader();

    header += transferFunc_.get()->getShaderDefines();
    return header;
}

} // namespace voreen
