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

#include "multislicerenderer.h"

#include "tgt/textureunit.h"

#include <sstream>

using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

void fillBlendingProp(StringOptionProperty& sop) {
    sop.addOption("maximumAlpha", "Maximum Alpha");
    sop.addOption("minimumAlpha", "Minimum Alpha");
    sop.addOption("replace", "Replace");
    sop.addOption("deactivate", "deactivate");
    sop.addOption("add", "Add");
    sop.addOption("subtract", "Subtract");
    sop.addOption("multiply", "Multiply");
    sop.addOption("divide", "Divide");
    sop.addOption("difference", "Difference");
    sop.addOption("decal", "Decal");
    sop.addOption("over", "Over");
    sop.addOption("setHue", "Set Hue");
    sop.addOption("blend", "Blend");
}

MultiSliceRenderer::MultiSliceRenderer()
    : VolumeRenderer()
    , raycastPrg_(0)
    , camera_("camera", "Camera", tgt::Camera())
    , interactionHandler_(0)
    , transferFunc1_("transferFunction1", "Transfer function 1")
    , transferFunc2_("transferFunction2", "Transfer function 2")
    , transferFunc3_("transferFunction3", "Transfer function 3")
    , transferFunc4_("transferFunction4", "Transfer function 4")
    , blendingMode1_("blendingMode1", "Blending 1", Processor::INVALID_PROGRAM)
    , blendingMode2_("blendingMode2", "Blending 2", Processor::INVALID_PROGRAM)
    , blendingMode3_("blendingMode3", "Blending 3", Processor::INVALID_PROGRAM)
    , blendingMode4_("blendingMode4", "Blending 4", Processor::INVALID_PROGRAM)
    , blendingFactor1_("blendingFactor1", "Blending Factor 1", 0.5f, 0.0f, 1.0f)
    , blendingFactor2_("blendingFactor2", "Blending Factor 2", 0.5f, 0.0f, 1.0f)
    , blendingFactor3_("blendingFactor3", "Blending Factor 3", 0.5f, 0.0f, 1.0f)
    , blendingFactor4_("blendingFactor4", "Blending Factor 4", 0.5f, 0.0f, 1.0f)
    , texFilterMode1_("textureFilterMode1_", "Texture Filtering 1")
    , texFilterMode2_("textureFilterMode2_", "Texture Filtering 2")
    , texFilterMode3_("textureFilterMode3_", "Texture Filtering 3")
    , texFilterMode4_("textureFilterMode4_", "Texture Filtering 4")
    , texClampMode1_("textureClampMode1_", "Texture Clamp 1")
    , texClampMode2_("textureClampMode2_", "Texture Clamp 2")
    , texClampMode3_("textureClampMode3_", "Texture Clamp 3")
    , texClampMode4_("textureClampMode4_", "Texture Clamp 4")
    , texBorderIntensity_("textureBorderIntensity", "Texture Border Intensity", 0.f)
    , volumeInport1_(Port::INPORT, "volume1", "Volume1 Input", false, Processor::INVALID_PROGRAM)
    , volumeInport2_(Port::INPORT, "volume2", "volume2 Input", false, Processor::INVALID_PROGRAM)
    , volumeInport3_(Port::INPORT, "volume3", "volume3 Input", false, Processor::INVALID_PROGRAM)
    , volumeInport4_(Port::INPORT, "volume4", "volume4 Input", false, Processor::INVALID_PROGRAM)
    , entryPort_(Port::INPORT, "image.entrypoints", "Entry-points Input")
    , outport_(Port::OUTPORT, "image.output", "Image Output", true, Processor::INVALID_PROGRAM)
{
    addProperty(camera_);

    addProperty(transferFunc1_);
    addProperty(transferFunc2_);
    addProperty(transferFunc3_);
    addProperty(transferFunc4_);

    fillBlendingProp(blendingMode1_);
    addProperty(blendingMode1_);
    addProperty(blendingFactor1_);

    fillBlendingProp(blendingMode2_);
    addProperty(blendingMode2_);
    addProperty(blendingFactor2_);

    fillBlendingProp(blendingMode3_);
    addProperty(blendingMode3_);
    addProperty(blendingFactor3_);

    fillBlendingProp(blendingMode4_);
    addProperty(blendingMode4_);
    addProperty(blendingFactor4_);

    // volume texture filtering
    texFilterMode1_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode1_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode1_.selectByKey("linear");
    addProperty(texFilterMode1_);
    texFilterMode2_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode2_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode2_.selectByKey("linear");
    addProperty(texFilterMode2_);
    texFilterMode3_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode3_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode3_.selectByKey("linear");
    addProperty(texFilterMode3_);
    texFilterMode4_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode4_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode4_.selectByKey("linear");
    addProperty(texFilterMode4_);

    // volume texture clamping
    //texClampMode1_.addOption("clamp",           "Clamp",             GL_CLAMP);
    texClampMode1_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode1_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode1_.selectByKey("clamp-to-edge");
    addProperty(texClampMode1_);
    //texClampMode2_.addOption("clamp",           "Clamp",             GL_CLAMP);
    texClampMode2_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode2_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode2_.selectByKey("clamp-to-edge");
    addProperty(texClampMode2_);
    //texClampMode3_.addOption("clamp",           "Clamp",             GL_CLAMP);
    texClampMode3_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode3_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode3_.selectByKey("clamp-to-edge");
    addProperty(texClampMode3_);
    //texClampMode4_.addOption("clamp",           "Clamp",             GL_CLAMP);
    texClampMode4_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode4_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode4_.selectByKey("clamp-to-edge");
    addProperty(texClampMode4_);
    addProperty(texBorderIntensity_);

    // assign texture access properties to property group
    texFilterMode1_.setGroupID("textureAccess");
    texFilterMode2_.setGroupID("textureAccess");
    texFilterMode3_.setGroupID("textureAccess");
    texFilterMode4_.setGroupID("textureAccess");
    texClampMode1_.setGroupID("textureAccess");
    texClampMode2_.setGroupID("textureAccess");
    texClampMode3_.setGroupID("textureAccess");
    texClampMode4_.setGroupID("textureAccess");
    texBorderIntensity_.setGroupID("textureAccess");
    setPropertyGroupGuiName("textureAccess", "Volume Texture Access");

    addPort(volumeInport1_);
    addPort(volumeInport2_);
    addPort(volumeInport3_);
    addPort(volumeInport4_);
    addPort(entryPort_);
    addPort(outport_);

    interactionHandler_ = new SliceCameraInteractionHandler("interactionHandler", "Camera Interaction", &camera_);
    addInteractionHandler(interactionHandler_);
}

MultiSliceRenderer::~MultiSliceRenderer() {
    delete interactionHandler_;
}

void MultiSliceRenderer::initialize() {
    VolumeRenderer::initialize();

    raycastPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "multislicerenderer.frag", generateHeader(), false);

    if (!raycastPrg_)
        throw VoreenException("Failed to load shaders");
}

void MultiSliceRenderer::deinitialize() {
    ShdrMgr.dispose(raycastPrg_);
    raycastPrg_ = 0;
    LGL_ERROR;

    VolumeRenderer::deinitialize();
}

void MultiSliceRenderer::compile() {
    raycastPrg_->setHeaders(generateHeader());
    raycastPrg_->rebuild();
}

bool MultiSliceRenderer::isReady() const {
    if(!entryPort_.isReady())
        return false;

    if(!volumeInport1_.isReady() && !volumeInport2_.isReady() && !volumeInport3_.isReady() && !volumeInport4_.isReady())
       return false;

    if (!outport_.isReady())
        return false;

    return true;
}

void MultiSliceRenderer::process() {
    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    // bind transfer function
    TextureUnit transferUnit, transferUnit2, transferUnit3, transferUnit4;
    transferUnit.activate();
    if (transferFunc1_.get())
        transferFunc1_.get()->getTexture()->bind();

    transferUnit2.activate();
    if (transferFunc2_.get())
        transferFunc2_.get()->getTexture()->bind();

    transferUnit3.activate();
    if (transferFunc3_.get())
        transferFunc3_.get()->getTexture()->bind();

    transferUnit4.activate();
    if (transferFunc4_.get())
        transferFunc4_.get()->getTexture()->bind();

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    transferFunc1_.setVolume(volumeInport1_.getData());
    transferFunc2_.setVolume(volumeInport2_.getData());
    transferFunc3_.setVolume(volumeInport3_.getData());
    transferFunc4_.setVolume(volumeInport4_.getData());

    TextureUnit entryUnit, entryDepthUnit;
    // bind entry params
    entryPort_.bindTextures(entryUnit.getEnum(), entryDepthUnit.getEnum());
    LGL_ERROR;

    // vector containing the volumes to bind; is passed to bindVolumes()
    std::vector<VolumeStruct> volumeTextures;

    TextureUnit volUnit, volUnit2, volUnit3, volUnit4;
    // add volumes
    if(volumeInport1_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport1_.getData(),
                    &volUnit,
                    "volume1_","volumeStruct1_",
                    texClampMode1_.getValue(),
                    tgt::vec4(texBorderIntensity_.get()),
                    texFilterMode1_.getValue())
                );
    }
    if(volumeInport2_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport2_.getData(),
                    &volUnit2,
                    "volume2_","volumeStruct2_",
                    texClampMode2_.getValue(),
                    tgt::vec4(texBorderIntensity_.get()),
                    texFilterMode2_.getValue())
                );
    }
    if(volumeInport3_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport3_.getData(),
                    &volUnit3,
                    "volume3_","volumeStruct3_",
                    texClampMode3_.getValue(),
                    tgt::vec4(texBorderIntensity_.get()),
                    texFilterMode3_.getValue())
                );
    }
    if(volumeInport4_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport4_.getData(),
                    &volUnit4,
                    "volume4_","volumeStruct4_",
                    texClampMode4_.getValue(),
                    tgt::vec4(texBorderIntensity_.get()),
                    texFilterMode4_.getValue())
                );
    }

    // initialize shader
    raycastPrg_->activate();

    // set common uniforms used by all shaders
    setGlobalShaderParameters(raycastPrg_);

    // bind the volumes and pass the necessary information to the shader
    bindVolumes(raycastPrg_, volumeTextures);

    // pass the remaining uniforms to the shader
    raycastPrg_->setUniform("entryPoints_", entryUnit.getUnitNumber());
    raycastPrg_->setUniform("entryPointsDepth_",entryDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(raycastPrg_, "entryParameters_");

    if(volumeInport1_.isReady()) {
        transferFunc1_.get()->setUniform(raycastPrg_, "transferFunc1_", "transferFuncTex1_", transferUnit.getUnitNumber());

        raycastPrg_->setIgnoreUniformLocationError(true);
        raycastPrg_->setUniform("blendingFactor1_", blendingFactor1_.get());
        raycastPrg_->setIgnoreUniformLocationError(false);
    }
    if(volumeInport2_.isReady()) {
        transferFunc2_.get()->setUniform(raycastPrg_, "transferFunc2_", "transferFuncTex2_", transferUnit2.getUnitNumber());

        raycastPrg_->setIgnoreUniformLocationError(true);
        raycastPrg_->setUniform("blendingFactor2_", blendingFactor2_.get());
        raycastPrg_->setIgnoreUniformLocationError(false);
    }
    if(volumeInport3_.isReady()) {
        transferFunc3_.get()->setUniform(raycastPrg_, "transferFunc3_", "transferFuncTex3_", transferUnit3.getUnitNumber());

        raycastPrg_->setIgnoreUniformLocationError(true);
        raycastPrg_->setUniform("blendingFactor3_", blendingFactor3_.get());
        raycastPrg_->setIgnoreUniformLocationError(false);
    }
    if(volumeInport4_.isReady()) {
        transferFunc4_.get()->setUniform(raycastPrg_, "transferFunc4_", "transferFuncTex4_", transferUnit4.getUnitNumber());

        raycastPrg_->setIgnoreUniformLocationError(true);
        raycastPrg_->setUniform("blendingFactor4_", blendingFactor4_.get());
        raycastPrg_->setIgnoreUniformLocationError(false);
    }

    LGL_ERROR;

    renderQuad();

    raycastPrg_->deactivate();

    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;

    outport_.deactivateTarget();
}

std::string MultiSliceRenderer::generateHeader(const tgt::GpuCapabilities::GlVersion*) {
    std::string headerSource = VolumeRenderer::generateHeader();

    if(volumeInport1_.isReady())
        headerSource += "#define VOLUME_1_ACTIVE\n";
    if(volumeInport2_.isReady())
        headerSource += "#define VOLUME_2_ACTIVE\n";
    if(volumeInport3_.isReady())
        headerSource += "#define VOLUME_3_ACTIVE\n";
    if(volumeInport4_.isReady())
        headerSource += "#define VOLUME_4_ACTIVE\n";

    headerSource += "#define TF_SAMPLER_TYPE_1 " + transferFunc1_.get()->getSamplerType() + "\n";
    headerSource += "#define TF_SAMPLER_TYPE_2 " + transferFunc2_.get()->getSamplerType() + "\n";
    headerSource += "#define TF_SAMPLER_TYPE_3 " + transferFunc3_.get()->getSamplerType() + "\n";
    headerSource += "#define TF_SAMPLER_TYPE_4 " + transferFunc4_.get()->getSamplerType() + "\n";

    headerSource += "#define BLENDINGMODE1(result, color, factor) " + blendingMode1_.get() + "(result, color, factor);\n";
    headerSource += "#define BLENDINGMODE2(result, color, factor) " + blendingMode2_.get() + "(result, color, factor);\n";
    headerSource += "#define BLENDINGMODE3(result, color, factor) " + blendingMode3_.get() + "(result, color, factor);\n";
    headerSource += "#define BLENDINGMODE4(result, color, factor) " + blendingMode4_.get() + "(result, color, factor);\n";

    return headerSource;
}

} // namespace
