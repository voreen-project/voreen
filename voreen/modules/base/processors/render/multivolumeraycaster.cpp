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

#include "multivolumeraycaster.h"

#include "tgt/textureunit.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/utils/classificationmodes.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include <sstream>

namespace voreen {

MultiVolumeRaycaster::MultiVolumeRaycaster()
    : VolumeRaycaster(true, true)
    //ports
    , volumeInport2_(Port::INPORT, "volume2", "Volume2 Input", false, Processor::INVALID_PROGRAM)
    , entryPort2_(Port::INPORT, "entryPort2", "Entry-points of Volume2", false, Processor::INVALID_RESULT)
    , exitPort2_(Port::INPORT, "exitPort2", "Exit-points of Volume2", false, Processor::INVALID_RESULT)
    , volumeInport3_(Port::INPORT, "volume3", "Volume3 Input", false, Processor::INVALID_PROGRAM)
    , entryPort3_(Port::INPORT, "entryPort3", "Entry-points of Volume3", false, Processor::INVALID_RESULT)
    , exitPort3_(Port::INPORT, "exitPort3", "Exit-points of Volume3", false, Processor::INVALID_RESULT)
    , volumeInport4_(Port::INPORT, "volume4", "Volume4 Input", false, Processor::INVALID_PROGRAM)
    , entryPort4_(Port::INPORT, "entryPort4", "Entry-points of Volume4", false, Processor::INVALID_RESULT)
    , exitPort4_(Port::INPORT, "exitPort4", "Exit-points of Volume4", false, Processor::INVALID_RESULT)
    , entryPortGlobal_(Port::INPORT, "entryPortGlobal", "Global Entry-points Input", false, Processor::INVALID_RESULT)
    , exitPortGlobal_(Port::INPORT, "exitPortGlobal", "Global Exit-points Input", false, Processor::INVALID_RESULT)
    //properties
    , shaderProp_("raycast.prg", "Raycasting Shader", "rc_multivolume.frag", "passthrough.vert")
    , classificationMode2_("classification2", "Classification 2", Processor::INVALID_PROGRAM)
    , classificationMode3_("classification3", "Classification 3", Processor::INVALID_PROGRAM)
    , classificationMode4_("classification4", "Classification 4", Processor::INVALID_PROGRAM)
    , shadeMode1_("shading1", "Shading 1", Processor::INVALID_PROGRAM)
    , shadeMode2_("shading2", "Shading 2", Processor::INVALID_PROGRAM)
    , shadeMode3_("shading3", "Shading 3", Processor::INVALID_PROGRAM)
    , shadeMode4_("shading4", "Shading 4", Processor::INVALID_PROGRAM)
    , transferFunc1_("transferFunction1", "Transfer Function 1")
    , transferFunc2_("transferFunction2", "Transfer Function 2")
    , transferFunc3_("transferFunction3", "Transfer Function 3")
    , transferFunc4_("transferFunction4", "Transfer Function 4")
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)), true)
{
    // ports
    volumeInport2_.addCondition(new PortConditionVolumeTypeGL());
    volumeInport3_.addCondition(new PortConditionVolumeTypeGL());
    volumeInport4_.addCondition(new PortConditionVolumeTypeGL());
    volumeInport_.showTextureAccessProperties(true);
    volumeInport2_.showTextureAccessProperties(true);
    volumeInport3_.showTextureAccessProperties(true);
    volumeInport4_.showTextureAccessProperties(true);
    addPort(volumeInport2_); addPort(entryPort2_); addPort(exitPort2_);
    addPort(volumeInport3_); addPort(entryPort3_); addPort(exitPort3_);
    addPort(volumeInport4_); addPort(entryPort4_); addPort(exitPort4_);
    addPort(entryPortGlobal_);
    addPort(exitPortGlobal_);


    // shader property
    addProperty(shaderProp_);

    addProperty(classificationMode_);
    ClassificationModes::fillProperty(&classificationMode2_);
    addProperty(classificationMode2_);
    ClassificationModes::fillProperty(&classificationMode3_);
    addProperty(classificationMode3_);
    ClassificationModes::fillProperty(&classificationMode4_);
    addProperty(classificationMode4_);

    // tf properties
    addProperty(transferFunc1_);
    addProperty(transferFunc2_);
    addProperty(transferFunc3_);
    addProperty(transferFunc4_);
    addProperty(camera_);

    // shading properties
    addProperty(gradientMode_);
    GLSL::fillShadingModesProperty(shadeMode1_);
    addProperty(shadeMode1_);
    GLSL::fillShadingModesProperty(shadeMode2_);
    addProperty(shadeMode2_);
    GLSL::fillShadingModesProperty(shadeMode3_);
    addProperty(shadeMode3_);
    GLSL::fillShadingModesProperty(shadeMode4_);
    addProperty(shadeMode4_);

    // compositing modes
    addProperty(compositingMode1_);
    addProperty(compositingMode2_);
    addProperty(compositingMode3_);
    addProperty(isoValue_);

    // lighting properties
    addProperty(lightPosition_);
    addProperty(lightAmbient_);
    addProperty(lightDiffuse_);
    addProperty(lightSpecular_);
    addProperty(materialShininess_);
    addProperty(applyLightAttenuation_);
    addProperty(lightAttenuation_);

    // assign lighting properties to property group
    lightPosition_.setGroupID("lighting");
    lightAmbient_.setGroupID("lighting");
    lightDiffuse_.setGroupID("lighting");
    lightSpecular_.setGroupID("lighting");
    materialShininess_.setGroupID("lighting");
    applyLightAttenuation_.setGroupID("lighting");
    lightAttenuation_.setGroupID("lighting");
    setPropertyGroupGuiName("lighting", "Lighting Parameters");

    // listen to changes of properties that influence the GUI state (i.e. visibility of other props)
    classificationMode_.onChange(MemberFunctionCallback<MultiVolumeRaycaster>(this, &MultiVolumeRaycaster::adjustPropertyVisibilities));
    shadeMode_.onChange(MemberFunctionCallback<MultiVolumeRaycaster>(this, &MultiVolumeRaycaster::adjustPropertyVisibilities));
    compositingMode1_.onChange(MemberFunctionCallback<MultiVolumeRaycaster>(this, &MultiVolumeRaycaster::adjustPropertyVisibilities));
    compositingMode1_.onChange(MemberFunctionCallback<MultiVolumeRaycaster>(this, &MultiVolumeRaycaster::adjustPropertyVisibilities));
    compositingMode2_.onChange(MemberFunctionCallback<MultiVolumeRaycaster>(this, &MultiVolumeRaycaster::adjustPropertyVisibilities));
    applyLightAttenuation_.onChange(MemberFunctionCallback<MultiVolumeRaycaster>(this, &MultiVolumeRaycaster::adjustPropertyVisibilities));
}

Processor* MultiVolumeRaycaster::create() const {
    return new MultiVolumeRaycaster();
}

void MultiVolumeRaycaster::initialize() {
    VolumeRaycaster::initialize();
    compile();

    LGL_ERROR;

    adjustPropertyVisibilities();
}

void MultiVolumeRaycaster::deinitialize() {

    LGL_ERROR;

    VolumeRaycaster::deinitialize();
}

void MultiVolumeRaycaster::compile() {
    shaderProp_.setHeader(generateHeader());
    shaderProp_.rebuild();
}

bool MultiVolumeRaycaster::isReady() const {
    if(!entryPortGlobal_.isReady() || !exitPortGlobal_.isReady())
        return false;

    if(!( ((volumeInport_.isReady() && entryPort_.isReady() && exitPort_.isReady()) ||
                (!volumeInport_.isReady() && !entryPort_.isReady() && !exitPort_.isReady()) ) &&
          ((volumeInport2_.isReady() && entryPort2_.isReady() && exitPort2_.isReady()) ||
                (!volumeInport2_.isReady() && !entryPort2_.isReady() && !exitPort2_.isReady()) ) &&
          ((volumeInport3_.isReady() && entryPort3_.isReady() && exitPort3_.isReady()) ||
                (!volumeInport3_.isReady() && !entryPort3_.isReady() && !exitPort3_.isReady()) ) &&
          ((volumeInport4_.isReady() && entryPort4_.isReady() && exitPort4_.isReady()) ||
            (!volumeInport4_.isReady() && !entryPort4_.isReady() && !exitPort4_.isReady())) ))
        return false;

    //check if at least one outport is connected:
    if (!outport1_.isReady() && !outport2_.isReady() && !outport3_.isReady())
        return false;

    return true;
}

void MultiVolumeRaycaster::beforeProcess() {
    VolumeRaycaster::beforeProcess();

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        PROFILING_BLOCK("compile");
        compile();
    }
    LGL_ERROR;

    transferFunc1_.setVolume(volumeInport_.getData());
    transferFunc2_.setVolume(volumeInport2_.getData());
    transferFunc3_.setVolume(volumeInport3_.getData());
    transferFunc4_.setVolume(volumeInport4_.getData());

    if(volumeInport_.hasChanged() || volumeInport2_.hasChanged() || volumeInport3_.hasChanged() || volumeInport4_.hasChanged()) {
        tgt::Bounds b;
        if(volumeInport_.hasData())
            b.addVolume(volumeInport_.getData()->getBoundingBox().getBoundingBox());
        if(volumeInport2_.hasData())
            b.addVolume(volumeInport2_.getData()->getBoundingBox().getBoundingBox());
        if(volumeInport3_.hasData())
            b.addVolume(volumeInport3_.getData()->getBoundingBox().getBoundingBox());
        if(volumeInport4_.hasData())
            b.addVolume(volumeInport4_.getData()->getBoundingBox().getBoundingBox());
        if(length(b.diagonal()) != 0.0)
            camera_.adaptInteractionToScene(b);
    }
}

float getVoxelSamplingStepSize(const VolumeBase* vol, float worldSamplingStepSize) {
    return tgt::min(worldSamplingStepSize * vol->getPhysicalToTextureMatrix().getScalingPart());
}

void MultiVolumeRaycaster::process() {
    // determine render size and activate internal port group
    const bool renderCoarse = QualityMode.isInteractionMode() && interactionCoarseness_.get() > 1;
    const tgt::svec2 renderSize = (renderCoarse ? (outport1_.getSize() / interactionCoarseness_.get()) : outport1_.getSize());
    internalPortGroup_.resize(renderSize);
    internalPortGroup_.activateTargets();
    internalPortGroup_.clearTargets();
    LGL_ERROR;

    // initialize shader
    tgt::Shader* raycastPrg = shaderProp_.getShader();
    raycastPrg->activate();

    // set common uniforms used by all shaders
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(raycastPrg, &cam, renderSize);
    LGL_ERROR;

    // bind entry/exit param textures
    tgt::TextureUnit entryGlobalUnit, entryGlobalDepthUnit, exitGlobalUnit, exitGlobalDepthUnit;
    entryPortGlobal_.bindTextures(entryGlobalUnit, entryGlobalDepthUnit, GL_NEAREST);
    raycastPrg->setUniform("entryGlobalPoints_", entryGlobalUnit.getUnitNumber());
    raycastPrg->setUniform("entryGlobalPointsDepth_", entryGlobalDepthUnit.getUnitNumber());
        //set entryExit parameters for each texture
    entryPortGlobal_.setTextureParameters(raycastPrg, "entryExitParameters_");
    exitPortGlobal_.bindTextures(exitGlobalUnit, exitGlobalDepthUnit, GL_NEAREST);
    raycastPrg->setUniform("exitGlobalPoints_", exitGlobalUnit.getUnitNumber());
    raycastPrg->setUniform("exitGlobalPointsDepth_", exitGlobalDepthUnit.getUnitNumber());
    LGL_ERROR;

    // bind the volumes and pass the necessary information to the shader
    tgt::TextureUnit volUnit1, volUnit2, volUnit3, volUnit4;
    tgt::TextureUnit entryUnit1, entryUnit2, entryUnit3, entryUnit4;
    tgt::TextureUnit exitUnit1, exitUnit2, exitUnit3, exitUnit4;
    std::vector<const VolumeBase*> volumes;
    std::vector<VolumeStruct> volumeTextures;
    if (volumeInport_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport_.getData(),
                    &volUnit1,
                    "volume1_","volumeStruct1_",
                    volumeInport_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport_.getTextureBorderIntensityProperty().get()),
                    volumeInport_.getTextureFilterModeProperty().getValue())
                );
        volumes.push_back(volumeInport_.getData());
        entryPort_.bindColorTexture(entryUnit1,GL_NEAREST);
        raycastPrg->setUniform("entryPoints1_", entryUnit1.getUnitNumber());
        exitPort_.bindColorTexture(exitUnit1,GL_NEAREST);
        raycastPrg->setUniform("exitPoints1_", exitUnit1.getUnitNumber());
    }
    if (volumeInport2_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport2_.getData(),
                    &volUnit2,
                    "volume2_","volumeStruct2_",
                    volumeInport2_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport2_.getTextureBorderIntensityProperty().get()),
                    volumeInport2_.getTextureFilterModeProperty().getValue())
                );
        volumes.push_back(volumeInport2_.getData());
        entryPort2_.bindColorTexture(entryUnit2,GL_NEAREST);
        raycastPrg->setUniform("entryPoints2_", entryUnit2.getUnitNumber());
        exitPort2_.bindColorTexture(exitUnit2,GL_NEAREST);
        raycastPrg->setUniform("exitPoints2_", exitUnit2.getUnitNumber());
    }
    if (volumeInport3_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport3_.getData(),
                    &volUnit3,
                    "volume3_","volumeStruct3_",
                    volumeInport3_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport3_.getTextureBorderIntensityProperty().get()),
                    volumeInport3_.getTextureFilterModeProperty().getValue())
                );
        volumes.push_back(volumeInport3_.getData());
        entryPort3_.bindColorTexture(entryUnit3,GL_NEAREST);
        raycastPrg->setUniform("entryPoints3_", entryUnit3.getUnitNumber());
        exitPort3_.bindColorTexture(exitUnit3,GL_NEAREST);
        raycastPrg->setUniform("exitPoints3_", exitUnit3.getUnitNumber());
    }
    if (volumeInport4_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport4_.getData(),
                    &volUnit4,
                    "volume4_","volumeStruct4_",
                    volumeInport4_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport4_.getTextureBorderIntensityProperty().get()),
                    volumeInport4_.getTextureFilterModeProperty().getValue())
                );
        volumes.push_back(volumeInport4_.getData());
        entryPort4_.bindColorTexture(entryUnit4,GL_NEAREST);
        raycastPrg->setUniform("entryPoints4_", entryUnit4.getUnitNumber());
        exitPort4_.bindColorTexture(exitUnit4,GL_NEAREST);
        raycastPrg->setUniform("exitPoints4_", exitUnit4.getUnitNumber());
    }
    bindVolumes(raycastPrg, volumeTextures, &cam, lightPosition_.get());
    LGL_ERROR;

    // determine ray step length in world coords. Must be done AFTER bindVolumes(), since that method also sets
    // the uniform "samplingStepSize_" using a texture space step length depending only on the first volume
    float samplingStepSizeWorld = 0.0f;
    if (volumeTextures.size() > 0) {
        float voxelSizeWorld = 999.f;
        float voxelSizeTexture = 999.f;
        for(size_t i=0; i<volumes.size(); ++i) {
            const VolumeBase* volume = volumes[i];
            tgtAssert(volume, "No volume");
            tgt::ivec3 volDim = volume->getDimensions();
            tgt::vec3 cubeSizeWorld = volume->getCubeSize() * volume->getPhysicalToWorldMatrix().getScalingPart();

            float tVoxelSizeWorld = tgt::max(cubeSizeWorld / tgt::vec3(volDim));
            if (tVoxelSizeWorld < voxelSizeWorld) {
                voxelSizeWorld = tVoxelSizeWorld;
                voxelSizeTexture = tgt::max(1.f / tgt::vec3(volDim));
            }
        }

        samplingStepSizeWorld = voxelSizeWorld / samplingRate_.get();
        float samplingStepSizeTexture = voxelSizeTexture / samplingRate_.get();

        if (QualityMode.isInteractionMode()) {
            samplingStepSizeWorld /= interactionQuality_.get();
            samplingStepSizeTexture /= interactionQuality_.get();
        }

        raycastPrg->setUniform("samplingStepSize_", samplingStepSizeWorld);
        if (compositingMode1_.isSelected("dvr")  ||
            (compositingMode2_.isSelected("dvr") && outport1_.isConnected()) ||
            (compositingMode3_.isSelected("dvr") && outport2_.isConnected()) ) {
            // adapts the compositing of the multivolume RC to the one of the singlevolume RC (see below).
            raycastPrg->setUniform("mvOpacityCorrectionFactor_", samplingStepSizeTexture / samplingStepSizeWorld);
        }
        LGL_ERROR;
    }
    LGL_ERROR;

    // bind transfer functions
    tgt::TextureUnit transferUnit1, transferUnit2, transferUnit3, transferUnit4;
    transferUnit1.activate();
    if (transferFunc1_.get() && volumeInport_.getData())
        ClassificationModes::bindTexture(classificationMode_.get(), transferFunc1_.get(), getVoxelSamplingStepSize(volumeInport_.getData(), samplingStepSizeWorld));

    transferUnit2.activate();
    if (transferFunc2_.get() && volumeInport2_.getData())
        ClassificationModes::bindTexture(classificationMode2_.get(), transferFunc2_.get(), getVoxelSamplingStepSize(volumeInport2_.getData(), samplingStepSizeWorld));

    transferUnit3.activate();
    if (transferFunc3_.get() && volumeInport3_.getData())
        ClassificationModes::bindTexture(classificationMode3_.get(), transferFunc3_.get(), getVoxelSamplingStepSize(volumeInport3_.getData(), samplingStepSizeWorld));

    transferUnit4.activate();
    if (transferFunc4_.get() && volumeInport4_.getData())
        ClassificationModes::bindTexture(classificationMode4_.get(), transferFunc4_.get(), getVoxelSamplingStepSize(volumeInport4_.getData(), samplingStepSizeWorld));

    LGL_ERROR;

    // pass raycaster specific uniforms to the shader
    if (compositingMode1_.get() ==  "iso" ||
        compositingMode2_.get() == "iso" ||
        compositingMode3_.get() == "iso")
        raycastPrg->setUniform("isoValue_", isoValue_.get());

    if(volumeInport_.isReady() && ClassificationModes::usesTransferFunction(classificationMode_.get()))
        transferFunc1_.get()->setUniform(raycastPrg, "transferFunc1_", "transferFuncTex1_", transferUnit1.getUnitNumber());
    if(volumeInport2_.isReady() && ClassificationModes::usesTransferFunction(classificationMode2_.get()))
        transferFunc2_.get()->setUniform(raycastPrg, "transferFunc2_", "transferFuncTex2_", transferUnit2.getUnitNumber());
    if(volumeInport3_.isReady() && ClassificationModes::usesTransferFunction(classificationMode3_.get()))
        transferFunc3_.get()->setUniform(raycastPrg, "transferFunc3_", "transferFuncTex3_", transferUnit3.getUnitNumber());
    if(volumeInport4_.isReady() && ClassificationModes::usesTransferFunction(classificationMode4_.get()))
        transferFunc4_.get()->setUniform(raycastPrg, "transferFunc4_", "transferFuncTex4_", transferUnit4.getUnitNumber());
    LGL_ERROR;


    // perform the actual raycasting by drawing a screen-aligned quad
    renderQuad();

    raycastPrg->deactivate();
    internalPortGroup_.deactivateTargets();
    LGL_ERROR;

    // copy over rendered images from internal port group to outports,
    // thereby rescaling them to outport dimensions
    if (outport1_.isConnected())
        rescaleRendering(internalRenderPort1_, outport1_);
    if (outport2_.isConnected())
        rescaleRendering(internalRenderPort2_, outport2_);
    if (outport2_.isConnected())
        rescaleRendering(internalRenderPort3_, outport3_);

    tgt::TextureUnit::setZeroUnit();
    LGL_ERROR;
}

std::string MultiVolumeRaycaster::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string headerSource = VolumeRaycaster::generateHeader();

    if(volumeInport_.isReady())
        headerSource += "#define VOLUME_1_ACTIVE\n";
    if(volumeInport2_.isReady())
        headerSource += "#define VOLUME_2_ACTIVE\n";
    if(volumeInport3_.isReady())
        headerSource += "#define VOLUME_3_ACTIVE\n";
    if(volumeInport4_.isReady())
        headerSource += "#define VOLUME_4_ACTIVE\n";

    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode_.get(), transferFunc1_.get(), "TF_SAMPLER_TYPE_1");
    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode2_.get(), transferFunc2_.get(), "TF_SAMPLER_TYPE_2");
    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode3_.get(), transferFunc3_.get(), "TF_SAMPLER_TYPE_3");
    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode4_.get(), transferFunc4_.get(), "TF_SAMPLER_TYPE_4");

    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode_.get(), "RC_APPLY_CLASSIFICATION");
    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode2_.get(), "RC_APPLY_CLASSIFICATION2");
    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode3_.get(), "RC_APPLY_CLASSIFICATION3");
    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode4_.get(), "RC_APPLY_CLASSIFICATION4");

    // configure shading mode
    headerSource += GLSL::getShaderDefine(shadeMode1_.get(), "APPLY_SHADING_1");
    headerSource += GLSL::getShaderDefine(shadeMode2_.get(), "APPLY_SHADING_2");
    headerSource += GLSL::getShaderDefine(shadeMode3_.get(), "APPLY_SHADING_3");
    headerSource += GLSL::getShaderDefine(shadeMode4_.get(), "APPLY_SHADING_4");

    // DVR opacity correction function adapting the MV compositing to the SVRC compositing,
    // used by the compositing macros below.
    // The adaption is necessary, because the multivolume RC samples in world space
    // instead of in texture space. Due to differing sampling base intervals, we would otherwise
    // still get correct compositing results, but the compositing would slightly differ from
    // the one performed by the SingleVolumeRaycaster.
    headerSource += "uniform float mvOpacityCorrectionFactor_;\n";
    headerSource += "vec4 mvOpacityCorrection(in vec4 color) {\n";
    headerSource += "  return vec4(color.rgb, 1.0 - pow(1.0-color.a, mvOpacityCorrectionFactor_));\n";
    headerSource += "}\n";

    internalPortGroup_.reattachTargets();
    headerSource += internalPortGroup_.generateHeader(shaderProp_.getShader());
    return headerSource;
}

void MultiVolumeRaycaster::adjustPropertyVisibilities() {
    bool useLighting = !shadeMode1_.isSelected("none") ||
                       !shadeMode2_.isSelected("none") ||
                       !shadeMode3_.isSelected("none") ||
                       !shadeMode4_.isSelected("none");
    setPropertyGroupVisible("lighting", useLighting);

    bool useIsovalue = (compositingMode1_.isSelected("iso")  ||
        compositingMode2_.isSelected("iso") ||
        compositingMode3_.isSelected("iso")   );
    isoValue_.setVisibleFlag(useIsovalue);

    lightAttenuation_.setVisibleFlag(applyLightAttenuation_.get());
}

} // namespace
