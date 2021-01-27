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

#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/utils/classificationmodes.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "tgt/textureunit.h"
#include "tgt/vector.h"

using tgt::mat4;
using tgt::TextureUnit;

namespace voreen {

const std::string VolumeRaycaster::loggerCat_("voreen.VolumeRaycaster");

VolumeRaycaster::VolumeRaycaster(bool addDefaultProps, bool usePortgroup, bool addPorts)
    : VolumeRenderer()
    , volumeInport_       (Port::INPORT,  "volumehandle.volumehandle", "Volume Input", false, Processor::INVALID_PROGRAM)
    , entryPort_          (Port::INPORT,  "image.entrypoints",         "Entry-points Input", false, Processor::INVALID_RESULT)
    , exitPort_           (Port::INPORT,  "image.exitpoints",          "Exit-points Input", false, Processor::INVALID_RESULT)
    , outport1_           (Port::OUTPORT, "image.output",              "Image Output", true, Processor::INVALID_PROGRAM)
    , outport2_           (Port::OUTPORT, "image.output1",             "Image1 Output", true, Processor::INVALID_PROGRAM)
    , outport3_           (Port::OUTPORT, "image.output2",             "Image2 Output", true, Processor::INVALID_PROGRAM)
    , internalRenderPort1_(Port::OUTPORT, "internalRenderPort1",       "Internal Render Port 1")
    , internalRenderPort2_(Port::OUTPORT, "internalRenderPort2",       "Internal Render Port 2")
    , internalRenderPort3_(Port::OUTPORT, "internalRenderPort3",       "Internal Render Port 3")
    , internalPortGroup_(true)
    , renderingQualityOption_  ("renderingQuality",      "Rendering Quality")
    , interactionQualityOption_("interQuality",          "Interaction Quality")
    , samplingRate_            ("samplingRate",          "Sampling Rate", 2.f, 0.01f, 20.f)
    , isoValue_                ("isoValue",              "Iso Value", 0.5f, 0.0f, 1.0f)
    , maskingMode_             ("masking",               "Masking", Processor::INVALID_PROGRAM)
    , gradientMode_            ("gradient",              "Gradient Calculation", Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , classificationMode_      ("classification",        "Classification", Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , shadeMode_               ("shading",               "Shading", Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , compositingMode1_        ("compositing",           "Compositing (Output 1)", Processor::INVALID_PROGRAM, false, Property::LOD_ADVANCED)
    , compositingMode2_        ("compositing2",          "Compositing (Output 2)", Processor::INVALID_PROGRAM, false, Property::LOD_ADVANCED)
    , compositingMode3_        ("compositing3",          "Compositing (Output 3)", Processor::INVALID_PROGRAM, false, Property::LOD_ADVANCED)
    , interactionCoarseness_   ("interactionCoarseness", "Interaction Coarseness", 3, 1, 8, Processor::VALID, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
    , interactionQuality_      ("interactionQuality",    "Interaction Step Size", 1.0f, 0.01f, 1.0f, Processor::VALID, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , gammaValue1_             ("gammaValue",            "Gamma Value (OP1)", 0, -1, 1, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , gammaValue2_             ("gammaValue1",           "Gamma Value (OP2)", 0, -1, 1, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , gammaValue3_             ("gammaValue2",           "Gamma Value (OP3)", 0, -1, 1, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , usePortgroup_(usePortgroup)
{
    // ports
    volumeInport_.addCondition(new PortConditionVolumeTypeGL());
    volumeInport_.showTextureAccessProperties(true);
    if (addPorts){
        addPort(volumeInport_);
        addPort(entryPort_);
        addPort(exitPort_);

        addPort(outport1_);
        addPrivateRenderPort(internalRenderPort1_);
        if (usePortgroup_){
            addPort(outport2_);
            addPort(outport3_);

            addPrivateRenderPort(internalRenderPort2_);
            addPrivateRenderPort(internalRenderPort3_);
        }
    }


    renderingQualityOption_.addSetting("Full",    &samplingRate_,                 10.f);
    renderingQualityOption_.addSetting<std::string>("Full",    &gradientMode_,                 "filtered");
    //renderingQualityOption_.addSetting<std::string>("Full",    &shadeMode_,                    "phong");

    renderingQualityOption_.addSetting("High",    &samplingRate_,                 4.f);
    renderingQualityOption_.addSetting<std::string>("High",    &gradientMode_,                 "sobel");
    //renderingQualityOption_.addSetting<std::string>("High",    &shadeMode_,                    "phong");

    renderingQualityOption_.addSetting("Medium",  &samplingRate_,                 2.f);
    renderingQualityOption_.addSetting<std::string>("Medium",  &gradientMode_,                 "central-differences");
    //renderingQualityOption_.addSetting<std::string>("Medium",  &shadeMode_,                    "phong");

    renderingQualityOption_.addSetting("Low",     &samplingRate_,                 1.f);
    renderingQualityOption_.addSetting<std::string>("Low",     &gradientMode_,                 "forward-differences");
    //renderingQualityOption_.addSetting<std::string>("Low",     &shadeMode_,                    "none");

    // quality settings are the same as for the rendering quality (and should be kept that way!)
    interactionQualityOption_.addSetting("Full",    &interactionQuality_,                      1.f);
    interactionQualityOption_.addSetting("Full",    &interactionCoarseness_,                   1);

    interactionQualityOption_.addSetting("High",    &interactionQuality_,                      1.f);
    interactionQualityOption_.addSetting("High",    &interactionCoarseness_,                   2);

    interactionQualityOption_.addSetting("Medium",  &interactionQuality_,                      0.75f);
    interactionQualityOption_.addSetting("Medium",  &interactionCoarseness_,                   4);

    interactionQualityOption_.addSetting("Low",    &interactionQuality_,                       0.5f);
    interactionQualityOption_.addSetting("Low",    &interactionCoarseness_,                    8);

    // initialization of the rendering properties
    // the properties are added in the respective subclasses
    maskingMode_.addOption("none", "none");
    maskingMode_.addOption("Segmentation", "Segmentation");

    gradientMode_.addOption("none",                 "none"                  );
    gradientMode_.addOption("forward-differences",  "Forward Differences"   );
    gradientMode_.addOption("central-differences",  "Central Differences"   );
    gradientMode_.addOption("sobel",                "Sobel"                 );
    gradientMode_.addOption("filtered",             "Filtered"              );
    gradientMode_.select("central-differences");

    ClassificationModes::fillProperty(&classificationMode_);

    GLSL::fillShadingModesProperty(shadeMode_);

    addCompositionModes(compositingMode1_);
    addCompositionModes(compositingMode2_);
    addCompositionModes(compositingMode3_);


    if(addDefaultProps) {
        addProperty(renderingQualityOption_);
        renderingQualityOption_.setGroupID("default_rc");
        addProperty(interactionQualityOption_);
        interactionQualityOption_.setGroupID("default_rc");
        addProperty(samplingRate_);
        samplingRate_.setGroupID("default_rc");
        addProperty(interactionCoarseness_);
        interactionCoarseness_.setGroupID("default_rc");
        addProperty(interactionQuality_);
        interactionQuality_.setGroupID("default_rc");
        renderingQualityOption_.selectMode("High");
        interactionQualityOption_.selectMode("Medium");
        setPropertyGroupGuiName("default_rc","Raycasting Settings");
    }
}

void VolumeRaycaster::initialize() {
    VolumeRenderer::initialize();
    QualityMode.addObserver(this);

    // load rescale shader program
    rescaleShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", generateHeader(), false);

    // if no camera property has been passed to the light property until now, try to find one automatically
    if(!lightPosition_.hasCamera())
        lightPosition_.getCamera();
    if (usePortgroup_){
        internalPortGroup_.initialize();
        internalPortGroup_.addPort(internalRenderPort1_);
        internalPortGroup_.addPort(internalRenderPort2_);
        internalPortGroup_.addPort(internalRenderPort3_);
    }
}

void VolumeRaycaster::deinitialize() {
    ShdrMgr.dispose(rescaleShader_);
    rescaleShader_ = 0;
    QualityMode.removeObserver(this);
    if (usePortgroup_){
        internalPortGroup_.removePort(internalRenderPort1_);
        internalPortGroup_.removePort(internalRenderPort2_);
        internalPortGroup_.removePort(internalRenderPort3_);
        internalPortGroup_.deinitialize();
    }
    VolumeRenderer::deinitialize();
}

/*
    further methods
*/
std::string VolumeRaycaster::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string headerSource = VolumeRenderer::generateHeader(version);

    // configure gradient calculation
    headerSource += "#define CALC_GRADIENT(volume, volumeStruct, samplePos) ";
    if (gradientMode_.isSelected("none"))
        headerSource += "(voxel.xyz-vec3(0.5))*2.0;\n";
    else if (gradientMode_.isSelected("forward-differences"))
        headerSource += "calcGradientAFD(volume, volumeStruct, samplePos);\n";
    else if (gradientMode_.isSelected("central-differences"))
        headerSource += "calcGradientA(volume, volumeStruct, samplePos);\n";
    else if (gradientMode_.isSelected("filtered"))
        headerSource += "calcGradientFiltered(volume, volumeStruct, samplePos);\n";
    else if (gradientMode_.isSelected("sobel"))
        headerSource += "calcGradientSobel(volume, volumeStruct, samplePos);\n";

    // configure classification
    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode_.get());

    // configure shading mode
    headerSource += GLSL::getShaderDefine(shadeMode_.get(), "APPLY_SHADING");

    // check for pre-integration
    bool usePreIntegrationClassification = false;
    const std::string& cString = classificationMode_.get();
    if (cString == "pre-integrated-gpu" || cString == "pre-integrated-cpu" || cString == "pre-integrated-approximation")
        usePreIntegrationClassification = true;

    std::string compositingMode1 =compositingMode1_.get();
    std::string compositingMode2 =compositingMode2_.get();
    std::string compositingMode3 =compositingMode3_.get();

    bool canUseEarlyTermination1 = ((compositingMode1 != "mip") && (compositingMode1 != "mop") && (compositingMode1 != "mida")) || !outport1_.isConnected();
    bool canUseEarlyTermination2 = ((compositingMode2 != "mip") && (compositingMode2 != "mop") && (compositingMode2 != "mida")) || !outport2_.isConnected();
    bool canUseEarlyTermination3 = ((compositingMode3 != "mip") && (compositingMode3 != "mop") && (compositingMode3 != "mida")) || !outport3_.isConnected();

    if (!usePortgroup_ && canUseEarlyTermination1){
        headerSource += "#define USE_EARLY_TERMINATION\n";
    }else if (usePortgroup_ && canUseEarlyTermination1 && canUseEarlyTermination2 && canUseEarlyTermination3){
        headerSource += "#define USE_EARLY_TERMINATION\n";
    }

    // configure compositing mode
    if (!usePortgroup_){
        headerSource += generateCompositionHeader("RC_APPLY_COMPOSITING", compositingMode1, usePreIntegrationClassification, outport1_.isConnected());
    }else{
        headerSource += generateCompositionHeader("RC_APPLY_COMPOSITING1", compositingMode1, usePreIntegrationClassification, outport1_.isConnected());
        headerSource += generateCompositionHeader("RC_APPLY_COMPOSITING2", compositingMode2, usePreIntegrationClassification, outport2_.isConnected());
        headerSource += generateCompositionHeader("RC_APPLY_COMPOSITING3", compositingMode3, usePreIntegrationClassification, outport3_.isConnected());
    }

    if (applyLightAttenuation_.get())
        headerSource += "#define PHONG_APPLY_ATTENUATION\n";

    return headerSource;
}

std::string VolumeRaycaster::generateCompositionHeader(std::string macroname, std::string compositionMode, bool usePreIntegrationClassification, bool isConnected){
    std::string headerSource;
    headerSource += "#define "+macroname+"(result, voxel, color, samplePos, gradient, t, samplingStepSize, tDepth, maxIntensity, gammaValue, isoValue) ";
    if (!isConnected)
        headerSource += "vec4(0, 0, 0, 1);\n"; // for early ray termination
    else if (compositionMode == "dvr" && usePreIntegrationClassification)
        headerSource += "compositeDVRNoOpacityCorrection(result, color, t, tDepth);\n";
    else if (compositionMode == "dvr")
        headerSource += "compositeDVR(result, color, t, samplingStepSize, tDepth);\n";
    else if (compositionMode == "mip")
        headerSource += "compositeMIP(result, voxel, color, t, tDepth, maxIntensity);\n";
    else if (compositionMode == "mop")
        headerSource += "compositeMOP(result, color, t, tDepth);\n";
    else if (compositionMode == "mida")
        headerSource += "compositeMIDA(result, voxel, color, maxIntensity, t, samplingStepSize, tDepth, gammaValue);\n";
    else if (compositionMode == "iso")
        headerSource += "compositeISO(result, color, t, tDepth, isoValue);\n";
    else if (compositionMode == "fhp")
        headerSource += "compositeFHP(samplePos, result, t, tDepth);\n";
    else if (compositionMode == "fhn")
        headerSource += "compositeFHN(gradient, result, t, tDepth);\n";
    return headerSource;
}

void VolumeRaycaster::setGlobalShaderParameters(tgt::Shader* shader, const tgt::Camera* camera, tgt::ivec2 screenDim) {
    VolumeRenderer::setGlobalShaderParameters(shader, camera, screenDim);

    shader->setIgnoreUniformLocationError(true);

    // provide values needed for correct depth value calculation
    if (camera) {
        float n = camera->getNearDist();
        float f = camera->getFarDist();
        shader->setUniform("const_to_z_e_1", 0.5f + 0.5f*((f+n)/(f-n)));
        shader->setUniform("const_to_z_e_2", ((f-n)/(f*n)));
        shader->setUniform("const_to_z_w_1", ((f*n)/(f-n)));
        shader->setUniform("const_to_z_w_2", 0.5f*((f+n)/(f-n))+0.5f);
    }

    shader->setIgnoreUniformLocationError(false);
}

void VolumeRaycaster::qualityModeChanged() {
    // make sure, we re-render with full resolution after switching out of interaction mode
    if (!QualityMode.isInteractionMode() && processedInInteraction_) {
        processedInInteraction_ = false;
        invalidate();
    }
}

void VolumeRaycaster::afterProcess() {
    VolumeRenderer::afterProcess();

    if (QualityMode.isInteractionMode()) {
        processedInInteraction_ = true;
    }
}

void VolumeRaycaster::invalidate(int inv) {
    VolumeRenderer::invalidate(inv);
}

float VolumeRaycaster::getSamplingStepSize(const VolumeBase* vol) {
    tgt::ivec3 dim = vol->getDimensions();

    // use dimension with the highest resolution for calculating the sampling step size
    float samplingStepSize = 1.f / (tgt::max(dim) * samplingRate_.get());

    if (QualityMode.isInteractionMode())
        samplingStepSize /= interactionQuality_.get();

    return samplingStepSize;
}

bool VolumeRaycaster::bindVolumes(tgt::Shader* shader, const std::vector<VolumeStruct> &volumes,
        const tgt::Camera* camera, const tgt::vec4& lightPosition) {

    if (!VolumeRenderer::bindVolumes(shader, volumes, camera, lightPosition))
        return false;

    shader->setIgnoreUniformLocationError(true);

    //TODO: This uses the first volume to set the step size. Could be changed so that step
    // size is calculated for each volume, but shaders need to be adapted as well to have volume
    // parameters available in ray setup and compositing. joerg
    if (volumes.size() > 0) {
        if (!volumes[0].volume_ || !volumes[0].volume_->getRepresentation<VolumeGL>() || !volumes[0].volume_->getRepresentation<VolumeGL>()->getTexture()) {
            LWARNING("No volume texture");
        }
        else {
            shader->setUniform("samplingStepSize_", getSamplingStepSize(volumes[0].volume_));
            LGL_ERROR;
        }
    }
    shader->setIgnoreUniformLocationError(false);

    return true;
}

void VolumeRaycaster::rescaleRendering(RenderPort& srcPort, RenderPort& destPort) {
    // activate and clear output render target
    destPort.activateTarget();
    destPort.clearTarget();

    // activate shader and set uniforms
    tgtAssert(rescaleShader_, "bypass shader not loaded");
    rescaleShader_->activate();
    setGlobalShaderParameters(rescaleShader_, 0, destPort.getSize());

    // bind input rendering to texture units
    tgt::TextureUnit colorUnit, depthUnit;
    srcPort.bindTextures(colorUnit.getEnum(), depthUnit.getEnum(), GL_LINEAR);
    srcPort.setTextureParameters(rescaleShader_, "texParams_");
    rescaleShader_->setUniform("colorTex_", colorUnit.getUnitNumber());
    rescaleShader_->setUniform("depthTex_", depthUnit.getUnitNumber());
    rescaleShader_->setUniform("useTexcoords", false);

    // render screen aligned quad
    renderQuad();

    // cleanup
    rescaleShader_->deactivate();
    destPort.deactivateTarget();
    TextureUnit::setZeroUnit();

    // check for OpenGL errors
    LGL_ERROR;
}

void VolumeRaycaster::addCompositionModes(StringOptionProperty &compositionProp)
{
    compositionProp.addOption("dvr", "DVR");
    compositionProp.addOption("mip", "MIP");
    compositionProp.addOption("mop", "MOP");
    compositionProp.addOption("mida", "MIDA");
    compositionProp.addOption("iso", "ISO");
    compositionProp.addOption("fhp", "FHP");
    compositionProp.addOption("fhn", "FHN");
}



} // namespace voreen
