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

#include "singlevolumeraycaster.h"

#include "tgt/textureunit.h"

#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "voreen/core/utils/classificationmodes.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include <sstream>

using tgt::vec3;
using tgt::TextureUnit;

#define LIGHT_POSITION_MAX_DIST_MULTIPLIER 10

namespace voreen {

const std::string SingleVolumeRaycaster::loggerCat_("voreen.SingleVolumeRaycaster");

SingleVolumeRaycaster::SingleVolumeRaycaster()
    : VolumeRaycaster(true, true)
    , transFuncTypeProp_("transFuncType", "TF Mode")
    , transFunc1DProp1_("transferFunction", "Transfer Function 1D (Channel 1)")
    , transFunc1DProp2_("transferFunction1", "Transfer Function 1D (Channel 2)")
    , transFunc1DProp3_("transferFunction2", "Transfer Function 1D (Channel 3)")
    , transFunc1DProp4_("transferFunction3", "Transfer Function 1D (Channel 4)")
    , transFunc2DProp1_("transFunc2D", "Transfer Function 2D (Channel 1)")
    , transFunc2DProp2_("transFunc2D_2", "Transfer Function 2D (Channel 2)")
    , transFunc2DProp3_("transFunc2D_3", "Transfer Function 2D (Channel 3)")
    , transFunc2DProp4_("transFunc2D_4", "Transfer Function 2D (Channel 4)")
    , transFunc1DGaussianProp1_("transFuncGaussian1", "Transfer Function 1D Gaussian (Channel 1)")
    , transFunc1DGaussianProp2_("transFuncGaussian2", "Transfer Function 1D Gaussian (Channel 2)")
    , transFunc1DGaussianProp3_("transFuncGaussian3", "Transfer Function 1D Gaussian (Channel 3)")
    , transFunc1DGaussianProp4_("transFuncGaussian4", "Transfer Function 1D Gaussian (Channel 4)")
    , channelShift1_("channelShift1", "Channel Shift (1)", tgt::vec3::zero, tgt::vec3(-100), tgt::vec3(+100))
    , channelShift2_("channelShift2", "Channel Shift (2)", tgt::vec3::zero, tgt::vec3(-100), tgt::vec3(+100))
    , channelShift3_("channelShift3", "Channel Shift (3)", tgt::vec3::zero, tgt::vec3(-100), tgt::vec3(+100))
    , channelShift4_("channelShift4", "Channel Shift (4)", tgt::vec3::zero, tgt::vec3(-100), tgt::vec3(+100))
    , enableChannelShift_("enableChannelShift", "Enable Channel Shift", false,  Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , shaderProp_("raycast.prg", "Raycasting Shader", "rc_singlevolume.frag", "passthrough.vert", "", Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)), true)

{

    volumeInport_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::volumeInportOnChange));
    addProperty(shaderProp_);

    // shading / classification props
    addProperty(transFuncTypeProp_);
    addProperty(transFunc1DProp1_);
    addProperty(transFunc1DProp2_);
    addProperty(transFunc1DProp3_);
    addProperty(transFunc1DProp4_);
    addProperty(transFunc2DProp1_);
    addProperty(transFunc2DProp2_);
    addProperty(transFunc2DProp3_);
    addProperty(transFunc2DProp4_);
        transFuncTypeProp_.addTransFuncProperty(&transFunc1DProp1_);
        transFuncTypeProp_.addTransFuncProperty(&transFunc2DProp1_);
        transFuncTypeProp_.addTransFuncProperty(&transFunc1DGaussianProp1_);
    transFuncTypeProp_.setGroupID("tf");
    transFunc1DProp1_.setGroupID("tf");
    transFunc1DProp2_.setGroupID("tf");
    transFunc1DProp3_.setGroupID("tf");
    transFunc1DProp4_.setGroupID("tf");
    transFunc2DProp1_.setGroupID("tf");
    transFunc2DProp2_.setGroupID("tf");
    transFunc2DProp3_.setGroupID("tf");
    transFunc2DProp4_.setGroupID("tf");

    addProperty(transFunc1DGaussianProp1_);
    addProperty(transFunc1DGaussianProp2_);
    addProperty(transFunc1DGaussianProp3_);
    addProperty(transFunc1DGaussianProp4_);
    transFunc1DGaussianProp1_.setGroupID("tf");
    transFunc1DGaussianProp2_.setGroupID("tf");
    transFunc1DGaussianProp3_.setGroupID("tf");
    transFunc1DGaussianProp4_.setGroupID("tf");

    setPropertyGroupGuiName("tf","Transfer Function");
    addProperty(camera_);
    addProperty(gradientMode_);
    addProperty(classificationMode_);
    addProperty(shadeMode_);

    gammaValue1_.setTracking(true);
    addProperty(gammaValue1_);
    gammaValue2_.setTracking(true);
    addProperty(gammaValue2_);
    gammaValue3_.setTracking(true);
    addProperty(gammaValue3_);

    // compositing modes
    addProperty(compositingMode1_);
    addProperty(compositingMode2_);
    addProperty(compositingMode3_);


    addProperty(isoValue_);

    // lighting
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

    addProperty(enableChannelShift_);
    addProperty(channelShift1_);
    addProperty(channelShift2_);
    addProperty(channelShift3_);
    addProperty(channelShift4_);

    enableChannelShift_.setGroupID("channelShift");
    channelShift1_.setGroupID("channelShift");
    channelShift2_.setGroupID("channelShift");
    channelShift3_.setGroupID("channelShift");
    channelShift4_.setGroupID("channelShift");
    setPropertyGroupGuiName("channelShift", "Channel Shift");


    // listen to changes of properties that influence the GUI state (i.e. visibility of other props)
    classificationMode_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    shadeMode_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    compositingMode1_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    compositingMode2_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    compositingMode3_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    applyLightAttenuation_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    transFuncTypeProp_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    enableChannelShift_.onChange(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
    volumeInport_.onNewData(MemberFunctionCallback<SingleVolumeRaycaster>(this, &SingleVolumeRaycaster::adjustPropertyVisibilities));
}

Processor* SingleVolumeRaycaster::create() const {
    return new SingleVolumeRaycaster();
}

void SingleVolumeRaycaster::initialize() {
    VolumeRaycaster::initialize();
    compile();

    adjustPropertyVisibilities();
}

void SingleVolumeRaycaster::deinitialize() {
    VolumeRaycaster::deinitialize();
}

void SingleVolumeRaycaster::compile() {
    shaderProp_.setHeader(generateHeader());
    shaderProp_.rebuild();
}

bool SingleVolumeRaycaster::isReady() const {
    //check if all inports are connected
    if(!entryPort_.isReady() || !exitPort_.isReady() || !volumeInport_.isReady())
        return false;

    //check if at least one outport is connected
    if (!outport1_.isReady() && !outport2_.isReady() && !outport3_.isReady())
        return false;

    if(!shaderProp_.hasValidShader() && !shaderProp_.requiresRebuild())
        return false;

    return true;
}

void SingleVolumeRaycaster::beforeProcess() {
    VolumeRaycaster::beforeProcess();

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        PROFILING_BLOCK("compile");
        compile();
    }
    LGL_ERROR;

    //HACK: 2D TFs use FBOs to update the texture, we trigger the update here to prevent conflicts in process()
    if(!transFunc2DProp1_.get()->isTextureValid())
        transFunc2DProp1_.get()->getTexture(); // this forces an update
    if(!transFunc2DProp2_.get()->isTextureValid())
        transFunc2DProp2_.get()->getTexture(); // this forces an update
    if(!transFunc2DProp3_.get()->isTextureValid())
        transFunc2DProp3_.get()->getTexture(); // this forces an update
    if(!transFunc2DProp4_.get()->isTextureValid())
        transFunc2DProp4_.get()->getTexture(); // this forces an update

    //if the pre-integration table is computed on the gpu, this must be done before the rendering process
    if (classificationMode_.get() == "pre-integrated-gpu"){
        transFunc1DProp1_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
        transFunc1DGaussianProp1_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
    }
    if (classificationMode_.get() == "pre-integrated-gpu"){
        transFunc1DProp2_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
        transFunc1DGaussianProp2_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
    }
    if (classificationMode_.get() == "pre-integrated-gpu"){
        transFunc1DProp3_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
        transFunc1DGaussianProp3_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
    }
    if (classificationMode_.get() == "pre-integrated-gpu"){
        transFunc1DProp4_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
        transFunc1DGaussianProp4_.get()->getPreIntegrationTable(getSamplingStepSize(volumeInport_.getData()), 0, false, true)->getTexture();
    }
}

void SingleVolumeRaycaster::process() {

    LGL_ERROR;

    int numChannels = volumeInport_.getData()->getNumChannels();

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
    LGL_ERROR;

    // set common uniforms used by all shaders
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(raycastPrg, &cam, renderSize);
    LGL_ERROR;

    if (enableChannelShift_.get()){
        vec3 dim = vec3(volumeInport_.getData()->getDimensions());

        raycastPrg->setUniform("channelShift1_", channelShift1_.get()/dim);
        if (numChannels >= 2)
            raycastPrg->setUniform("channelShift2_", channelShift2_.get()/dim);
        if (numChannels >= 3)
            raycastPrg->setUniform("channelShift3_", channelShift3_.get()/dim);
        if (numChannels >= 4)
            raycastPrg->setUniform("channelShift4_", channelShift4_.get()/dim);
    }

    // bind entry/exit param textures
    tgt::TextureUnit entryUnit, entryDepthUnit, exitUnit, exitDepthUnit;
    entryPort_.bindTextures(entryUnit, entryDepthUnit, GL_NEAREST);
    raycastPrg->setUniform("entryPoints_", entryUnit.getUnitNumber());
    raycastPrg->setUniform("entryPointsDepth_", entryDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(raycastPrg, "entryParameters_");

    exitPort_.bindTextures(exitUnit, exitDepthUnit, GL_NEAREST);
    raycastPrg->setUniform("exitPoints_", exitUnit.getUnitNumber());
    raycastPrg->setUniform("exitPointsDepth_", exitDepthUnit.getUnitNumber());
    exitPort_.setTextureParameters(raycastPrg, "exitParameters_");
    LGL_ERROR;

    // bind the volumes and pass the necessary information to the shader
    TextureUnit volUnit;
    std::vector<VolumeStruct> volumeTextures;
    volumeTextures.push_back(VolumeStruct(
        volumeInport_.getData(),
        &volUnit,
        "volume_","volumeStruct_",
        volumeInport_.getTextureClampModeProperty().getValue(),
        tgt::vec4(volumeInport_.getTextureBorderIntensityProperty().get()),
        volumeInport_.getTextureFilterModeProperty().getValue())
    );
    bindVolumes(raycastPrg, volumeTextures, &cam, lightPosition_.get());
    LGL_ERROR;

    // bind transfer function
    tgt::TextureUnit transferUnit1, transferUnit2, transferUnit3, transferUnit4;
    //
    switch(transFuncTypeProp_.getValue()) {

    case TFT_1DKEYS:
        transferUnit1.activate();
        ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DProp1_.get(), getSamplingStepSize(volumeInport_.getData()));
        if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
            transFunc1DProp1_.get()->setUniform(raycastPrg, "transferFunc1_", "transferFuncTex1_", transferUnit1.getUnitNumber());
        }
        if (numChannels >= 2){
            transferUnit2.activate();
            ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DProp2_.get(), getSamplingStepSize(volumeInport_.getData()));
            if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
                transFunc1DProp2_.get()->setUniform(raycastPrg, "transferFunc2_", "transferFuncTex2_", transferUnit2.getUnitNumber());
            }
        }
        if (numChannels >= 3){
            transferUnit3.activate();
            ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DProp3_.get(), getSamplingStepSize(volumeInport_.getData()));
            if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
                transFunc1DProp3_.get()->setUniform(raycastPrg, "transferFunc3_", "transferFuncTex3_", transferUnit3.getUnitNumber());
            }
        }
        if (numChannels >= 4){
            transferUnit4.activate();
            ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DProp4_.get(), getSamplingStepSize(volumeInport_.getData()));
            if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
                transFunc1DProp4_.get()->setUniform(raycastPrg, "transferFunc4_", "transferFuncTex4_", transferUnit4.getUnitNumber());
            }
        }
        break;
    case TFT_1DGAUSSIAN:
        transferUnit1.activate();
        ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DGaussianProp1_.get(), getSamplingStepSize(volumeInport_.getData()));
        if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
            transFunc1DGaussianProp1_.get()->setUniform(raycastPrg, "transferFunc1_", "transferFuncTex1_", transferUnit1.getUnitNumber());
        }
        if (numChannels >= 2){
            transferUnit2.activate();
            ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DGaussianProp2_.get(), getSamplingStepSize(volumeInport_.getData()));
            if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
                transFunc1DGaussianProp2_.get()->setUniform(raycastPrg, "transferFunc2_", "transferFuncTex2_", transferUnit2.getUnitNumber());
            }
        }
        if (numChannels >= 3){
            transferUnit3.activate();
            ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DGaussianProp3_.get(), getSamplingStepSize(volumeInport_.getData()));
            if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
                transFunc1DGaussianProp3_.get()->setUniform(raycastPrg, "transferFunc3_", "transferFuncTex3_", transferUnit3.getUnitNumber());
            }
        }
        if (numChannels >= 4){
            transferUnit4.activate();
            ClassificationModes::bindTexture(classificationMode_.get(), transFunc1DGaussianProp4_.get(), getSamplingStepSize(volumeInport_.getData()));
            if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
                transFunc1DGaussianProp4_.get()->setUniform(raycastPrg, "transferFunc4_", "transferFuncTex4_", transferUnit4.getUnitNumber());
            }
        }
        break;
    case TFT_2DPRIMITIVES:
        /************************************************************************/
        /*                                                                      */
        /************************************************************************/
        transferUnit1.activate();
        ClassificationModes::bindTexture(classificationMode_.get(), transFunc2DProp1_.get(), getSamplingStepSize(volumeInport_.getData()));
        if (ClassificationModes::usesTransferFunction(classificationMode_.get())){
            transFunc2DProp1_.get()->setUniform(raycastPrg, "transferFunc1_", "transferFuncTex1_", transferUnit1.getUnitNumber());
        }
        break;
    default:
        tgtAssert(false,"Unknown trans func type!");
        LERROR("Unknown trans func type!");
    }

    LGL_ERROR;

    // pass remaining uniforms to shader
    if (compositingMode1_.isSelected("iso")  ||
        compositingMode1_.isSelected("iso") ||
        compositingMode2_.isSelected("iso") )
        raycastPrg->setUniform("isoValue_", isoValue_.get());

    if (compositingMode1_.isSelected("mida"))
        raycastPrg->setUniform("gammaValue1_", gammaValue1_.get());

    if (compositingMode1_.isSelected("mida"))
        raycastPrg->setUniform("gammaValue2_", gammaValue2_.get());

    if (compositingMode2_.isSelected("mida"))
        raycastPrg->setUniform("gammaValue3_", gammaValue3_.get());

    LGL_ERROR;

    // perform the actual raycasting by drawing a screen-aligned quad
    {
        PROFILING_BLOCK("raycasting");
        renderQuad();
    }

    raycastPrg->deactivate();
    internalPortGroup_.deactivateTargets();
    LGL_ERROR;

    // copy over rendered images from internal port group to outports,
    // thereby rescaling them to outport dimensions
    if (outport1_.isConnected())
        rescaleRendering(internalRenderPort1_, outport1_);
    if (outport2_.isConnected())
        rescaleRendering(internalRenderPort2_, outport2_);
    if (outport3_.isConnected())
        rescaleRendering(internalRenderPort3_, outport3_);

    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

std::string SingleVolumeRaycaster::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string headerSource = VolumeRaycaster::generateHeader();
    //switch between 1D and 2D
    switch(transFuncTypeProp_.getValue()) {
    case TFT_1DKEYS:
        headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode_.get(), transFunc1DProp1_.get(), "TF_SAMPLER_TYPE");
        break;
    case TFT_1DGAUSSIAN:
        headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode_.get(), transFunc1DGaussianProp1_.get());
        break;
    case TFT_2DPRIMITIVES:
        headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode_.get(), transFunc2DProp1_.get());
        break;
    default:
        tgtAssert(false,"Unknown trans func type!");
        LERROR("Unknown trans func type!");
    }
    {
        int numChannels = volumeInport_.hasData() ? volumeInport_.getData()->getNumChannels() : 1;
        std::stringstream numChannels_str;
        numChannels_str << numChannels;
        headerSource += "#define NUM_CHANNELS "+numChannels_str.str()+"\n";
    }
    if (compositingMode1_.get() == "mip"){
        headerSource += "#define SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT1\n";
    }
    if (compositingMode2_.get() == "mip"){
        headerSource += "#define SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT2\n";
    }
    if (compositingMode3_.get() == "mip"){
        headerSource += "#define SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT3\n";
    }
    if (enableChannelShift_.get()){
        std::stringstream ss;
        headerSource += "#define ENABLE_CHANNEL_SHIFT\n";
    }
    // check for pre-integration
    /*bool usePreIntegrationClassification = false;
    const std::string& cString = classificationMode_.get();
    if (cString == "pre-integrated-gpu" || cString == "pre-integrated-cpu" || cString == "pre-integrated-approximation")
        usePreIntegrationClassification = true;*/

    internalPortGroup_.reattachTargets();
    headerSource += internalPortGroup_.generateHeader(shaderProp_.getShader());
    return headerSource;
}

//---------------------------------------------------------------------
//      Callbacks
//---------------------------------------------------------------------
void SingleVolumeRaycaster::adjustPropertyVisibilities() {
    bool useLighting = !shadeMode_.isSelected("none");
    setPropertyGroupVisible("lighting", useLighting);

    bool useIsovalue = (compositingMode1_.isSelected("iso")  ||
                        compositingMode2_.isSelected("iso") ||
                        compositingMode3_.isSelected("iso")   );
    isoValue_.setVisibleFlag(useIsovalue);

    bool tf1d = transFuncTypeProp_.getValue() == TFT_1DKEYS;
    bool tf2d = transFuncTypeProp_.getValue() == TFT_2DPRIMITIVES;
    bool tf1dgausian = transFuncTypeProp_.getValue() == TFT_1DGAUSSIAN;

    int channels = (volumeInport_.getData()!=nullptr) ? volumeInport_.getData()->getNumChannels() : 1;
    transFunc1DProp2_.setVisibleFlag(channels >= 2 && tf1d);
    transFunc1DProp3_.setVisibleFlag(channels >= 3 && tf1d);
    transFunc1DProp4_.setVisibleFlag(channels >= 4 && tf1d);

    transFunc2DProp2_.setVisibleFlag(channels >= 2 && tf2d);
    transFunc2DProp3_.setVisibleFlag(channels >= 3 && tf2d);
    transFunc2DProp4_.setVisibleFlag(channels >= 4 && tf2d);

    transFunc1DGaussianProp2_.setVisibleFlag(channels >= 2 && tf1dgausian);
    transFunc1DGaussianProp3_.setVisibleFlag(channels >= 3 && tf1dgausian);
    transFunc1DGaussianProp4_.setVisibleFlag(channels >= 4 && tf1dgausian);

    bool channelShift = enableChannelShift_.get();
    channelShift1_.setVisibleFlag(channelShift && channels >= 1);
    channelShift2_.setVisibleFlag(channelShift && channels >= 2);
    channelShift3_.setVisibleFlag(channelShift && channels >= 3);
    channelShift4_.setVisibleFlag(channelShift && channels >= 4);


    lightAttenuation_.setVisibleFlag(applyLightAttenuation_.get());

    gammaValue1_.setVisibleFlag(compositingMode1_.isSelected("mida"));
    gammaValue2_.setVisibleFlag(compositingMode2_.isSelected("mida"));
    gammaValue3_.setVisibleFlag(compositingMode3_.isSelected("mida"));

    //2d tf does not support pre-integration
    if(transFuncTypeProp_.getValue() == TFT_2DPRIMITIVES
        && !(classificationMode_.get() == "transfer-function")) {
        LWARNING("Pre-integration is not allowed with 2d transfer functions. Switching to normal transfer function mode.");
        classificationMode_.select("transfer-function");
    }


    if((compositingMode1_.get() == "mip" || compositingMode2_.get() == "mip" || compositingMode3_.get() == "mip")
        && !(classificationMode_.get() == "transfer-function")) {
            LWARNING("Pre-integration is not allowed with Mip rendering. Switching to normal transfer function mode.");
            classificationMode_.select("transfer-function");
    }
}

void SingleVolumeRaycaster::volumeInportOnChange() {
    //adjust camera to scene
    if(volumeInport_.hasData()) {
        auto bb = volumeInport_.getData()->getBoundingBox().getBoundingBox();
        camera_.adaptInteractionToScene(bb, tgt::min(volumeInport_.getData()->getSpacing()));

        lightPosition_.setMaxDist(LIGHT_POSITION_MAX_DIST_MULTIPLIER * tgt::length(bb.diagonal()));
    }

    //update tf properties
    int channels = volumeInport_.hasData() ? volumeInport_.getData()->getNumChannels() : 1;

    transFunc1DProp1_.setVolume(volumeInport_.getData(), 0);
    if (channels >= 2) transFunc1DProp2_.setVolume(volumeInport_.getData(), 1);
    if (channels >= 3) transFunc1DProp3_.setVolume(volumeInport_.getData(), 2);
    if (channels >= 4) transFunc1DProp4_.setVolume(volumeInport_.getData(), 3);
    transFunc2DProp1_.setVolume(volumeInport_.getData(), 0);
    if (channels >= 2) transFunc2DProp2_.setVolume(volumeInport_.getData(), 1);
    if (channels >= 3) transFunc2DProp3_.setVolume(volumeInport_.getData(), 2);
    if (channels >= 4) transFunc2DProp4_.setVolume(volumeInport_.getData(), 3);

    transFunc1DGaussianProp1_.setVolume(volumeInport_.getData(), 0);
    if (channels >= 2) transFunc1DGaussianProp2_.setVolume(volumeInport_.getData(), 1);
    if (channels >= 3) transFunc1DGaussianProp3_.setVolume(volumeInport_.getData(), 2);
    if (channels >= 4) transFunc1DGaussianProp4_.setVolume(volumeInport_.getData(), 3);
}

} // namespace
