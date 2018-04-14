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

#include "segmentationraycaster.h"

#include "voreen/core/datastructures/volume/modality.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "tgt/gpucapabilities.h"
#include "tgt/textureunit.h"

#include <sstream>
#include <cstdio>
#include <algorithm>

namespace voreen {

    /** Currently supported segments. */
    const size_t     MAX_SEGMENTATIONS = 256;
    /** Currently used GUI color for defined segments. */
    const tgt::col4 HIGHLIGHT_COLOR   = tgt::col4(255,0,0,255);

SegmentationRaycaster::SegmentationRaycaster()
    : VolumeRaycaster(true, true) //add default properties
    //ports
    , segmentationInport_(Port::INPORT, "volumehandle.segmentation", "Segmentation Volume Input")
    //properties
        //rendering
    , defaultTransFuncProp_("SegmentationRaycaster.TransFunc", "Default Transfer Function")
     //segmentation
    , applySegmentationProp_("SegmentationRaycatser.applySegmentation", "Apply Segmentation", false, Processor::INVALID_PROGRAM)
    , segmentIndexProp_("segmentIndexProp","Current Segment",Processor::INVALID_RESULT,true)
    , isDefaultTFProp_("isDefaultTFProp","Default Function", true)
    , resetCurrentSegmentProp_("resetCurrentSegmentProp", "Reset Current Function")
    , segmentTransFuncProp_("segmentTransFuncProp","Segment Transfer Function")
    , resetAllSegmentsProp_("resetAllSegmentsProp", "Reset All Functions")
        //debug
    , shaderProp_("raycast.prg", "Raycasting Shader", "rc_segmentation.frag", "passthrough.vert","",Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , cameraProp_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)), true)
    //handle segmentation
    , segmentTransFuncVector_(MAX_SEGMENTATIONS)
    , segmentationTransFuncTex_(0)
    , segmentationTransFuncTexValid_(false)
{
    // initialize segment tf vector with zeros     // not needed default ist zero
    //TransFunc1DKeys* segmentTransFuncArray_[MAX_SEGMENTATIONS];
    //std::memset(segmentTransFuncArray_, 0, sizeof(segmentTransFuncArray_));
    //segmentTransFuncVector_.assign(&segmentTransFuncArray_[0],&segmentTransFuncArray_[0]+MAX_SEGMENTATIONS);

    // ports
    volumeInport_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::volumeInportOnChange));
    PortConditionLogicalOr* groupCondition = new PortConditionLogicalOr();
    groupCondition->addLinkedCondition(new PortConditionVolumeTypeUInt8());
    groupCondition->addLinkedCondition(new PortConditionVolumeTypeUInt12());
    groupCondition->addLinkedCondition(new PortConditionVolumeTypeUInt16());
    segmentationInport_.addCondition(groupCondition);
    //segmentationInport_.showTextureAccessProperties(true); // is always nearest filtering... hard-coded in process()
    addPort(segmentationInport_);



    //properties
        //default (raycasting properties added by volumeraycaster)
        //tf and compositing
    addProperty(defaultTransFuncProp_);
        defaultTransFuncProp_.setGroupID("tf");
    addProperty(isoValue_); //VolumeRaycaster
        isoValue_.setGroupID("tf");
    addProperty(gammaValue1_);
        gammaValue1_.setGroupID("tf");
    addProperty(gammaValue2_);
        gammaValue2_.setGroupID("tf");
     addProperty(gammaValue3_);
        gammaValue3_.setGroupID("tf");
    addProperty(shadeMode_);    //VolumeRaycaster
        shadeMode_.setGroupID("tf");
    setPropertyGroupGuiName("tf", "Rendering Settings");
        //segment
    addProperty(applySegmentationProp_);
    applySegmentationProp_.setGroupID("segment");
    addProperty(segmentIndexProp_);
        for(int i = 0; i < (int) MAX_SEGMENTATIONS; i++)
            segmentIndexProp_.addOption(std::to_string(i),"Segment " + std::to_string(i), i);
        segmentIndexProp_.setGroupID("segment");
    addProperty(isDefaultTFProp_);
        isDefaultTFProp_.setReadOnlyFlag(true);
        isDefaultTFProp_.setGroupID("segment");
    addProperty(resetCurrentSegmentProp_);
        resetCurrentSegmentProp_.setGroupID("segment");
    addProperty(segmentTransFuncProp_);
        segmentTransFuncProp_.setGroupID("segment");
    addProperty(resetAllSegmentsProp_);
        resetAllSegmentsProp_.setGroupID("segment");
    setPropertyGroupGuiName("segment", "Segment Settings");
        //light settings
    addProperty(gradientMode_); //VolumeRaycaster
        gradientMode_.setGroupID("light");
    addProperty(lightPosition_); //VolumeRenderer
        lightPosition_.setGroupID("light");
    addProperty(lightAmbient_);  //VolumeRenderer
        lightAmbient_.setGroupID("light");
    addProperty(lightDiffuse_);  //VolumeRenderer
        lightDiffuse_.setGroupID("light");
    addProperty(lightSpecular_); //VolumeRenderer
        lightSpecular_.setGroupID("light");
    setPropertyGroupGuiName("light", "Light Settings");
        //debug and unimportant
    addProperty(shaderProp_);
    addProperty(cameraProp_);

    //callbacks
    defaultTransFuncProp_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::defaultTransFuncOnChange));
    applySegmentationProp_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::applySegmentationOnChange));
    segmentIndexProp_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::segmentIndexOnChange));
    resetCurrentSegmentProp_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::resetCurrentSegment));
    segmentTransFuncProp_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::segmentTransFuncOnChange));
    resetAllSegmentsProp_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::resetAllSegments));

    // listen to changes of properties that influence the GUI state (i.e. visibility of other props)
    shadeMode_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::adjustPropertyVisibilities));
    compositingMode1_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::adjustPropertyVisibilities));
    compositingMode2_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::adjustPropertyVisibilities));
    compositingMode3_.onChange(MemberFunctionCallback<SegmentationRaycaster>(this, &SegmentationRaycaster::adjustPropertyVisibilities));

    std::fill(segmentTransFuncVector_.begin(), segmentTransFuncVector_.end(), nullptr);
}

SegmentationRaycaster::~SegmentationRaycaster() {
}

bool SegmentationRaycaster::isReady() const {
    //check if all inports are connected:
    if(!entryPort_.isReady() || !exitPort_.isReady() || !volumeInport_.isReady())
        return false;

    if (!outport1_.isReady() && !outport2_.isReady() && !outport3_.isReady())
        return false;

    if (applySegmentationProp_.get() && !segmentationInport_.isReady())
        return false;

    return true;
}

void SegmentationRaycaster::initialize() {

    VolumeRaycaster::initialize();
    compile();

    //create/initialize segmentation Function
    segmentationTransFuncTex_ = new tgt::Texture(tgt::ivec3(TRANSFUNC_TEXTURE_SIZE, MAX_SEGMENTATIONS * 3, 1), GL_RGBA,
                                                 GL_RGBA, GL_UNSIGNED_BYTE, tgt::Texture::LINEAR, tgt::Texture::CLAMP);
    segmentationTransFuncTex_->alloc(true);
    updateEntireSegmentationTransFuncTexture();

    //adjust visibilities
    applySegmentationOnChange();
    adjustPropertyVisibilities();
    //HACK: if this processor replaces another (drag/drop) the volume from the inport is assigned to the
    //      TFs before it is initialized. This results in trouble. The solution is to call the callback after the initialization.
    volumeInportOnChange();
    //update flags
    segmentIndexOnChange();
    ON_CHANGE(applySegmentationProp_, SegmentationRaycaster, updateEntireSegmentationTransFuncTexture);
}

void SegmentationRaycaster::deinitialize() {
    //delete and clear the combined segmentation texture
    delete segmentationTransFuncTex_;
    segmentationTransFuncTex_ = 0;
    LGL_ERROR;


    VolumeRaycaster::deinitialize();

    //clear transfer functions
    for(auto it = segmentTransFuncVector_.begin(); it != segmentTransFuncVector_.end(); it++)
        delete *it;
}

void SegmentationRaycaster::compile() {
    shaderProp_.setHeader(generateHeader());
    shaderProp_.rebuild();
}

std::string SegmentationRaycaster::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string headerSource = VolumeRaycaster::generateHeader();

    if (defaultTransFuncProp_.get())
        headerSource += defaultTransFuncProp_.get()->getShaderDefines();

    if (applySegmentationProp_.get() && segmentationInport_.getData()) {
        headerSource += "#define MOD_APPLY_SEGMENTATION\n";
        headerSource += "#define SEGMENTATION_TRANSFUNC_HEIGHT " + std::to_string(MAX_SEGMENTATIONS*3) + "\n";
    }



    internalPortGroup_.reattachTargets();
    headerSource += internalPortGroup_.generateHeader(shaderProp_.getShader());
    return headerSource;
}

void SegmentationRaycaster::beforeProcess() {
    VolumeRaycaster::beforeProcess();

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;
}

void SegmentationRaycaster::process() {
    compile();
    const VolumeBase* volumeHandle = volumeInport_.getData();

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
    tgt::Camera cam = cameraProp_.get();
    setGlobalShaderParameters(raycastPrg, &cam, renderSize);
    LGL_ERROR;

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

    // vector containing the volumes to bind; is passed to bindVolumes()
    std::vector<VolumeStruct> volumeTextures;

    // add main volume
    tgt::TextureUnit volTexUnit;
    volumeTextures.push_back(VolumeStruct(
        volumeHandle,
        &volTexUnit,
        "volume_","volumeStruct_")
    );

    // segmentation volume
    const VolumeGL* segVolume = 0;
    if (segmentationInport_.getData()) {
        segVolume = segmentationInport_.getData()->getRepresentation<VolumeGL>();
    }

    tgt::TextureUnit segUnit, segTransferUnit, transferUnit;
    if (segVolume) {
        volumeTextures.push_back(VolumeStruct(
            segmentationInport_.getData(),
            &segUnit,
            "segmentation_","segmentationParameters_",
            GL_CLAMP_TO_EDGE,
            tgt::vec4(0.f),
            GL_NEAREST)
        );
    }

    if (segVolume && applySegmentationProp_.get()) {
        segTransferUnit.activate();
        segmentationTransFuncTex_->uploadTexture();
        segmentationTransFuncTex_->bind();
        if (!segmentationTransFuncTexValid_) {
            segmentationTransFuncTex_->uploadTexture();
            segmentationTransFuncTexValid_ = true;
        }
    }
    else {
        transferUnit.activate();
        defaultTransFuncProp_.get()->getTexture()->bind();
    }
    // bind the volumes and pass the necessary information to the shader
    bindVolumes(raycastPrg, volumeTextures, &cam, lightPosition_.get());

    // assign segment tf to texture units
    if (segVolume && applySegmentationProp_.get()) {
        raycastPrg->setUniform("segmentationTransferFunc_", segTransferUnit.getUnitNumber());
        raycastPrg->setUniform("segmentationTransferFuncDomains_",&segmentationTransferFuncDomains_[0],MAX_SEGMENTATIONS); //HACK: max segemntations must be adapted in the fragment shader
    } else {
        defaultTransFuncProp_.get()->setUniform(raycastPrg, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());
    }

    // pass remaining uniforms to shader
    if (compositingMode1_.isSelected("iso")  ||
        compositingMode2_.isSelected("iso") ||
        compositingMode3_.isSelected("iso") )
        raycastPrg->setUniform("isoValue_", isoValue_.get());

    if (compositingMode1_.isSelected("mida"))
        raycastPrg->setUniform("gammaValue_", gammaValue1_.get());

    if (compositingMode2_.isSelected("mida"))
        raycastPrg->setUniform("gammaValue1_", gammaValue2_.get());

    if (compositingMode3_.isSelected("mida"))
        raycastPrg->setUniform("gammaValue2_", gammaValue3_.get());


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
    if (outport3_.isConnected())
        rescaleRendering(internalRenderPort3_, outport3_);

    tgt::TextureUnit::setZeroUnit();

    LGL_ERROR;
}

//-------------------------------------------------------------------//
//  Handle Segmentations                                             //
//-------------------------------------------------------------------//
void SegmentationRaycaster::serialize(Serializer& s) const {
    VolumeRaycaster::serialize(s);
    s.serialize("segmentTransFuncVector",segmentTransFuncVector_);
}

void SegmentationRaycaster::deserialize(Deserializer& d) {
    VolumeRaycaster::deserialize(d);
    d.deserialize("segmentTransFuncVector",segmentTransFuncVector_);
}

void SegmentationRaycaster::updateEntireSegmentationTransFuncTexture() {
    for (int i=0; i < (int) MAX_SEGMENTATIONS; ++i)
        updateSegmentationTransFuncTexture(i);
}

void SegmentationRaycaster::updateSegmentationTransFuncTexture(int segment) {

    TransFunc1DKeys* intensityTF = nullptr;

    if(segmentTransFuncVector_[segment] == nullptr) { //use default function
        intensityTF = defaultTransFuncProp_.get();
    } else {
        intensityTF = segmentTransFuncVector_[segment];
    }

    tgtAssert(intensityTF, "1D transfer function expected");
    tgtAssert(segmentationTransFuncTex_, "No segmentation transfer function texture");
    tgtAssert(segment >= 0 && segment < tgt::ifloor(segmentationTransFuncTex_->getHeight() / 3.f), "Invalid segment id");

    // A segment's 1D transfer function is stored in the 2D segmentation tf texture as a
    // 3-row wide stripe which is centered around the row 3*i+1.
    int tfRow = 3*segment+1;
    size_t bytesPerRow = segmentationTransFuncTex_->getDimensions().x * segmentationTransFuncTex_->getBpt();

    uint32_t *out = reinterpret_cast<uint32_t*>(&segmentationTransFuncTex_->getCpuTextureData()[(tfRow-1)*bytesPerRow]);
    int outsize = segmentationTransFuncTex_->getDimensions().x;
    uint32_t *in = reinterpret_cast<uint32_t*>(intensityTF->getTexture()->getCpuTextureData());
    int insize = intensityTF->getDimensions().x;

    resampleTexture1D(out, outsize, in, insize);
    out = reinterpret_cast<uint32_t*>(&segmentationTransFuncTex_->getCpuTextureData()[(tfRow)*bytesPerRow]);
    resampleTexture1D(out, outsize, in, insize);
    out = reinterpret_cast<uint32_t*>(&segmentationTransFuncTex_->getCpuTextureData()[(tfRow+1)*bytesPerRow]);
    resampleTexture1D(out, outsize, in, insize);

    segmentationTransferFuncDomains_[segment] = intensityTF->getDomain();

    // transfer func texture has to be uploaded in the next rendering pass
    segmentationTransFuncTexValid_ = false;
}

//-------------------------------------------------------------------//
//  Callbacks                                                        //
//-------------------------------------------------------------------//
void SegmentationRaycaster::adjustPropertyVisibilities() {
    bool useLighting = !shadeMode_.isSelected("none");
    setPropertyGroupVisible("light", useLighting);

    bool useIsovalue = (compositingMode1_.isSelected("iso")  ||
                        compositingMode2_.isSelected("iso") ||
                        compositingMode3_.isSelected("iso")   );
    isoValue_.setVisibleFlag(useIsovalue);

    gammaValue1_.setVisibleFlag(compositingMode1_.isSelected("mida"));
    gammaValue2_.setVisibleFlag(compositingMode2_.isSelected("mida"));
    gammaValue3_.setVisibleFlag(compositingMode3_.isSelected("mida"));
}

void SegmentationRaycaster::volumeInportOnChange() {
    //adjust camera to scene
    if(volumeInport_.hasData())
        cameraProp_.adaptInteractionToScene(volumeInport_.getData()->getBoundingBox().getBoundingBox(), tgt::min(volumeInport_.getData()->getSpacing()));

    //update tf properties
    defaultTransFuncProp_.setVolume(volumeInport_.getData());
    segmentTransFuncProp_.setVolume(volumeInport_.getData());
}

void SegmentationRaycaster::defaultTransFuncOnChange() {
    //ignore on initialize
    if(!isInitialized()) return;
    //update segment function
    if(isDefaultTFProp_.get()) {
        segmentTransFuncProp_.blockCallbacks(true);
        segmentTransFuncProp_.get()->setMemberValuesFrom(defaultTransFuncProp_.get());
        segmentTransFuncProp_.updateWidgets();
        segmentTransFuncProp_.blockCallbacks(false);
    }
    //texture must be updated
    if (segmentationTransFuncTex_)
        updateEntireSegmentationTransFuncTexture();
}

void SegmentationRaycaster::applySegmentationOnChange() {
    segmentIndexProp_.setVisibleFlag(applySegmentationProp_.get());
    isDefaultTFProp_.setVisibleFlag(applySegmentationProp_.get());
    resetCurrentSegmentProp_.setVisibleFlag(applySegmentationProp_.get());
    segmentTransFuncProp_.setVisibleFlag(applySegmentationProp_.get());
    resetAllSegmentsProp_.setVisibleFlag(applySegmentationProp_.get());
}

void SegmentationRaycaster::segmentIndexOnChange() {
    //set default tf if null
    if(segmentTransFuncVector_.at(segmentIndexProp_.getValue()) == nullptr) {
        segmentTransFuncProp_.blockCallbacks(true);
        segmentTransFuncProp_.set1DKeys(defaultTransFuncProp_.get()->clone());
        segmentTransFuncProp_.blockCallbacks(false);
        isDefaultTFProp_.set(true);
    } else { //else clone function
        segmentTransFuncProp_.blockCallbacks(true);
        segmentTransFuncProp_.set1DKeys(segmentTransFuncVector_.at(segmentIndexProp_.getValue())->clone());
        segmentTransFuncProp_.blockCallbacks(false);
        isDefaultTFProp_.set(false);
    }
}

void SegmentationRaycaster::resetCurrentSegment() {
    //set default
    isDefaultTFProp_.set(true);
    segmentTransFuncProp_.blockCallbacks(true);
    segmentTransFuncProp_.set1DKeys(defaultTransFuncProp_.get()->clone());
    segmentTransFuncProp_.blockCallbacks(false);
    //update option color
    const_cast<voreen::Option<int>*>(&segmentIndexProp_.getOptions()[segmentIndexProp_.getValue()])->guiColor_ = tgt::col4::zero;
    segmentIndexProp_.updateWidgets(); //force color update
    //clear old
    delete segmentTransFuncVector_[segmentIndexProp_.getValue()];
    segmentTransFuncVector_[segmentIndexProp_.getValue()] = nullptr;
    //update segmentation TF
    if (segmentationTransFuncTex_)
        updateSegmentationTransFuncTexture(segmentIndexProp_.getValue());
}

void SegmentationRaycaster::segmentTransFuncOnChange() {
    //ignore on initialize
    if(!isInitialized()) return;
    //function has changed
    isDefaultTFProp_.set(false);
    const_cast<voreen::Option<int>*>(&segmentIndexProp_.getOptions()[segmentIndexProp_.getValue()])->guiColor_ = HIGHLIGHT_COLOR;
    segmentIndexProp_.updateWidgets(); //force color update
    //set function in vector
    if(segmentTransFuncVector_.at(segmentIndexProp_.getValue()) == nullptr) {
        segmentTransFuncVector_[segmentIndexProp_.getValue()] = segmentTransFuncProp_.get()->clone();
    } else {
        segmentTransFuncVector_[segmentIndexProp_.getValue()]->setMemberValuesFrom(segmentTransFuncProp_.get());
    }
    //update segmentation TF
    if (segmentationTransFuncTex_)
        updateSegmentationTransFuncTexture(segmentIndexProp_.getValue());
}

void SegmentationRaycaster::resetAllSegments() {
    //set default
    isDefaultTFProp_.set(true);
    segmentTransFuncProp_.blockCallbacks(true);
    segmentTransFuncProp_.set1DKeys(defaultTransFuncProp_.get()->clone());
    segmentTransFuncProp_.blockCallbacks(false);
    //clear old
    for(size_t i = 0; i < MAX_SEGMENTATIONS; i++) {
        delete segmentTransFuncVector_[i];
        segmentTransFuncVector_[i] = nullptr;
        const_cast<voreen::Option<int>*>(&segmentIndexProp_.getOptions()[i])->guiColor_ = tgt::col4::zero;
    }
    segmentIndexProp_.updateWidgets(); //force color update
    //update segmentation TF
    if (segmentationTransFuncTex_)
        updateEntireSegmentationTransFuncTexture();
}

void SegmentationRaycaster::resampleTexture1D( uint32_t *out, int outsize, const uint32_t *in, int insize ){
    float fac = 1.0f*insize/outsize;
    for(int i = 0; i != outsize; i++){
        int idx = static_cast<int>(fac*i);
        if (idx >= insize) idx = insize-1;
        out[i] = in[idx];
    }

}

} // namespace voreen
