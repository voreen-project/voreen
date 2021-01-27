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

#include "thresholdtotransfunclinker.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"

namespace voreen {

const std::string ThresholdToTransFuncLinker::loggerCat_("voreen.octreeregiongrowing.ThresholdToTransFuncLinker");

ThresholdToTransFuncLinker::ThresholdToTransFuncLinker()
    : Processor()
    //ports
    , volumePort_(Port::INPORT, "volumePort", "Octree Volume")
    //current mode
    , thresholdChannel_("thresholdChannel", "Threshold Visualization")
    //tf
    , channel1TransFuncProp_("channel1TransFuncProp","To be linked with \"SliceViewer.Transfer Function\"")
    , color1Prop_("color1Prop", "Color Map of Channel 1",tgt::vec4(1.f,1.f,1.f,1.f))
    , channel2TransFuncProp_("channel2TransFuncProp","To be linked with \"SliceViewer.Transfer Function 2\"")
    , color2Prop_("color2Prop", "Color Map of Channel 2",tgt::vec4(1.f,0.f,0.f,1.f))
    , channel3TransFuncProp_("channel3TransFuncProp","To be linked with \"SliceViewer.Transfer Function 3\"")
    , color3Prop_("color3Prop", "Color Map of Channel 3",tgt::vec4(0.f,1.f,0.f,1.f))
    , channel4TransFuncProp_("channel4TransFuncProp","To be linked with \"SliceViewer.Transfer Function 4\"")
    , color4Prop_("color4Prop", "Color Map of Channel 4",tgt::vec4(0.f,0.f,1.f,1.f))
    //links, e.g. to octree quantification
    , channelAProp_("channelAProp", "Channel A (to be linked)", 1, 1, 4)
    , thresholdAProp_("thresholdAProp", "Channel A: (absolute) Threshold Region", 1000, 1, INT_MAX)
    , maxValueAProp_("maxValueAProp", "Channel A: Max Value (to be linked)", 1, 1, INT_MAX)
    , channelBProp_("channelBProp", "Channel B (to be linked)", 1, 1, 4)
    , thresholdBProp_("thresholdBProp", "Channel B: (absolute) Threshold Region", 1000, 1, INT_MAX)
    , maxValueBProp_("maxValueBProp", "Channel B: Max Value (to be linked)", 1, 1, INT_MAX)
    , channelCProp_("channelCProp", "Channel C (to be linked)", 1, 1, 4)
    , thresholdCProp_("thresholdCProp", "Channel C: (absolute) Threshold Region", 1000, 1, INT_MAX)
    , maxValueCProp_("maxValueCProp", "Channel C: Max Value (to be linked)", 1, 1, INT_MAX)
{
    //ports
    volumePort_.onChange(MemberFunctionCallback<ThresholdToTransFuncLinker>(this, &ThresholdToTransFuncLinker::adjustToInputVolume));
    addPort(volumePort_);

    //current mode
    thresholdChannel_.addOption("channelA", "Channel A",STTTFL_QUANTIFICATION_A);
    thresholdChannel_.addOption("channelB", "Channel B",STTTFL_QUANTIFICATION_B);
    thresholdChannel_.addOption("channelC", "Channel C",STTTFL_QUANTIFICATION_C);
    addProperty(thresholdChannel_);

    //color values
    addProperty(color1Prop_);
    addProperty(color2Prop_);
    addProperty(color3Prop_);
    addProperty(color4Prop_);

    //tb be linked properties
    addProperty(channel1TransFuncProp_);
    addProperty(channel2TransFuncProp_);
    addProperty(channel3TransFuncProp_);
    addProperty(channel4TransFuncProp_);
    addProperty(channelAProp_);
    addProperty(thresholdAProp_);
    addProperty(maxValueAProp_);
    addProperty(channelBProp_);
    addProperty(thresholdBProp_);
    addProperty(maxValueBProp_);    
    addProperty(channelCProp_);
    addProperty(thresholdCProp_);
    addProperty(maxValueCProp_);        
    
    //hide linkes properties
    channel1TransFuncProp_.setVisibleFlag(false);
    channel2TransFuncProp_.setVisibleFlag(false);
    channel3TransFuncProp_.setVisibleFlag(false);
    channel4TransFuncProp_.setVisibleFlag(false);
    channelAProp_.setVisibleFlag(false);
    thresholdAProp_.setVisibleFlag(false);
    maxValueAProp_.setVisibleFlag(false);
    channelBProp_.setVisibleFlag(false);
    thresholdBProp_.setVisibleFlag(false);
    maxValueBProp_.setVisibleFlag(false);
    channelCProp_.setVisibleFlag(false);
    thresholdCProp_.setVisibleFlag(false);
    maxValueCProp_.setVisibleFlag(false);

    //add on change functions
    color1Prop_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,size_t>(this,
                                        &ThresholdToTransFuncLinker::onChangeCreateOpaqueTransFunc,1));
    color2Prop_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,size_t>(this,
                                        &ThresholdToTransFuncLinker::onChangeCreateOpaqueTransFunc,2));
    color3Prop_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,size_t>(this,
                                        &ThresholdToTransFuncLinker::onChangeCreateOpaqueTransFunc,3));
    color4Prop_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,size_t>(this,
                                        &ThresholdToTransFuncLinker::onChangeCreateOpaqueTransFunc,4));

    thresholdChannel_.onChange(MemberFunctionCallback<ThresholdToTransFuncLinker>(this,&ThresholdToTransFuncLinker::onChangeThresholdChannel));
    channelCProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_C));
    thresholdCProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_C));
    maxValueCProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_C));
    channelAProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_A));
    thresholdAProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_A));
    maxValueAProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_A));
    channelBProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_B));
    thresholdBProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_B));
    maxValueBProp_.onChange(MemberFunctionOneParameterCallback<ThresholdToTransFuncLinker,StttflMode>(this,
                                      &ThresholdToTransFuncLinker::onChangeOtherProperty,STTTFL_QUANTIFICATION_B));

    //tf should always fit to domain
    channel1TransFuncProp_.setDomainFittingStrategy(TransFuncPropertyBase::FIT_DOMAIN_ALWAYS);
    channel2TransFuncProp_.setDomainFittingStrategy(TransFuncPropertyBase::FIT_DOMAIN_ALWAYS);
    channel3TransFuncProp_.setDomainFittingStrategy(TransFuncPropertyBase::FIT_DOMAIN_ALWAYS);
    channel4TransFuncProp_.setDomainFittingStrategy(TransFuncPropertyBase::FIT_DOMAIN_ALWAYS);
}

Processor* ThresholdToTransFuncLinker::create() const {
    return new ThresholdToTransFuncLinker();
}

bool ThresholdToTransFuncLinker::isReady() const {
    return (isInitialized() && volumePort_.hasData());
}

void ThresholdToTransFuncLinker::adjustToInputVolume() {
    //adjust properties
    if(volumePort_.hasData()) {
        color1Prop_.setVisibleFlag(false);
        color2Prop_.setVisibleFlag(false);
        color3Prop_.setVisibleFlag(false);
        color4Prop_.setVisibleFlag(false);
        switch(volumePort_.getData()->getNumChannels()) {
        case 4:
            channel4TransFuncProp_.setVolume(volumePort_.getData(),3);
            color4Prop_.setVisibleFlag(true);
        case 3:
            channel3TransFuncProp_.setVolume(volumePort_.getData(),2);
            color3Prop_.setVisibleFlag(true);
        case 2:
            channel2TransFuncProp_.setVolume(volumePort_.getData(),1);
            color2Prop_.setVisibleFlag(true);
        case 1:
            channel1TransFuncProp_.setVolume(volumePort_.getData(),0);
            color1Prop_.setVisibleFlag(true);
        }
    }
}

void ThresholdToTransFuncLinker::process() {
    // do nothing 
}

void ThresholdToTransFuncLinker::initialize() {
    //init tf properties
    Processor::initialize();
    //set tfs on first start
    onChangeCreateOpaqueTransFunc(1);
    onChangeCreateOpaqueTransFunc(2);
    onChangeCreateOpaqueTransFunc(3);
    onChangeCreateOpaqueTransFunc(4);
    //adjust visibilities
    onChangeThresholdChannel();
}

//-------------------------------------------------------------------
//  onChange functions
//-------------------------------------------------------------------
void ThresholdToTransFuncLinker::onChangeCreateOpaqueTransFunc(size_t channel) {
    tgtAssert(channel >= 1 && channel <= 4, "channel must be between 1-4");

    //get properties associated with channel
    TransFunc1DKeysProperty* tmpTF = 0; ColorProperty* tmpCol = 0;
    switch(channel) {
    case 1:
        tmpTF = &channel1TransFuncProp_;
        tmpCol = &color1Prop_;
        break;
    case 2:
        tmpTF = &channel2TransFuncProp_;
        tmpCol = &color2Prop_;
        break;
    case 3:
        tmpTF = &channel3TransFuncProp_;
        tmpCol = &color3Prop_;
        break;
    case 4:
        tmpTF = &channel4TransFuncProp_;
        tmpCol = &color4Prop_;
        break;
    default:
        tgtAssert(false,"wrong channel");
    }

    //create opaque keys
    std::vector<TransFuncMappingKey*> vec;
    vec.push_back(new TransFuncMappingKey(0.f, tmpCol->getRGBA()));
    vec.push_back(new TransFuncMappingKey(1.f, tmpCol->getRGBA()));

    //set tf
    TransFunc1DKeys* tf = dynamic_cast<TransFunc1DKeys*>(tmpTF->get());
    tgtAssert(tf, "no tf");
    //tf gets ownership
    tf->setKeys(vec);
    tmpTF->invalidate();
}

void ThresholdToTransFuncLinker::onChangeOtherProperty(StttflMode associatedMode) {
    //check, if value change is relevant atm
    if(thresholdChannel_.getValue() == associatedMode)
        onChangeThresholdChannel();
}

void ThresholdToTransFuncLinker::onChangeThresholdChannel() {
    //get current values
    int channel = 0; /*int absValue = 0;*/ tgt::ivec2 absThreshold(static_cast<int>(0)); int maxValue = 0;
    switch(thresholdChannel_.getValue()) {
    case STTTFL_QUANTIFICATION_A:
        channel = channelAProp_.get();
        absThreshold = thresholdAProp_.get();
        maxValue = maxValueAProp_.get();
        break;
    case STTTFL_QUANTIFICATION_B:
        channel = channelBProp_.get();
        absThreshold = thresholdBProp_.get();
        maxValue = maxValueBProp_.get();
        break;
    case STTTFL_QUANTIFICATION_C:
        channel = channelCProp_.get();
        absThreshold = thresholdCProp_.get();
        maxValue = maxValueCProp_.get();
        break;
    default:
        tgtAssert(false,"unknown stttf");
    }

    //help vec (element 0 is empty to use indices between 1-4)
    std::vector<TransFunc1DKeysProperty*> vec(5,static_cast<TransFunc1DKeysProperty*>(0));
    vec[1] = &channel1TransFuncProp_; vec[2] = &channel2TransFuncProp_;
    vec[3] = &channel3TransFuncProp_; vec[4] = &channel4TransFuncProp_;
    //set threshold
    for(int i = 1; i <= 4; i++) {
        if(i == channel) {
            vec[i]->get()->setAlphaMode(TransFuncBase::TF_ONE_ALPHA);
            if(TransFunc1DKeys* tf = dynamic_cast<TransFunc1DKeys*>(vec[i]->get())) {
                //LERROR("rwm: " <<tf->normalizedToRealWorld(static_cast<float>(absValue)/static_cast<float>(maxValue)));
                //LERROR("normal: " << static_cast<float>(absValue)/static_cast<float>(maxValue));
                tf->setThreshold(static_cast<float>(absThreshold.x)/static_cast<float>(maxValue), 
                                 std::min(1.f, static_cast<float>(absThreshold.y) / static_cast<float>(maxValue)));
            } else {
                tgtAssert(false, "no 1d keys tf");
            }
        } else {
            //set alpha zero for non active transfunc
            vec[i]->get()->setAlphaMode(TransFuncBase::TF_ZERO_ALPHA);
        }
        //invalidate property to execute links
        vec[i]->invalidate();
    }
}


} // namespace voreen
