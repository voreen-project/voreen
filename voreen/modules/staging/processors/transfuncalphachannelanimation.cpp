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

#include "transfuncalphachannelanimation.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"

namespace voreen {

const std::string TransFuncAlphaChannelAnimation::loggerCat_("voreen.staging.TransFuncAlphaChannelAnimation");

TransFuncAlphaChannelAnimation::TransFuncAlphaChannelAnimation()
    : Processor()
    // tf
    , channel1TransFuncProp_("channel1TransFuncProp","Transfer Function 1 (for linking)")
    , channel2TransFuncProp_("channel2TransFuncProp","Transfer Function 2 (for linking)")
    , channel3TransFuncProp_("channel3TransFuncProp","Transfer Function 3 (for linking)")
    , channel4TransFuncProp_("channel4TransFuncProp","Transfer Function 4 (for linking)")
    // options
    , channel1AlphaMode_("channel1alphamode", "Channel 1 Alpha Mode") 
    , channel2AlphaMode_("channel2alphamode", "Channel 2 Alpha Mode")
    , channel3AlphaMode_("channel3alphamode", "Channel 3 Alpha Mode")
    , channel4AlphaMode_("channel4alphamode", "Channel 4 Alpha Mode")
{
    // add options and option properties
    channel1AlphaMode_.addOption("alpha", "Use Alpha", 0);
    channel1AlphaMode_.addOption("transparent", "Transparent", 1);
    channel1AlphaMode_.addOption("opaque", "Opaque", 2);
    addProperty(channel1AlphaMode_);

    channel2AlphaMode_.addOption("alpha", "Use Alpha", 0);
    channel2AlphaMode_.addOption("transparent", "Transparent", 1);
    channel2AlphaMode_.addOption("opaque", "Opaque", 2);
    addProperty(channel2AlphaMode_);

    channel3AlphaMode_.addOption("alpha", "Use Alpha", 0);
    channel3AlphaMode_.addOption("transparent", "Transparent", 1);
    channel3AlphaMode_.addOption("opaque", "Opaque", 2);
    addProperty(channel3AlphaMode_);

    channel4AlphaMode_.addOption("alpha", "Use Alpha", 0);
    channel4AlphaMode_.addOption("transparent", "Transparent", 1);
    channel4AlphaMode_.addOption("opaque", "Opaque", 2);
    addProperty(channel4AlphaMode_);

    // tf linking properties
    addProperty(channel1TransFuncProp_);
    addProperty(channel2TransFuncProp_);
    addProperty(channel3TransFuncProp_);
    addProperty(channel4TransFuncProp_);        

    // register callbacks
    channel1TransFuncProp_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onTransFuncChange1));
    channel2TransFuncProp_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onTransFuncChange2));
    channel3TransFuncProp_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onTransFuncChange3));
    channel4TransFuncProp_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onTransFuncChange4));

    channel1AlphaMode_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onModeChange1));
    channel2AlphaMode_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onModeChange2));
    channel3AlphaMode_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onModeChange3));
    channel4AlphaMode_.onChange(MemberFunctionCallback<TransFuncAlphaChannelAnimation>(this, &TransFuncAlphaChannelAnimation::onModeChange4));
}

Processor* TransFuncAlphaChannelAnimation::create() const {
    return new TransFuncAlphaChannelAnimation();
}

bool TransFuncAlphaChannelAnimation::isReady() const {
    return isInitialized();
}

TransFuncBase::AlphaMode TransFuncAlphaChannelAnimation::intToAlphaMode(int i) {
    switch (i) {
        case 0 : { return TransFuncBase::TF_USE_ALPHA;  }
        case 1 : { return TransFuncBase::TF_ZERO_ALPHA; }
        case 2 : { return TransFuncBase::TF_ONE_ALPHA;  }
        default: { tgtAssert(false, "Should not get here"); 
                   return TransFuncBase::TF_USE_ALPHA; // fixes compiler warning 
                 }
    }
}

int TransFuncAlphaChannelAnimation::alphaModeToInt(TransFuncBase::AlphaMode mode) {
    switch (mode) {
        case TransFuncBase::TF_USE_ALPHA :  {   return 0; }
        case TransFuncBase::TF_ZERO_ALPHA : {   return 1; }
        case TransFuncBase::TF_ONE_ALPHA :  {   return 2; }
        default: { tgtAssert(false, "Should not get here"); 
                   return TransFuncBase::TF_USE_ALPHA; // fixes compiler warning 
                 }
    }
}

//-------------------------------------------------------------------
//  onChange functions
//-------------------------------------------------------------------

void TransFuncAlphaChannelAnimation::onTransFuncChange1() {
    if (channel1TransFuncProp_.get())
        channel1AlphaMode_.selectByValue(alphaModeToInt(channel1TransFuncProp_.get()->getAlphaMode()));
}

void TransFuncAlphaChannelAnimation::onTransFuncChange2() {
    if (channel2TransFuncProp_.get())
        channel2AlphaMode_.selectByValue(alphaModeToInt(channel2TransFuncProp_.get()->getAlphaMode()));
}

void TransFuncAlphaChannelAnimation::onTransFuncChange3() {
    if (channel3TransFuncProp_.get())
        channel3AlphaMode_.selectByValue(alphaModeToInt(channel3TransFuncProp_.get()->getAlphaMode()));
}

void TransFuncAlphaChannelAnimation::onTransFuncChange4() {
    if (channel4TransFuncProp_.get())
        channel4AlphaMode_.selectByValue(alphaModeToInt(channel4TransFuncProp_.get()->getAlphaMode()));
}


void TransFuncAlphaChannelAnimation::onModeChange1() {
    if (channel1TransFuncProp_.get()) {
        channel1TransFuncProp_.get()->setAlphaMode(intToAlphaMode(channel1AlphaMode_.getValue()));
        channel1TransFuncProp_.invalidate();
    }
}

void TransFuncAlphaChannelAnimation::onModeChange2() {
    if (channel2TransFuncProp_.get()) {
        channel2TransFuncProp_.get()->setAlphaMode(intToAlphaMode(channel2AlphaMode_.getValue()));
        channel2TransFuncProp_.invalidate();
    }
}

void TransFuncAlphaChannelAnimation::onModeChange3() {
    if (channel3TransFuncProp_.get()) {
        channel3TransFuncProp_.get()->setAlphaMode(intToAlphaMode(channel3AlphaMode_.getValue()));
        channel3TransFuncProp_.invalidate();
    }
}

void TransFuncAlphaChannelAnimation::onModeChange4() {
    if (channel4TransFuncProp_.get()) {
        channel4TransFuncProp_.get()->setAlphaMode(intToAlphaMode(channel4AlphaMode_.getValue()));
        channel4TransFuncProp_.invalidate();
    }
}

} // namespace voreen
