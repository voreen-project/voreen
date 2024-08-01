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

#include "transfuncchannellinker.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"

namespace voreen {

const std::string TransFuncChannelLinker::loggerCat_("voreen.octreequantification.TransFuncChannelLinker");

TransFuncChannelLinker::TransFuncChannelLinker()
    : Processor()
    //tf
    , channel1TransFuncProp_("channel1TransFuncProp","First input transfer function")
    , channel2TransFuncProp_("channel2TransFuncProp","Second input transfer function")
    , channel3TransFuncProp_("channel3TransFuncProp","Third input transfer function")
    , channel4TransFuncProp_("channel4TransFuncProp","Fourth input transfer function")
    //links, e.g. to octree quantification
    , channelProp_("channel", "Selected output channel", 1, 1, 4)
    , outputTransFuncProp_("outputTransFuncProp","Output transfer function")
{
    addProperty(channel1TransFuncProp_);
    addProperty(channel2TransFuncProp_);
    addProperty(channel3TransFuncProp_);
    addProperty(channel4TransFuncProp_);
    addProperty(channelProp_);

    addProperty(outputTransFuncProp_);
    outputTransFuncProp_.setReadOnlyFlag(true);

    ON_CHANGE(channel1TransFuncProp_, TransFuncChannelLinker, onChangeCallback);
    ON_CHANGE(channel2TransFuncProp_, TransFuncChannelLinker, onChangeCallback);
    ON_CHANGE(channel3TransFuncProp_, TransFuncChannelLinker, onChangeCallback);
    ON_CHANGE(channel4TransFuncProp_, TransFuncChannelLinker, onChangeCallback);
    ON_CHANGE(channelProp_, TransFuncChannelLinker, onChangeCallback);
}

Processor* TransFuncChannelLinker::create() const {
    return new TransFuncChannelLinker();
}


void TransFuncChannelLinker::process() {
    // do nothing 
}


void TransFuncChannelLinker::onChangeCallback() {
    if (!isInitialized())
        return;

    if (channelProp_.get() == 1) {
        outputTransFuncProp_.get()->setMemberValuesFrom(channel1TransFuncProp_.get());
        outputTransFuncProp_.setDomainFittingStrategy(channel1TransFuncProp_.getDomainFittingStrategy());
    }
    else if (channelProp_.get() == 2) {
        outputTransFuncProp_.get()->setMemberValuesFrom(channel2TransFuncProp_.get());
        outputTransFuncProp_.setDomainFittingStrategy(channel2TransFuncProp_.getDomainFittingStrategy());
    }
    else if (channelProp_.get() == 3) {
        outputTransFuncProp_.get()->setMemberValuesFrom(channel3TransFuncProp_.get());
        outputTransFuncProp_.setDomainFittingStrategy(channel3TransFuncProp_.getDomainFittingStrategy());
    }
    else if (channelProp_.get() == 4) {
        outputTransFuncProp_.get()->setMemberValuesFrom(channel4TransFuncProp_.get());
        outputTransFuncProp_.setDomainFittingStrategy(channel4TransFuncProp_.getDomainFittingStrategy());
    }

    outputTransFuncProp_.invalidate();
}

} // namespace voreen
