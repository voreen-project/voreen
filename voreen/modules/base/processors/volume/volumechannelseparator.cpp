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

#include "volumechannelseparator.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

const std::string VolumeChannelSeparator::loggerCat_("voreen.VolumeChannelSeparator");

VolumeChannelSeparator::VolumeChannelSeparator()
    : CachingVolumeProcessor()
    , volumeInport_(Port::INPORT, "volumeInput", "Volume Input")
    , volumeOutport_(Port::OUTPORT, "volumeOutport", "Volume Output")
    , volumeOutport2_(Port::OUTPORT, "volumeInput2", "Volume Output 2")
    , volumeOutport3_(Port::OUTPORT, "volumeInput3", "Volume Output 3")
    , volumeOutport4_(Port::OUTPORT, "volumeInput4", "Volume Output 4")
    , autoRefresh_("autorefresh", "Auto Refresh", false)
    , refresh_("refreshbutton", "Refresh\n(clear and re-compute)")
    , refreshNecessary_(false)
{
    addPort(volumeInport_);
    addPort(volumeOutport_);
    addPort(volumeOutport2_);
    addPort(volumeOutport3_);
    addPort(volumeOutport4_);

    addProperty(autoRefresh_);
    addProperty(refresh_);
    refresh_.onClick(MemberFunctionCallback<VolumeChannelSeparator>(this, &VolumeChannelSeparator::refresh));
}

Processor* VolumeChannelSeparator::create() const {
    return new VolumeChannelSeparator();
}

void VolumeChannelSeparator::initialize() {
    CachingVolumeProcessor::initialize();
}

void VolumeChannelSeparator::deinitialize() {
    volumeOutport_.clear();
    volumeOutport2_.clear();
    volumeOutport3_.clear();
    volumeOutport4_.clear();
    CachingVolumeProcessor::deinitialize();
}

bool VolumeChannelSeparator::isReady() const {
    return volumeInport_.isReady() && (volumeOutport_.isReady() || volumeOutport2_.isReady() || volumeOutport3_.isReady() || volumeOutport4_.isReady());
}

void VolumeChannelSeparator::beforeProcess() {
    if (autoRefresh_.get())
        refreshNecessary_ = true;

    // if refreshing: look for results in cache
    if (refreshNecessary_)
        CachingVolumeProcessor::beforeProcess();
}

void VolumeChannelSeparator::process() {

    setProgress(0.f);

    volumeOutport_.clear();
    volumeOutport2_.clear();
    volumeOutport3_.clear();
    volumeOutport4_.clear();

    if (!volumeInport_.getData() || !refreshNecessary_) {
        //setProgress(0.f);
        return;
    }

    if (volumeInport_.getData()->getNumChannels() == 1) {
        LINFO("Input volume is single-channel volume... passing through to first outport.");
        volumeOutport_.setData(volumeInport_.getData(), false);
        setProgress(1.f);
    }
    else if (volumeInport_.getData()->getNumChannels() > 4) {
        LERROR("Input volume with more than 4 channels currently not supported. No output generated.");
        setProgress(0.f);
    }
    else {
        //2-4 channels
        std::vector<Volume*> result = separateChannels(volumeInport_.getData());

        if (result.empty()) {
            LERROR("Could not create output volumes");
            return;
        }

        LINFO("Created " << result.size() << " single-channel volumes.");

        //set volumes with outport as owner so that they are correctly destroyed
        volumeOutport_.setData(result.at(0));

        if (result.size() > 1)
            volumeOutport2_.setData(result.at(1));

        if (result.size() > 2)
            volumeOutport3_.setData(result.at(2));

        if (result.size() > 3)
            volumeOutport4_.setData(result.at(3));
    }
}

void VolumeChannelSeparator::afterProcess() {
    VolumeProcessor::afterProcess();

    // do not store the configuration if nothing happened to avoid problems with caching empty volumes
    if (refreshNecessary_ && useCaching_.get() && VoreenApplication::app()->useCaching()) {
        cache_.store();
    }
    refreshNecessary_ = false;
}

std::vector<Volume*> VolumeChannelSeparator::separateChannels(const VolumeBase* volume) {

    std::string format = volume->getBaseType();

    //volume representations to be created
    std::vector<VolumeRAM*> separatedData;

    //create volume fusion
    if (format == "uint8")
        separatedData = separateVolumeChannels<uint8_t>(volume);
    else if (format == "int8")
        separatedData = separateVolumeChannels<int8_t>(volume);
    else if (format == "uint16")
        separatedData = separateVolumeChannels<uint16_t>(volume);
    else if (format == "int16")
        separatedData = separateVolumeChannels<int16_t>(volume);
    else if (format == "uint32")
        separatedData = separateVolumeChannels<uint32_t>(volume);
    else if (format == "int32")
        separatedData = separateVolumeChannels<int32_t>(volume);
    else if (format == "float")
        separatedData = separateVolumeChannels<float>(volume);
    else if (format == "double")
        separatedData = separateVolumeChannels<double>(volume);
    else {
        LERROR("Format currently not supported: " << format);
    }

    // create volumes and copy meta data
    std::vector<Volume*> result;
    for (size_t i = 0; i < separatedData.size(); ++i) {
        result.push_back(new Volume(separatedData.at(i), volume));
    }

    return result;
}

void VolumeChannelSeparator::refresh() {
    // only refresh when ready (do not buffer the refreshNecessary_ flag)
    if (!isInitialized() || !isReady())
        return;

    volumeOutport_.clear();
    volumeOutport2_.clear();
    volumeOutport3_.clear();
    volumeOutport4_.clear();
    refreshNecessary_ = true;
    invalidate();
}

}   // namespace
