/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "volumechannelmerger.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

const std::string VolumeChannelMerger::loggerCat_("voreen.VolumeChannelMerger");

VolumeChannelMerger::VolumeChannelMerger()
    : CachingVolumeProcessor()
    , volumeInport_(Port::INPORT, "volumeInput", "Volume Input")
    , volumeInport2_(Port::INPORT, "volumeInput2", "Volume Input 2")
    , volumeInport3_(Port::INPORT, "volumeInput3", "Volume Input 3")
    , volumeInport4_(Port::INPORT, "volumeInput4", "Volume Input 4")
    , volumeOutport_(Port::OUTPORT, "volumeOutport", "Volume Output")
    , autoRefresh_("autorefresh", "Auto Refresh", false)
    , refresh_("refreshbutton", "Refresh\n(clear and re-compute)")
    , refreshNecessary_(false)
{
    addPort(volumeInport_);
    addPort(volumeInport2_);
    addPort(volumeInport3_);
    addPort(volumeInport4_);
    addPort(volumeOutport_);

    addProperty(autoRefresh_);
    addProperty(refresh_);
    refresh_.onChange(MemberFunctionCallback<VolumeChannelMerger>(this, &VolumeChannelMerger::refresh));
}

Processor* VolumeChannelMerger::create() const {
    return new VolumeChannelMerger();
}

void VolumeChannelMerger::initialize() {
    CachingVolumeProcessor::initialize();
}

void VolumeChannelMerger::deinitialize() {
    volumeOutport_.clear();
    CachingVolumeProcessor::deinitialize();
}

bool VolumeChannelMerger::isReady() const {
    return ((volumeInport_.isReady() || volumeInport2_.isReady() || volumeInport3_.isReady() || volumeInport4_.isReady()) && volumeOutport_.isReady());
}

void VolumeChannelMerger::beforeProcess() {
    if (autoRefresh_.get())
        refreshNecessary_ = true;

    // if refreshing: look for results in cache
    if (refreshNecessary_)
        CachingVolumeProcessor::beforeProcess();
}

void VolumeChannelMerger::process() {

    setProgress(0.f);

    volumeOutport_.clear();

    //get the available volumes in input order and push into vector
    std::vector<const VolumeBase*> volumes;

    if (volumeInport_.getData())
        volumes.push_back(volumeInport_.getData());

    if (volumeInport2_.getData())
        volumes.push_back(volumeInport2_.getData());

    if (volumeInport3_.getData())
        volumes.push_back(volumeInport3_.getData());

    if (volumeInport4_.getData())
        volumes.push_back(volumeInport4_.getData());

    if(volumes.empty()) {
        return;
    }

    if (volumes.size() == 1) {
        //only one volume: just pass it through without further tests
        volumeOutport_.setData(volumes.at(0), false);
        LINFO("Only one input volume... passing through.");
        setProgress(1.f);
    }
    else {

        //check if all volumes are single-channel
        for (size_t i = 0; i < volumes.size(); ++i) {
            if (volumes.at(i)->getNumChannels() != 1) {
                LERROR("Volume at inport " << (i+1) << " is a multi-channel volume. Only single-channel volumes allowed.");
                return;
            }
        }

        //check data type dimensions, spacing, offset, transformation
        std::string type = volumes.at(0)->getFormat();
        tgt::svec3 dimensions = volumes.at(0)->getDimensions();
        tgt::vec3 spacing = volumes.at(0)->getSpacing();
        tgt::vec3 offset = volumes.at(0)->getOffset();
        tgt::mat4 transformation = volumes.at(0)->getPhysicalToWorldMatrix();

        for (size_t i = 1; i < volumes.size(); ++i) {

            if (type != volumes.at(i)->getFormat()) {
                LERROR("Input volumes do not have the same format.");
                return;
            }

            if (dimensions != volumes.at(i)->getDimensions()) {
                LERROR("Dimensions of input volumes differ.");
                return;
            }

            if (spacing != volumes.at(i)->getSpacing()) {
                LERROR("Spacing of input volumes differs.");
                return;
            }

            if (offset != volumes.at(i)->getOffset()) {
                LERROR("Offset of input volumes differs.");
                return;
            }

            if (transformation != volumes.at(i)->getPhysicalToWorldMatrix()) {
                LERROR("Transformations of input volumes differ.");
                return;
            }
        }

        //volumes are single-channel and compatible -> merge
        Volume* result = mergeVolumes(volumes);

        if (result)
            LINFO("Created multi-channel volume.");

        //set volume with ouport as owner so that it is correctly destroyed
        volumeOutport_.setData(result);

    }
}

void VolumeChannelMerger::afterProcess() {
    VolumeProcessor::afterProcess();

    // do not store the configuration if nothing happened to avoid problems with caching empty volumes
    if (refreshNecessary_ && useCaching_.get() && VoreenApplication::app()->useCaching()) {
        cache_.store();
    }
    refreshNecessary_ = false;
}

Volume* VolumeChannelMerger::mergeVolumes(std::vector<const VolumeBase*> volumes) {

    std::string format = volumes.at(0)->getFormat();

    //volume representation to be created
    VolumeRAM* mergedData = nullptr;

    //create volume fusion
    if (format == "uint8")
        mergedData = mergeVolumeRepresentations<uint8_t>(volumes);
    else if (format == "int8")
        mergedData = mergeVolumeRepresentations<int8_t>(volumes);
    else if (format == "uint16")
        mergedData = mergeVolumeRepresentations<uint16_t>(volumes);
    else if (format == "int16")
        mergedData = mergeVolumeRepresentations<int16_t>(volumes);
    else if (format == "uint32")
        mergedData = mergeVolumeRepresentations<uint32_t>(volumes);
    else if (format == "int32")
        mergedData = mergeVolumeRepresentations<int32_t>(volumes);
    else if (format == "float")
        mergedData = mergeVolumeRepresentations<float>(volumes);
    else if (format == "double")
        mergedData = mergeVolumeRepresentations<double>(volumes);
    else {
        LERROR("Format currently not supported: " << format);
        return nullptr;
    }

    if(!mergedData) {
        return nullptr;
    }

    // create volume and copy meta data
    Volume* result = new Volume(mergedData, volumes.at(0));
    return result;
}

void VolumeChannelMerger::refresh() {
    // only refresh when ready (do not buffer the refreshNecessary_ flag)
    if (!isInitialized() || !isReady())
        return;

    //volumeOutport_.clear();
    refreshNecessary_ = true;
    invalidate();
}

}   // namespace
