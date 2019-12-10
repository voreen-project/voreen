/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "volumeformatconversion.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"

#include "tgt/vector.h"

namespace voreen {

const std::string VolumeFormatConversion::loggerCat_("voreen.base.VolumeFormatConversion");

VolumeFormatConversion::VolumeFormatConversion()
    : CachingVolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output",false),
    enableProcessing_("enabled", "Enable", false),
    targetFormat_("targetFormat", "Target Data Type"),
    numChannels_("numChannels", "Num Channels", 0, 0, 4)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);

    targetFormat_.addOption("uint8_t",  "8 Bit Unsigned Integer (uint8)");
    targetFormat_.addOption("int8_t",   "8 Bit Signed Integer (int8)");
    targetFormat_.addOption("uint16_t", "16 Bit Unsigned Integer (uint16)");
    targetFormat_.addOption("int16_t",  "16 Bit Signed Integer (int16)");
    targetFormat_.addOption("uint32_t", "32 Bit Unsigned Integer (uint32)");
    targetFormat_.addOption("int32_t",  "32 Bit Signed Integer (int32)");
    targetFormat_.addOption("uint64_t", "64 Bit Unsigned Integer (uint64)");
    targetFormat_.addOption("int64_t",  "64 Bit Signed Integer (int64)");
    targetFormat_.addOption("float",    "Float");
    targetFormat_.addOption("double",   "Double");
    addProperty(targetFormat_);

    numChannels_.setReadOnlyFlag(true);
    addProperty(numChannels_);
}

VolumeFormatConversion::~VolumeFormatConversion() {}

Processor* VolumeFormatConversion::create() const {
    return new VolumeFormatConversion();
}

void VolumeFormatConversion::process() {
    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    tgtAssert(inputVolume, "no volume");

    size_t numChannels = inputVolume->getNumChannels();
    numChannels_.set((int)numChannels);

    if (!enableProcessing_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
        return;
    }

    Volume* outputVolume = 0;

    VolumeOperatorConvert voConvert;
    voConvert.setProgressReporter(this);
    if (numChannels == 1) {
        if (targetFormat_.isSelected("uint8_t")) {
            outputVolume = voConvert.apply<uint8_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("int8_t")) {
            outputVolume = voConvert.apply<int8_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint16_t")) {
            outputVolume = voConvert.apply<uint16_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("int16_t")) {
            outputVolume = voConvert.apply<int16_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint32_t")) {
            outputVolume = voConvert.apply<uint32_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("int32_t")) {
            outputVolume = voConvert.apply<int32_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint64_t")) {
            outputVolume = voConvert.apply<uint64_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("int64_t")) {
            outputVolume = voConvert.apply<int64_t>(inport_.getData());
        }
        else if (targetFormat_.isSelected("float")) {
            outputVolume = voConvert.apply<float>(inport_.getData());
        }
        else if (targetFormat_.isSelected("double")) {
            outputVolume = voConvert.apply<double>(inport_.getData());
        }
        else {
            LERROR("Unknown target format: " << targetFormat_.get());
        }
    }
    else if (numChannels == 2) {
        if (targetFormat_.isSelected("uint8_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<uint8_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int8_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<int8_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint16_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<uint16_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int16_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<int16_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint32_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<uint32_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int32_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<int32_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint64_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<uint64_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int64_t")) {
            outputVolume = voConvert.apply<tgt::Vector2<int64_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("float")) {
            outputVolume = voConvert.apply<tgt::Vector2<float> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("double")) {
            outputVolume = voConvert.apply<tgt::Vector2<double> >(inport_.getData());
        }
        else {
            LERROR("Unknown target format: " << targetFormat_.get());
        }
    }
    else if (numChannels == 3) {
        if (targetFormat_.isSelected("uint8_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<uint8_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int8_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<int8_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint16_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<uint16_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int16_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<int16_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint32_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<uint32_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int32_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<int32_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint64_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<uint64_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int64_t")) {
            outputVolume = voConvert.apply<tgt::Vector3<int64_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("float")) {
            outputVolume = voConvert.apply<tgt::Vector3<float> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("double")) {
            outputVolume = voConvert.apply<tgt::Vector3<double> >(inport_.getData());
        }
        else {
            LERROR("Unknown target format: " << targetFormat_.get());
        }
    }
    else if (numChannels == 4) {
        if (targetFormat_.isSelected("uint8_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<uint8_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int8_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<int8_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint16_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<uint16_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int16_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<int16_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint32_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<uint32_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int32_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<int32_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("uint64_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<uint64_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("int64_t")) {
            outputVolume = voConvert.apply<tgt::Vector4<int64_t> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("float")) {
            outputVolume = voConvert.apply<tgt::Vector4<float> >(inport_.getData());
        }
        else if (targetFormat_.isSelected("double")) {
            outputVolume = voConvert.apply<tgt::Vector4<double> >(inport_.getData());
        }
        else {
            LERROR("Unknown target format: " << targetFormat_.get());
        }
    }
    else {
        LERROR("Unsupported channel count: " << numChannels);
    }

    // assign computed volume to outport
    outport_.setData(outputVolume);
}

}   // namespace
