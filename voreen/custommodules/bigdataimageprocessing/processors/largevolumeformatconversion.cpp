/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "largevolumeformatconversion.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include "tgt/vector.h"

namespace voreen {

const std::string LargeVolumeFormatConversion::loggerCat_("voreen.bigdataimageprocessing.LargeVolumeFormatConversion");

LargeVolumeFormatConversion::LargeVolumeFormatConversion()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output",false)
    , enableProcessing_("enabled", "Enable", false)
    , targetBaseType_("targetFormat", "Target Data Type")
    , numChannels_("numChannels", "Num Channels", 0, 0, 4)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output volume file path", "Output volume file path", "", LZ4SliceVolumeBase::FILE_EXTENSION)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);

    targetBaseType_.addOption("uint8",  "8 Bit Unsigned Integer (uint8)");
    targetBaseType_.addOption("int8",   "8 Bit Signed Integer (int8)");
    targetBaseType_.addOption("uint16", "16 Bit Unsigned Integer (uint16)");
    targetBaseType_.addOption("int16",  "16 Bit Signed Integer (int16)");
    targetBaseType_.addOption("uint32", "32 Bit Unsigned Integer (uint32)");
    targetBaseType_.addOption("int32",  "32 Bit Signed Integer (int32)");
    targetBaseType_.addOption("float",    "Float");
    targetBaseType_.addOption("double",   "Double");
    addProperty(targetBaseType_);


    numChannels_.setReadOnlyFlag(true);
    addProperty(numChannels_);

    addProperty(outputVolumeFilePath_);
}

LargeVolumeFormatConversion::~LargeVolumeFormatConversion() {}

Processor* LargeVolumeFormatConversion::create() const {
    return new LargeVolumeFormatConversion();
}

template<typename OutputFormat>
void processDispatch(const VolumeBase& input, std::unique_ptr<Volume>& output, const std::string& outputPath, ProgressReporter& progressReporter) {
    float scale;
    float offset;

    {
        VolumeAtomic<OutputFormat> outputMetadata(tgt::svec3::one);

        // determine input mapping range
        auto vmm = input.getDerivedData<VolumeMinMax>();
        float min =  std::numeric_limits<float>::infinity();
        float max = -std::numeric_limits<float>::infinity();
        for(size_t i = 0; i < vmm->getNumChannels(); ++i) {
            min = std::min(min, vmm->getMinNormalized(i));
            max = std::max(max, vmm->getMaxNormalized(i));
        }

        // map [min:max] to output range
        if (outputMetadata.isInteger() && outputMetadata.isSigned()) {
            // map [min:max] to [-1.0:1.0]
            offset = -min - 1.f;
            scale = 1.f / (max - min + 1.f);
        } else {
            // map [min:max] to [0.0:1.0]
            offset = -min;
            scale = 1.f / (max - min);
        }
    }

    tgt::svec3 dim = input.getDimensions();

    // real world mapping has to revert the applied value transformation
    RealWorldMapping destMapping = RealWorldMapping::combine(RealWorldMapping(1.f/scale, -offset/scale, ""), input.getRealWorldMapping());

    LZ4SliceVolumeBuilder<OutputFormat> builder(outputPath,
            LZ4SliceVolumeMetadata(dim)
            .withOffset(input.getOffset())
            .withSpacing(input.getSpacing())
            .withPhysicalToWorldTransformation(input.getPhysicalToWorldMatrix())
            .withRealWorldMapping(destMapping));
    size_t numChannels = input.getNumChannels();

    for(size_t z = 0; z < dim.z; ++z) {
        progressReporter.setProgress(static_cast<float>(z)/dim.z);
        std::unique_ptr<VolumeRAM> inputSlice(input.getSlice(z));
        auto outputSlice = builder.getNextWriteableSlice();
        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::svec3 slicePos(x,y,0);
                for(size_t c = 0; c < numChannels; ++c) {
                    float val = inputSlice->getVoxelNormalized(slicePos, c)*scale + offset;
                    outputSlice->setVoxelNormalized(val, slicePos, c);
                }
            }
        }
    }
    progressReporter.setProgress(1.f);

    output = std::move(builder).finalize().toVolume();
}

LargeVolumeFormatConversion::ComputeInput LargeVolumeFormatConversion::prepareComputeInput() {
    if(!enableProcessing_.get()) {
        return LargeVolumeFormatConversion::ComputeInput { "", "", nullptr };
    }

    if(outputVolumeFilePath_.get().empty()) {
        throw InvalidInputException("No volume path specified", InvalidInputException::S_ERROR);
    }

    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    return LargeVolumeFormatConversion::ComputeInput {
        outputVolumeFilePath_.get(),
        targetBaseType_.get(),
        inport_.getData()
    };
}

LargeVolumeFormatConversion::ComputeOutput LargeVolumeFormatConversion::compute(LargeVolumeFormatConversion::ComputeInput input, ProgressReporter& progressReporter) const {
    if(!enableProcessing_.get()) {
        return { nullptr };
    }
    tgtAssert(input.inputVolume_, "No volume");

    std::string outputFormat = getFormatFromBaseTypeAndChannels(targetBaseType_.get(), input.inputVolume_->getNumChannels());

    std::unique_ptr<Volume> outputVolume(nullptr);

    DISPATCH_FOR_FORMAT(outputFormat, processDispatch, *input.inputVolume_, outputVolume, input.outputPath_, progressReporter);

    return {
        std::move(outputVolume)
    };
}
void LargeVolumeFormatConversion::processComputeOutput(LargeVolumeFormatConversion::ComputeOutput output) {
    if(!enableProcessing_.get()) {
        outport_.setData(inport_.getData(), false);
    } else {
        outport_.setData(output.outputVolume_.release());
    }
}

}   // namespace
