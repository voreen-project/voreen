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

#include "volumemultiplier.h"

#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"

#include "tgt/bounds.h"
#include "tgt/filesystem.h"

#include <vector>
#include <memory>

namespace voreen {

const std::string VolumeMultiplier::loggerCat_("voreen.bigdataimageprocessing.binarymedian");

VolumeMultiplier::VolumeMultiplier()
    : AsyncComputeProcessor<VolumeMultiplierInput, VolumeMultiplierOutput>()
    , inport_(Port::INPORT, "connectedcomponentanalysis.inport", "Binary Volume Input")
    , outport_(Port::OUTPORT, "connectedcomponentanalysis.outport", "Label Volume Output", false, Processor::VALID)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , outputVolumeDeflateLevel_("outputVolumeDeflateLevel", "Deflate Level", 1, 0, 9, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEFAULT)
    , multiplicationFactor_("multiplicationFactor", "Multiplication Factor", tgt::ivec3(1), tgt::ivec3(1), tgt::ivec3(100))
    , outputSizeDisplay_("outputSize", "Output Size", "")
{
    addPort(inport_);
        inport_.onChange(LambdaFunctionCallback([this] () {
                    updateOutputSizeDisplay();
                    outport_.setData(nullptr);
                    }));
    addPort(outport_);

        addProperty(outputVolumeFilePath_);
            outputVolumeFilePath_.setGroupID("output");
        addProperty(outputVolumeDeflateLevel_);
            outputVolumeDeflateLevel_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

        addProperty(multiplicationFactor_);
            multiplicationFactor_.setGroupID("multiplication");
            ON_CHANGE(multiplicationFactor_, VolumeMultiplier, updateOutputSizeDisplay);
        addProperty(outputSizeDisplay_);
            outputSizeDisplay_.setGroupID("multiplication");
            outputSizeDisplay_.setReadOnlyFlag(true);
    setPropertyGroupGuiName("multiplication", "Multiplication");
}


void VolumeMultiplier::updateOutputSizeDisplay() {
    std::stringstream output;
    if(inport_.hasData()) {
        const VolumeBase& inputVolume = *inport_.getData();
        auto voxels = inputVolume.getDimensions() * tgt::svec3(multiplicationFactor_.get());
        output << voxels.x << "x" << voxels.y << "x" << voxels.z << "voxels, ";
        auto size = tgt::hmul(voxels) * inputVolume.getBytesPerVoxel();
        output << formatMemorySize(size) << " (in memory)";
    } else {
        output << "No input volume";
    }
    outputSizeDisplay_.set(output.str());
}

VolumeMultiplier::~VolumeMultiplier() {
}
VoreenSerializableObject* VolumeMultiplier::create() const {
    return new VolumeMultiplier();
}

void VolumeMultiplier::initialize() {
    AsyncComputeProcessor<VolumeMultiplierInput,VolumeMultiplierOutput>::initialize();
    updateOutputSizeDisplay();
}
bool VolumeMultiplier::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if(!inport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }
    return true;
}

VolumeMultiplierInput VolumeMultiplier::prepareComputeInput() {
    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto inputVolumePtr = inport_.getThreadSafeData();
    const VolumeBase& inputVolume = *inputVolumePtr;

    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outport_.setData(nullptr);

    const std::string volumeFilePath = outputVolumeFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const std::string baseType = inputVolume.getBaseType();
    const tgt::svec3 dim = inputVolume.getDimensions() * tgt::svec3(multiplicationFactor_.get());
    const size_t numChannels = inputVolume.getNumChannels();

    RealWorldMapping rwm;
    if(inputVolume.hasMetaData("RealWorldMapping")) {
        rwm = inputVolume.getRealWorldMapping();
    }

    if(volumeFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, dim, numChannels, true, outputVolumeDeflateLevel_.get(), tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    outputVolume->writeSpacing(inputVolume.getSpacing());
    outputVolume->writeOffset(inputVolume.getOffset());
    outputVolume->writeRealWorldMapping(inputVolume.getRealWorldMapping());

    return VolumeMultiplierInput(
        inputVolume,
        std::move(outputVolume)
    );
}

/*
int mirror(int i, int d) {
    int even_odd_i = (i < 0) ? std::abs(i+1) : i;
    if ((((even_odd_i/d)%2)==0) ^ (i < 0)) {
        return std::abs(i)%d;
    } else {
        return d-1 - (std::abs(i)%d);
    }
}
*/
size_t mirror(size_t i, size_t d) {
    if (((i/d)%2)==0) {
        return i%d;
    } else {
        return d-1 - (i%d);
    }
}

VolumeMultiplierOutput VolumeMultiplier::compute(VolumeMultiplierInput input, ProgressReporter& progressReporter) const {
    tgtAssert(input.outputVolume, "No outputVolume");

    const auto outputdim = input.outputVolume->getDimensions();
    const auto inputdim = input.inputVolume.getDimensions();
    const auto numChannels = input.inputVolume.getNumChannels();
    const auto format = input.inputVolume.getFormat();
    VolumeFactory factory;

    for(size_t z=0; z < outputdim.z; ++z) {
        progressReporter.setProgress(static_cast<float>(z)/outputdim.z);

        std::unique_ptr<VolumeRAM> inputslice(input.inputVolume.getSlice(mirror(z, inputdim.z)));
        std::unique_ptr<VolumeRAM> outputslice(factory.create(format, tgt::svec3(outputdim.x, outputdim.y, 1))); //This will probably not work, need a factory

        for(size_t y=0; y < outputdim.y; ++y) {
            for(size_t x=0; x < outputdim.x; ++x) {
                for(size_t c=0; c < numChannels; ++c) {
                    float val = inputslice->getVoxelNormalized(tgt::svec3(mirror(x, inputdim.x), mirror(y, inputdim.y), 0), c);
                    outputslice->setVoxelNormalized(val, tgt::svec3(x, y, 0), c);
                }
            }
        }
        input.outputVolume->writeSlices(outputslice.get(), z);
    }

    // As the size of the multiplied volume is larger than that of this one,
    // it's worth it to get the VMM of the input volume even if it has to be calculated.
    VolumeMinMax* vmm = input.inputVolume.getDerivedData<VolumeMinMax>(); //NOT owned, so we do not need to delete it
    input.outputVolume->writeVolumeMinMax(vmm);

    progressReporter.setProgress(1.0f);

    return {
        input.outputVolume->getFileName()
    };
    //outputVolume will be destroyed and thus closed now.
}
void VolumeMultiplier::processComputeOutput(VolumeMultiplierOutput output) {
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    std::unique_ptr<VolumeList> volumes(HDF5VolumeReader().read(output.outputVolumeFilePath));
    const VolumeBase* vol = volumes->at(0);
    outport_.setData(vol);
}

} // namespace voreen
