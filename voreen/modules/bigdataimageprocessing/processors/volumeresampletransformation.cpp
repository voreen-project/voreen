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

#include "volumeresampletransformation.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresample.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"

#include <climits>

namespace voreen {

const std::string VolumeResampleTransformation::loggerCat_("voreen.base.VolumeResampleTransformation");

VolumeResampleTransformation::VolumeResampleTransformation()
    : AsyncComputeProcessor<VolumeResampleTransformationInput, VolumeResampleTransformationOutput>()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , spacingHandling_("aspacingHandling", "Spacing Handling") //"a"-prefix is a hack so that spacingHandling_ is serialized before outputSpacing_ in XML-Serializer (which sorts the entriey alphabetically)
    , outputSpacing_("outputSpacing", "Spacing", tgt::vec3(1.0f), tgt::vec3(std::numeric_limits<float>::epsilon()), tgt::vec3(10.0f))
    , transformMatrix_("transformMatrix", "Transformation Matrix", tgt::mat4::identity, tgt::mat4(-1e10), tgt::mat4(1e10))
    , filteringMode_("filteringMode", "Filtering")
    , outsideVolumeHandling_("outsideVolumeHandling", "Outside Volume Handling")
    , outsideVolumeValue_("outsideVolumeValue", "Outside Volume Value", 0, std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max())
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , outputDimensions_("outputDimensions", "Resulting Dimensions", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(std::numeric_limits<int>::max()))
    , outputSizeMB_("outputSizeMB", "Output size (MB)")
{
    addPort(inport_);
    addPort(outport_);

        addProperty(spacingHandling_);
            spacingHandling_.setGroupID("config");
            spacingHandling_.addOption("keep", "Keep", KEEP);
            spacingHandling_.addOption("set", "Set", SET);
            spacingHandling_.addOption("multiply", "Multiply", MULTIPLY);
            spacingHandling_.selectByValue(KEEP);
            ON_CHANGE_LAMBDA(spacingHandling_, [this] () {
                        updateSpacingProperties();
                        adjustPropertiesToInput();
                    });
            updateSpacingProperties();

        addProperty(outputSpacing_);
            outputSpacing_.setTracking(false);
            outputSpacing_.setGroupID("config");
            ON_CHANGE(outputSpacing_, VolumeResampleTransformation, recalculateOutputProperties);
        addProperty(transformMatrix_);
            transformMatrix_.setGroupID("config");
            ON_CHANGE(transformMatrix_, VolumeResampleTransformation, recalculateOutputProperties);
        addProperty(filteringMode_);
            filteringMode_.setGroupID("config");
            filteringMode_.addOption("nearest", "Nearest", NEAREST);
            filteringMode_.addOption("linear",  "Linear", LINEAR);
            //filteringMode_.addOption("cubic",   "Cubic");
            filteringMode_.selectByValue(LINEAR);
        addProperty(outsideVolumeHandling_);
            outsideVolumeHandling_.setGroupID("config");
            outsideVolumeHandling_.addOption("clamp", "Clamp", CLAMP);
            outsideVolumeHandling_.addOption("constant_value",  "Constant Value", CONSTANT_VALUE);
            ON_CHANGE_LAMBDA(outsideVolumeHandling_, [this] () {
                    outsideVolumeValue_.setVisibleFlag(outsideVolumeHandling_.getValue() == CONSTANT_VALUE);
                    });
            outsideVolumeHandling_.selectByValue(CLAMP);
        addProperty(outsideVolumeValue_);
            outsideVolumeValue_.setGroupID("config");
        addProperty(outputVolumeFilePath_);
            outputVolumeFilePath_.setGroupID("config");
    setPropertyGroupGuiName("config", "Configuration");

        addProperty(outputDimensions_);
            outputDimensions_.setReadOnlyFlag(true);
            outputDimensions_.setGroupID("output");
        addProperty(outputSizeMB_);
            outputSizeMB_.setReadOnlyFlag(true);
            outputSizeMB_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output Information");
}

VolumeResampleTransformationInput VolumeResampleTransformation::prepareComputeInput() {
    const VolumeBase* inputPtr = inport_.getData();
    if(!inputPtr) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }
    const VolumeBase& input = *inputPtr;

    const std::string volumeFilePath = outputVolumeFilePath_.get();
    if(volumeFilePath.empty()) {
        throw InvalidInputException("No output volume file path specified!", InvalidInputException::S_ERROR);
    }

    tgt::Bounds inputBoundsPhysical(input.getLLF(), input.getURB());


    // Delete volume so that we can overwrite the .hdf5 file
    outport_.setData(nullptr);

    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const std::string baseType = input.getBaseType();
    const tgt::svec3 outputDim(getOutputDimensions(input));
    const int deflateLevel = 1;

    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, outputDim, 1, true, deflateLevel, tgt::svec3(outputDim.xy(), 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    tgt::vec3 spacing = getCurrentSpacing();
    tgt::vec3 offset = getOutputOffset(input);

    outputVolume->writeSpacing(spacing);
    outputVolume->writeOffset(offset);
    outputVolume->writePhysicalToWorldTransformation(input.getPhysicalToWorldMatrix());
    outputVolume->writeVolumeMinMax(input.getDerivedData<VolumeMinMax>());

    RealWorldMapping rwm = input.getRealWorldMapping();
    outputVolume->writeRealWorldMapping(rwm);

    tgt::mat4 invertedTransform;
    if(!transformMatrix_.get().invert(invertedTransform)) {
        throw InvalidInputException("Singular Transformation matrix!", InvalidInputException::S_ERROR);
    }

    tgt::mat4 outputVoxelToWorld = input.getPhysicalToWorldMatrix() * tgt::mat4::createTranslation(offset) * tgt::mat4::createScale(spacing);

    const tgt::mat4 voxelOutToVoxelIn = input.getWorldToVoxelMatrix() * invertedTransform * outputVoxelToWorld;

    tgt::ivec3 inputMaxVoxelPos(tgt::ivec3(input.getDimensions()) - tgt::ivec3(1));
    float outsideVolumeValueNormalized = rwm.realWorldToNormalized(outsideVolumeValue_.get());

    return VolumeResampleTransformationInput(
            input,
            std::move(outputVolume),
            voxelOutToVoxelIn,
            filteringMode_.getValue(),
            outsideVolumeHandling_.getValue(),
            outsideVolumeValueNormalized
            );
}
VolumeResampleTransformationOutput VolumeResampleTransformation::compute(VolumeResampleTransformationInput input, ProgressReporter& progress) const {
    tgt::ivec3 inputMaxVoxelPos(tgt::ivec3(input.input.getDimensions()) - tgt::ivec3(1));
    std::string outputFilename = input.outputVolume->getFileName();
    switch(input.filteringMode) {
        case NEAREST:
            switch(input.outsideVolumeHandling) {
                case CLAMP:
                    resample(input.input, std::move(input.outputVolume), input.voxelOutToVoxelIn, NearestFiltering<Clamping>(Clamping(inputMaxVoxelPos)), progress);
                    break;
                case CONSTANT_VALUE:
                    resample(input.input, std::move(input.outputVolume), input.voxelOutToVoxelIn, NearestFiltering<ConstantValue>(ConstantValue(inputMaxVoxelPos, input.outsideVolumeValue)), progress);
                    break;
            }
            break;
        case LINEAR:
            switch(outsideVolumeHandling_.getValue()) {
                case CLAMP:
                    resample(input.input, std::move(input.outputVolume), input.voxelOutToVoxelIn, LinearFiltering<Clamping>(Clamping(inputMaxVoxelPos)), progress);
                    break;
                case CONSTANT_VALUE:
                    resample(input.input, std::move(input.outputVolume), input.voxelOutToVoxelIn, LinearFiltering<ConstantValue>(ConstantValue(inputMaxVoxelPos, input.outsideVolumeValue)), progress);
                    break;
            }
            break;
    }

    return VolumeResampleTransformationOutput {
        outputFilename
    };
    //outputVolume will be destroyed and thus closed now.
}
void VolumeResampleTransformation::processComputeOutput(VolumeResampleTransformationOutput output) {
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    const VolumeBase* vol = HDF5VolumeReader().read(output.outputVolumeFilePath)->at(0);
    outport_.setData(vol);
}

bool VolumeResampleTransformation::isReady() const {
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

Processor* VolumeResampleTransformation::create() const {
    return new VolumeResampleTransformation();
}
tgt::vec3 VolumeResampleTransformation::getCurrentSpacing() const {
    switch(spacingHandling_.getValue()) {
        case KEEP:
            return outputSpacing_.get();
        case SET:
            return outputSpacing_.get();
        case MULTIPLY:
            if(inport_.hasData()) {
                return inport_.getData()->getSpacing() * outputSpacing_.get();
            } else {
                return outputSpacing_.get();
            }
    }
    return tgt::ivec3(0); //Should not happen
}

void VolumeResampleTransformation::adjustPropertiesToInput() {
    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }
    if(tgt::hmul(outputSpacing_.get()) == 0) {
        LERROR("Output spacing cannot be zero");
    }
    const VolumeMinMax* vmm = input->getDerivedData<VolumeMinMax>();
    tgtAssert(vmm, "No VolumeMinMax");
    outsideVolumeValue_.setMinValue(vmm->getMin());
    outsideVolumeValue_.setMaxValue(vmm->getMax());

    if(spacingHandling_.getValue() == KEEP) {
        outputSpacing_.set(input->getSpacing());
    }

    recalculateOutputProperties();
}

// private methods
//

void VolumeResampleTransformation::recalculateOutputProperties() {
    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }

    tgt::svec3 outputDim(getOutputDimensions(*input));

    outputDimensions_.set(tgt::ivec3(outputDim));
    size_t numVoxels = tgt::hmul(outputDim);
    outputSizeMB_.set(formatMemorySize(numVoxels*input->getBytesPerVoxel()));
}

void VolumeResampleTransformation::updateSpacingProperties() {
    switch(spacingHandling_.getValue()) {
        case KEEP:
            outputSpacing_.setReadOnlyFlag(true);
            if(inport_.hasData()) {
                outputSpacing_.set(inport_.getData()->getSpacing());
            }
            break;
        case SET:
            outputSpacing_.setReadOnlyFlag(false);
            if(inport_.hasData()) {
                outputSpacing_.set(inport_.getData()->getSpacing());
            }
            break;
        case MULTIPLY:
            outputSpacing_.setReadOnlyFlag(false);
            outputSpacing_.set(tgt::vec3(1,1,1));
            break;
    }
}

tgt::Bounds VolumeResampleTransformation::getPhysicalOutputBounds(const VolumeBase& input) const {
    // Bounding box of the input volume in voxel coordinates:
    // Notice the added 0.5 voxel wide border to account for the fact
    // that voxels are often considered cubes around their center
    // (i.e. for rendering).
    tgt::Bounds inputVoxelBB(tgt::vec3(-0.5), tgt::vec3(input.getDimensions()) - tgt::vec3(0.5));

    tgt::Bounds outputVoxelBB = inputVoxelBB.transform(
              input.getWorldToVoxelMatrix()
            * transformMatrix_.get()
            * input.getVoxelToWorldMatrix()
            );
    // Ensure that the (at least) first voxel of the volume is aligned with the original grid
    tgt::Bounds adjustetOutputVoxelBB(tgt::ifloor(outputVoxelBB.getLLF()), tgt::iceil(outputVoxelBB.getURB()));
    return adjustetOutputVoxelBB.transform(
              input.getVoxelToPhysicalMatrix()
            );
}
tgt::ivec3 VolumeResampleTransformation::getOutputDimensions(const VolumeBase& input) const {
    tgt::Bounds outputPhysicalBB = getPhysicalOutputBounds(input);
    return tgt::iceil(outputPhysicalBB.diagonal()/getCurrentSpacing()) + tgt::ivec3::one;
}
tgt::vec3 VolumeResampleTransformation::getOutputOffset(const VolumeBase& input) const {
    tgt::Bounds outputPhysicalBB = getPhysicalOutputBounds(input);
    return outputPhysicalBB.getLLF();
}

}   // namespace
