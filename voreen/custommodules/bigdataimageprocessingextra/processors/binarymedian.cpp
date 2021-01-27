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

#include "binarymedian.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5filevolume.h"

#include "tgt/bounds.h"
#include "tgt/filesystem.h"

#include "modules/bigdataimageprocessing/volumefiltering/binarymedianfilter.h"

#include <vector>
#include <memory>

namespace voreen {

const std::string BinaryMedian::loggerCat_("voreen.bigdataimageprocessing.binarymedian");

BinaryMedian::BinaryMedian()
    : AsyncComputeProcessor<BinaryMedianInput, BinaryMedianOutput>()
    , inport_(Port::INPORT, "connectedcomponentanalysis.inport", "Binary Volume Input")
    , outport_(Port::OUTPORT, "connectedcomponentanalysis.outport", "Label Volume Output", false, Processor::VALID)
    , enabled_("enabled", "Enabled", true)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , outputVolumeDeflateLevel_("outputVolumeDeflateLevel", "Deflate Level", 1, 0, 9, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEFAULT)
    , binarizationThreshold_("binarizationThreshold", "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , samplingStrategyType_("samplingStrategyType", "Sampling Strategy", SamplingStrategyType::CLAMP_T)
    , outsideVolumeValue_("outsideVolumeValue", "Outside Volume Value", 0, 0, 1)
    , useUniformExtent_("useUniformExtent", "Uniform Extent", true)
    , extentX_("extentx", "Extent X", 1, 1, 10)
    , extentY_("extenty", "Extent Y", 1, 1, 10)
    , extentZ_("extentz", "Extent Z", 1, 1, 10)
    , kernelSizeDisplay_("kernelSizeDisplay", "Kernel Size", "")
    , forceMedian_("forceMedian", "Force Median", true)
    , objectVoxelThreshold_("objectVoxelThreshold", "Object Voxel Threshold", 0, 0, std::numeric_limits<int>::max())
{
    addPort(inport_);
        inport_.onChange(LambdaFunctionCallback([this] () {
                    outport_.setData(nullptr);
                    }));
    addPort(outport_);

        addProperty(enabled_);
            enabled_.setGroupID("output");
        addProperty(outputVolumeFilePath_);
            outputVolumeFilePath_.setGroupID("output");
        addProperty(outputVolumeDeflateLevel_);
            outputVolumeDeflateLevel_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

        addProperty(binarizationThreshold_);
            binarizationThreshold_.setGroupID("sampling");
        addProperty(samplingStrategyType_);
            samplingStrategyType_.setGroupID("sampling");
            samplingStrategyType_.addOption("clamp", "Clamp", SamplingStrategyType::CLAMP_T);
            samplingStrategyType_.addOption("mirror", "Mirror", SamplingStrategyType::MIRROR_T);
            samplingStrategyType_.addOption("set", "Set", SamplingStrategyType::SET_T);
            ON_CHANGE_LAMBDA(samplingStrategyType_, [this] () {
                    outsideVolumeValue_.setVisibleFlag(samplingStrategyType_.getValue() == SamplingStrategyType::SET_T);
                    });
        addProperty(outsideVolumeValue_);
            outsideVolumeValue_.setGroupID("sampling");
    setPropertyGroupGuiName("sampling", "Sampling");

        addProperty(useUniformExtent_);
            useUniformExtent_.setGroupID("filter");
            ON_CHANGE_LAMBDA(useUniformExtent_, [this] () {
                    if(useUniformExtent_.get()) {
                        extentY_.set(extentX_.get());
                        extentZ_.set(extentX_.get());
                    }
                    });
        addProperty(extentX_);
            extentX_.setGroupID("filter");
            ON_CHANGE_LAMBDA(extentX_, [this] () {
                    if(useUniformExtent_.get()) {
                        extentY_.set(extentX_.get());
                        extentZ_.set(extentX_.get());
                    }
                    updateObjectVoxelThreshold();
                    updateKernelSizeDisplay();
                    });
        addProperty(extentY_);
            extentY_.setGroupID("filter");
            ON_CHANGE_LAMBDA(extentY_, [this] () {
                    if(useUniformExtent_.get()) {
                        extentX_.set(extentY_.get());
                        extentZ_.set(extentY_.get());
                    }
                    updateObjectVoxelThreshold();
                    updateKernelSizeDisplay();
                    });
        addProperty(extentZ_);
            extentZ_.setGroupID("filter");
            ON_CHANGE_LAMBDA(extentZ_, [this] () {
                    if(useUniformExtent_.get()) {
                        extentX_.set(extentZ_.get());
                        extentY_.set(extentZ_.get());
                    }
                    updateObjectVoxelThreshold();
                    updateKernelSizeDisplay();
                    });
        addProperty(kernelSizeDisplay_);
            kernelSizeDisplay_.setGroupID("filter");
            kernelSizeDisplay_.setReadOnlyFlag(true);
        addProperty(forceMedian_);
            forceMedian_.setGroupID("filter");
            ON_CHANGE(forceMedian_, BinaryMedian, updateObjectVoxelThreshold);
        addProperty(objectVoxelThreshold_);
            objectVoxelThreshold_.setGroupID("filter");
    setPropertyGroupGuiName("filter", "Filter");

    // Update Outside Volume Value
    samplingStrategyType_.invalidate();
    // Kernel Size Display
    extentX_.invalidate();
}

void BinaryMedian::updateObjectVoxelThreshold() {
    bool medianForced = forceMedian_.get();
    objectVoxelThreshold_.setReadOnlyFlag(medianForced);
    objectVoxelThreshold_.setMaxValue((2*extentX_.get()+1)*(2*extentY_.get()+1)*(2*extentZ_.get()+1));
    if(medianForced) {
        objectVoxelThreshold_.set(objectVoxelThreshold_.getMaxValue()/2);
    }
}

void BinaryMedian::updateKernelSizeDisplay() {
    kernelSizeDisplay_.set(std::to_string(2*extentX_.get()+1)
            + "x" + std::to_string(2*extentY_.get()+1)
            + "x" + std::to_string(2*extentZ_.get()+1));
}

BinaryMedian::~BinaryMedian() {
}
VoreenSerializableObject* BinaryMedian::create() const {
    return new BinaryMedian();
}

void BinaryMedian::initialize() {
    AsyncComputeProcessor<BinaryMedianInput,BinaryMedianOutput>::initialize();
    updateObjectVoxelThreshold();
}
bool BinaryMedian::isReady() const {
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

BinaryMedianInput BinaryMedian::prepareComputeInput() {
    if(!enabled_.get()) {
        return BinaryMedianInput(
                nullptr,
                nullptr
                );
    }
    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto inputVolPtr = inport_.getThreadSafeData();
    const VolumeBase& inputVolume = *inputVolPtr;

    if(inputVolume.getNumChannels() != 1) {
        throw InvalidInputException("Input volume has multiple channels, but a single channel volume is expected!", InvalidInputException::S_ERROR);
    }


    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outport_.setData(nullptr);

    const std::string volumeFilePath = outputVolumeFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const std::string baseType = "uint8";
    const tgt::svec3 dim = inputVolume.getDimensions();

    RealWorldMapping rwm;
    if(inputVolume.hasMetaData("RealWorldMapping")) {
        rwm = inputVolume.getRealWorldMapping();
    }
    float normalizedBinarizationThreshold = rwm.realWorldToNormalized(binarizationThreshold_.get());

    if(volumeFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, dim, 1, true, outputVolumeDeflateLevel_.get(), tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    outputVolume->writeSpacing(inputVolume.getSpacing());
    outputVolume->writeOffset(inputVolume.getOffset());
    outputVolume->writeRealWorldMapping(RealWorldMapping(1,0,""));
    // For all zero or all one volumes the following is not correct,
    // and we cannot easily get the real min/max values without iterating
    // through the whole resulting volume.
    //const VolumeMinMax vmm(0, 1, 0, 1);
    //outputVolume->writeVolumeMinMax(&vmm);

    std::unique_ptr<SliceReader> sliceReader = VolumeFilterStackBuilder(inputVolume).addLayer(
            std::unique_ptr<VolumeFilter>(
                new BinaryMedianFilter(
                    tgt::ivec3(extentX_.get(), extentY_.get(), extentZ_.get()),
                    normalizedBinarizationThreshold,
                    objectVoxelThreshold_.get(),
                    SamplingStrategy<float>(samplingStrategyType_.getValue(), static_cast<float>(outsideVolumeValue_.get()))
                    )
                )
            ).build(0);

    return BinaryMedianInput(
        std::move(sliceReader),
        std::move(outputVolume)
    );
}
BinaryMedianOutput BinaryMedian::compute(BinaryMedianInput input, ProgressReporter& progressReporter) const {
    if(!enabled_.get()) {
        return { "" };
    }
    tgtAssert(input.sliceReader, "No sliceReader");
    tgtAssert(input.outputVolume, "No outputVolume");

    writeSlicesToHDF5File(*input.sliceReader, *input.outputVolume, &progressReporter);

    return {
        input.outputVolume->getFileName()
    };
    //outputVolume will be destroyed and thus closed now.
}
void BinaryMedian::processComputeOutput(BinaryMedianOutput output) {
    if(!enabled_.get()) {
        outport_.setData(inport_.getData(), false);
    } else {
        // outputVolume has been destroyed and thus closed by now.
        // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
        const VolumeBase* vol = HDF5VolumeReader().read(output.outputVolumeFilePath)->at(0);
        outport_.setData(vol);
    }
}

void BinaryMedian::adjustPropertiesToInput() {
    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }
    if(!input->hasDerivedData<VolumeMinMax>()) {
        LINFO("Calculating VolumeMinMax. This may take a while...");
    }
    const VolumeMinMax* mm = input->getDerivedData<VolumeMinMax>();

    binarizationThreshold_.setMinValue(mm->getMin());
    binarizationThreshold_.setMaxValue(mm->getMax());
    binarizationThreshold_.adaptDecimalsToRange(2);
}

} // namespace voreen
