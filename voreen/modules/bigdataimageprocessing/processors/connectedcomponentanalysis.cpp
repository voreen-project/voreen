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

#include "connectedcomponentanalysis.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "modules/hdf5/io/hdf5volumewriter.h"

#include "tgt/bounds.h"
#include "tgt/filesystem.h"

#include <vector>
#include <memory>

namespace voreen {

CCANodeMetaData::CCANodeMetaData()
    : volume_(0)
    , bounds_()
{
}
CCANodeMetaData::CCANodeMetaData(tgt::svec2 yzPos, size_t lowerBound, size_t upperBound)
    : volume_(upperBound-lowerBound)
    , bounds_(tgt::svec3(lowerBound, yzPos.x, yzPos.y),
                tgt::svec3(upperBound-1, yzPos.x, yzPos.y))
{
}
CCANodeMetaData& CCANodeMetaData::operator+=(const CCANodeMetaData& rhs) {
    volume_ += rhs.volume_;
    bounds_.addVolume(rhs.bounds_);
    return *this;
}

const std::string ConnectedComponentAnalysis::loggerCat_("voreen.bigdataimageprocessing.connectedcomponentanalysis");

ConnectedComponentAnalysis::ConnectedComponentAnalysis()
    : AsyncComputeProcessor()
    , inport_(Port::INPORT, "connectedcomponentanalysis.inport", "Binary Volume Input")
    , outport_(Port::OUTPORT, "connectedcomponentanalysis.outport", "Label Volume Output", false, Processor::VALID)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , componentStatFilePath_("componentStatFilePath", "Component Stat File", "Path", "", "CSV (*.csv)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , writeComponentStatFile_("writeComponentStatFile", "Write Component Stat File", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , outputVolumeDeflateLevel_("outputVolumeDeflateLevel", "Deflate Level", 1, 0, 9, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEFAULT)
    , neighbourhoodMode_("neighbourhoodMode", "Neighborhood")
    , invertBinarization_("invertBinarization", "Invert", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , binarizationThreshold_("binarizationThreshold", "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , minBoundsDiagonal_("minBoundsDiagonal", "Min Bounds Diagonal (mm)", 0.0f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEFAULT)
    , minBoundsDiagonalRelative_("minBoundsDiagonalRelative", "Min Bounds Diagonal (relative)", 0.0f, 0.0f, 1.0f, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEFAULT)
    , minVoxelVolume_("minVoxelVolume", "Min Voxel Volume", 0, 0, std::numeric_limits<int>::max(), Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEFAULT)
    , applyLabeling_("applyLabeling", "Apply Labeling", true, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
{
    addPort(inport_);
        inport_.onChange(LambdaFunctionCallback([this] () {
                    outport_.setData(nullptr);
                    }));
    addPort(outport_);

        addProperty(outputVolumeFilePath_);
            outputVolumeFilePath_.setGroupID("output");
        addProperty(outputVolumeDeflateLevel_);
            outputVolumeDeflateLevel_.setGroupID("output");
        addProperty(componentStatFilePath_);
            componentStatFilePath_.setGroupID("output");
        addProperty(writeComponentStatFile_);
            writeComponentStatFile_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

        addProperty(invertBinarization_);
            invertBinarization_.setGroupID("binarization");
        addProperty(binarizationThreshold_);
            binarizationThreshold_.setGroupID("binarization");
    setPropertyGroupGuiName("binarization", "Binarization");

        addProperty(neighbourhoodMode_);
            neighbourhoodMode_.addOption("n6", "6", N_6);
            neighbourhoodMode_.addOption("n18", "18", N_18);
            neighbourhoodMode_.addOption("n26", "26", N_26);
            neighbourhoodMode_.setGroupID("componentfinding");
        addProperty(applyLabeling_);
            applyLabeling_.setGroupID("componentfinding");
        addProperty(minBoundsDiagonalRelative_);
            minBoundsDiagonalRelative_.setGroupID("componentfinding");
        addProperty(minBoundsDiagonal_);
            minBoundsDiagonal_.setGroupID("componentfinding");
        addProperty(minVoxelVolume_);
            minVoxelVolume_.setGroupID("componentfinding");
    setPropertyGroupGuiName("componentfinding", "Component Finding");
}

ConnectedComponentAnalysis::~ConnectedComponentAnalysis() {
}
VoreenSerializableObject* ConnectedComponentAnalysis::create() const {
    return new ConnectedComponentAnalysis();
}

bool ConnectedComponentAnalysis::isReady() const {
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


CCAComputeInput ConnectedComponentAnalysis::prepareComputeInput() {
    if (!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto inputVolPtr = inport_.getThreadSafeData();
    const VolumeBase& inputVolume = *inputVolPtr;

    if(inputVolume.getNumChannels() != 1) {
        throw InvalidInputException("Input volume has multiple channels, but a single channel volume is expected!", InvalidInputException::S_ERROR);
    }

    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outport_.setData(nullptr);

    const std::string baseType = applyLabeling_.get() ? "uint32" : "uint8";
    const std::string volumeFilePath = outputVolumeFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    if(volumeFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    // We need the statWriter to live as long as the lambda,
    // but we cannot move it into the lambda because we need c++14 to do that.
    // Workaround: Keep an (possible null) unique_ptr in the function scope
    // that lives longer than writeMetaData
    std::unique_ptr<CCAWriterType> statWriter = nullptr;
    std::function<void(int id, const CCANodeMetaData&)> writeMetaData;

    if(writeComponentStatFile_.get()) {
        if(!applyLabeling_.get()) {
            throw InvalidInputException("Cannot write stats if labeling is disabled!", InvalidInputException::S_ERROR);
        }
        if(componentStatFilePath_.get().empty()) {
            throw InvalidInputException("No component stat file path specified!", InvalidInputException::S_ERROR);
        }

        try {
            statWriter = std::unique_ptr<CCAWriterType>(new CCAWriterType(componentStatFilePath_.get(), ','));
        } catch (tgt::IOException& e) {
            std::string error = std::string("Could not create component stat file: ") + e.what();
            throw InvalidInputException(error, InvalidInputException::S_ERROR);
        }
        //statWriter->writeHeader("id", "volume", "llf.x", "llf.y", "llf.z", "urb.x", "urb.y", "urb.z");

        // Workaround (see above) pass the pointer to the statwriter into the lambda.
        CCAWriterType* statWriterPtr = statWriter.get();
        writeMetaData = [statWriterPtr] (uint32_t id, const CCANodeMetaData& m) {
            tgt::vec3 llf = m.bounds_.getLLF();
            tgt::vec3 urb = m.bounds_.getURB();
            statWriterPtr->write(id, m.volume_, llf.x, llf.y, llf.z, urb.x, urb.y, urb.z);
        };
    } else {
        writeMetaData = [] (uint32_t, const CCANodeMetaData&) {};
    }

    const tgt::svec3 dim = inputVolume.getDimensions();
    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, dim, 1, true, outputVolumeDeflateLevel_.get(), tgt::svec3(dim.x, dim.y, 1), false);
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    return CCAComputeInput(
        std::move(statWriter),
        writeMetaData,
        inputVolume,
        std::move(outputVolume),
        neighbourhoodMode_.getValue()
    );
}

CCAComputeOutput ConnectedComponentAnalysis::compute(CCAComputeInput input, ProgressReporter& progressReporter) const {
    tgtAssert(input.outputVolume, "No output volume");
    HDF5FileVolume& outputVolume = *input.outputVolume;

    CCAComputeOutput output;

    switch(input.neighbourhoodMode) {
        case N_26:
            output.stats = runCCA<0>(input.inputVolume, outputVolume, input.writeMetaData, progressReporter);
            break;
        case N_18:
            output.stats = runCCA<1>(input.inputVolume, outputVolume, input.writeMetaData, progressReporter);
            break;
        case N_6:
            output.stats = runCCA<2>(input.inputVolume, outputVolume, input.writeMetaData, progressReporter);
            break;
        default:
            tgtAssert(false, "Invalid Neighborhood Mode");
    }
    return output;
}

static float getDiagonalLength(const VolumeBase& vol) {
    return tgt::length(tgt::vec3(vol.getDimensions())*vol.getSpacing());
}

void ConnectedComponentAnalysis::processComputeOutput(CCAComputeOutput output) {
    if(applyLabeling_.get()) {
        LINFO("Used " << output.stats.numComponents << " ids for " << output.stats.numVoxels << " voxels.");
    } else {
        LINFO("Marked " << output.stats.numVoxels << " voxels.");
    }

    const std::string volumeFilePath = outputVolumeFilePath_.get();
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    const VolumeBase* vol = HDF5VolumeReader().read(volumeFilePath)->at(0);
    outport_.setData(vol);
}
std::function<bool(const CCANodeMetaData&)> ConnectedComponentAnalysis::generateComponentConstraintTest(const VolumeBase& originVolume) const {
    const tgt::vec3 spacing = originVolume.getSpacing();
    const float minBoundsDiagonal = std::max(
            minBoundsDiagonal_.get(),
            minBoundsDiagonalRelative_.get() * getDiagonalLength(originVolume)
            );
    const size_t minVoxelVolume = static_cast<size_t>(minVoxelVolume_.get());

    return [spacing, minBoundsDiagonal, minVoxelVolume] (const CCANodeMetaData& metaData) {
        tgt::vec3 voxelDiagonal(metaData.bounds_.diagonal());
        return tgt::lengthSq(voxelDiagonal*spacing) >= minBoundsDiagonal*minBoundsDiagonal
            && metaData.volume_ >= minVoxelVolume;
    };
}


void ConnectedComponentAnalysis::adjustPropertiesToInput() {
    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }
    if(!input->hasDerivedData<VolumeMinMax>()) {
        LINFO("Calculating VolumeMinMax. This may take a while...");
    }
    const VolumeMinMax* mm = input->getDerivedData<VolumeMinMax>();

    // We cannot display values larger than int::max in IntProperties.
    // In most cases it also should not be necessary to set such a large threshold.
    // If you do find yourself in need of it, consider changing minVoxelVolume_
    // to an UInt64Property (if it is implemented already).
    minVoxelVolume_.setMaxValue(static_cast<int>(std::min(static_cast<size_t>(std::numeric_limits<int>::max()), input->getNumVoxels())));

    minBoundsDiagonal_.setMaxValue(getDiagonalLength(*input));
    minBoundsDiagonal_.adaptDecimalsToRange(2);

    binarizationThreshold_.setMinValue(mm->getMin());
    binarizationThreshold_.setMaxValue(mm->getMax());
    binarizationThreshold_.adaptDecimalsToRange(4);

    // Volumes without an a RWM are expected to have normalized max and min values (and normalized values in general).
    tgtAssert(input->hasMetaData("RealWorldMapping") || mm->getMin() == mm->getMinNormalized() && mm->getMax() == mm->getMaxNormalized(), "Missing RealWorldMapping or invalid VolumeMinMax");
}

} // namespace voreen
