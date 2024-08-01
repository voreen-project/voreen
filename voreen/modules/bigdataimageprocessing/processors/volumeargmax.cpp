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

#include "volumeargmax.h"
#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"

namespace voreen {


const std::string VolumeArgMax::loggerCat_("voreen.RandomWalker.VolumeArgMax");

VolumeArgMax::VolumeArgMax()
    : AsyncComputeProcessor<VolumeArgMaxInput, VolumeArgMaxOutput>()
    , inportVolume0_(Port::INPORT, "volume.vol0")
    , inportVolume1_(Port::INPORT, "volume.vol1")
    , inportVolume2_(Port::INPORT, "volume.vol2")
    , inportVolume3_(Port::INPORT, "volume.vol3")
    , outportIds_(Port::OUTPORT, "volume.ids", "volume.ids", false)
    , outputVolumePath_("outputVolumePath", "ID Volume Output", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
{
    // ports
    addPort(inportVolume0_);
    addPort(inportVolume1_);
    addPort(inportVolume2_);
    addPort(inportVolume3_);
    addPort(outportIds_);

    // random walker properties
    addProperty(outputVolumePath_);
}

VolumeArgMax::~VolumeArgMax() {
}

Processor* VolumeArgMax::create() const {
    return new VolumeArgMax();
}

bool VolumeArgMax::isReady() const {
    return inportVolume0_.isReady();
}

VolumeArgMax::ComputeInput VolumeArgMax::prepareComputeInput() {
    tgtAssert(inportVolume0_.hasData(), "no input volume 0");

    const VolumeBase& vol0 = *inportVolume0_.getData();
    const tgt::svec3 dim = vol0.getDimensions();

    const VolumeBase* vol1 = inportVolume1_.getData();
    const VolumeBase* vol2 = inportVolume2_.getData();
    const VolumeBase* vol3 = inportVolume3_.getData();

    if(vol1 && vol1->getDimensions() != dim) {
        throw InvalidInputException("Dimension mismatch between volume 0 and volume 1", InvalidInputException::S_ERROR);
    }
    if(vol2 && vol2->getDimensions() != dim) {
        throw InvalidInputException("Dimension mismatch between volume 0 and volume 2", InvalidInputException::S_ERROR);
    }
    if(vol3 && vol3->getDimensions() != dim) {
        throw InvalidInputException("Dimension mismatch between volume 0 and volume 3", InvalidInputException::S_ERROR);
    }

    outportIds_.setData(nullptr);

    const std::string volumePath = outputVolumePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;

    const std::string baseType = "uint8";
    const tgt::vec3 spacing = vol0.getSpacing();
    const tgt::vec3 offset = vol0.getOffset();
    const RealWorldMapping rwm = RealWorldMapping::createDenormalizingMapping(baseType);
    const int deflateLevel = 1;

    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumePath, volumeLocation, baseType, dim, 1, true, deflateLevel, tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    outputVolume->writeSpacing(spacing);
    outputVolume->writeOffset(offset);
    outputVolume->writeRealWorldMapping(rwm);

    return VolumeArgMaxInput {
        vol0,
        vol1,
        vol2,
        vol3,
        std::move(outputVolume),
    };
}

static void sliceArgmax(const VolumeBase& vol, size_t z, uint8_t id, VolumeAtomic<uint8_t>& ids, VolumeAtomic<float>& maxVals) {
    const tgt::svec3 dim = vol.getDimensions();
    std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));
    RealWorldMapping rwm = vol.getRealWorldMapping();
    for(size_t y=0; y<dim.y; ++y) {
        for(size_t x=0; x<dim.x; ++x) {
            float curVal = rwm.normalizedToRealWorld(inputSlice->getVoxelNormalized(x,y,0));
            float maxVal = maxVals.voxel(x,y,0);
            if(curVal > maxVal) {
                maxVals.setVoxelNormalized(curVal, x,y,0);
                ids.voxel(x,y,0) = id;
            }
        }
    }
}

static VolumeAtomic<float> getSliceWithRWM(const VolumeBase& vol, size_t z) {
    RealWorldMapping rwm = vol.getRealWorldMapping();
    std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));
    VolumeAtomic<float> sliceOut(inputSlice->getDimensions());

    size_t numVoxels = sliceOut.getNumVoxels();
    for(int i=0; i<numVoxels; ++i) {
        float v = rwm.normalizedToRealWorld(inputSlice->getVoxelNormalized(i));
        sliceOut.voxel(i) = v;
    }
    return sliceOut;
}

VolumeArgMax::ComputeOutput VolumeArgMax::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    progressReporter.setProgress(0.0);

    const tgt::svec3 dim = input.vol0_.getDimensions();
    tgtAssert(!input.vol1_ || input.vol1_->getDimensions() == dim, "dim mismatch");
    tgtAssert(!input.vol2_ || input.vol2_->getDimensions() == dim, "dim mismatch");
    tgtAssert(!input.vol3_ || input.vol3_->getDimensions() == dim, "dim mismatch");

    for(size_t z=0; z<dim.z; ++z) {
        VolumeAtomic<uint8_t> sliceId(tgt::vec3(dim.x, dim.y, 1));
        sliceId.fill(0);
        VolumeAtomic<float> sliceMaxVals = getSliceWithRWM(input.vol0_, z);

        if(input.vol1_) {
            sliceArgmax(*input.vol1_, z, 1, sliceId, sliceMaxVals);
        }
        if(input.vol2_) {
            sliceArgmax(*input.vol2_, z, 2, sliceId, sliceMaxVals);
        }
        if(input.vol3_) {
            sliceArgmax(*input.vol3_, z, 3, sliceId, sliceMaxVals);
        }

        input.outputVolume_->writeSlices(&sliceId, z);
        progressReporter.setProgress(static_cast<float>(z+1)/dim.z);
    }

    return {
        input.outputVolume_->getFileName(),
    };
}

void VolumeArgMax::processComputeOutput(ComputeOutput output) {
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    std::unique_ptr<VolumeList> volumes(HDF5VolumeReaderOriginal().read(output.outputVolumePath_));
    const VolumeBase* vol = volumes->at(0);

    outportIds_.setData(vol);
}

}   // namespace
