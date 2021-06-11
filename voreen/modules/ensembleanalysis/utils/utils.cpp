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

#include "utils.h"

#include "voreen/core/io/volumeserializer.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumeurl.h"

#include "../ensembleanalysismodule.h"

#ifdef VRN_MODULE_HDF5
#include "modules/hdf5/io/hdf5volumereader.h"
#endif

namespace voreen {

EnsembleVolumeReaderPopulator::EnsembleVolumeReaderPopulator(ProgressBar* progressBar)
    : volumeSerializerPopulator_(progressBar)
{
}

VolumeReader* EnsembleVolumeReaderPopulator::getVolumeReader(const std::string& path) const {

#ifdef VRN_MODULE_HDF5
    // For HDF5 files we first try to use the multi-channel reader.
    std::string ext = tgt::FileSystem::fileExtension(path);
    if (ext == "h5" || ext == "hdf5") {
        static HDF5VolumeReaderOriginal hdf5Reader;
        return &hdf5Reader;
    }
#endif

    // TODO: Add NIfTI multi-channel volume reader.

    try {
        return volumeSerializerPopulator_.getVolumeSerializer()->getReaders(path).front();
    }
    catch (tgt::UnsupportedFormatException&) {
    }

    return nullptr;
}


class VolumeRAMSwapDisk : public VolumeDisk {
public:

    VolumeRAMSwapDisk(Volume* volume)
        : VolumeDisk(volume->getFormat(), volume->getDimensions())
        , hash_(volume->getHash())
        , url_(volume->getOrigin())
        , owner_(volume)
    {
    }

    virtual std::string getHash() const {
        return hash_;
    }

    VolumeRAM* loadVolume() const {
        // The Disk -> RAM converter calls this function, hence we directly return the representation.
        LINFOC("VolumeRAMSwap", "Requesting entire volume...");
        return loadFromDisk();
    }

    VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
        return loadBrick(tgt::svec3(0, 0, firstZSlice),
                         tgt::svec3(getDimensions().xy(), lastZSlice-firstZSlice+1));
    }
    VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
        LINFOC("VolumeRAMSwap", "Requesting Brick...");

        // This function is typically called by derived data threads, after the checked for a VolumeRAM representation
        // to be present. Since the RAM representation was not present, but we don't want to read the entire volume
        // again and again, we add the representation here manually to the owner.

        if(!owner_->hasRepresentation<VolumeRAM>()) {
            owner_->addRepresentation(loadFromDisk());
        }

        VolumeRAMRepresentationLock ram(owner_);
        return ram->getSubVolume(dimensions, offset);
    }

protected:

    VolumeRAM* loadFromDisk() const {
        LINFOC("VolumeRAMSwap", "Reloading volume " << url_.getURL() << ". This may take a while...");

        EnsembleVolumeReaderPopulator populator;
        VolumeReader* reader = populator.getVolumeReader(url_.getPath());

        std::unique_ptr<Volume> volume(dynamic_cast<Volume*>(reader->read(url_)));
        tgtAssert(volume, "volume expected");
        tgtAssert(volume->getNumRepresentations() == 1, "Single RAM representation expected");

        VolumeRAMRepresentationLock lock(volume.get()); // Required to trigger memory manager to clean up.
        VolumeRAM* ramRepresentation = volume->getWritableRepresentation<VolumeRAM>();
        volume->releaseAllRepresentations(); // Give away ownership.
        return ramRepresentation;
    }

private:

    mutable Volume* owner_;
    std::string hash_;
    VolumeURL url_;
};


bool VolumeRAMSwap::tryAddVolumeSwap(VolumeBase* volumeBase) {
    EnsembleAnalysisModule* instance = EnsembleAnalysisModule::getInstance();
    if(!instance || !instance->getForceDiskRepresentation()) {
        return false;
    }

    if(!volumeBase || volumeBase->hasRepresentation<VolumeDisk>()) {
        return false;
    }

    Volume* volume = dynamic_cast<Volume*>(volumeBase);
    if(!volume || volume->getOrigin().getURL().empty()) {
        return false;
    }

    std::unique_ptr<VolumeRAMSwapDisk> swap;
    try {
        swap.reset(new VolumeRAMSwapDisk(volume));
    } catch(tgt::IOException&) {
        return false;
    }

    volume->addRepresentation(swap.release());
    volume->removeRepresentation<VolumeRAM>(); // TODO: Why is this necessary?

    LDEBUGC("VolumeRAMSwap", "Added Swap Disk for " << volume->getOrigin().getURL());
    return true;
}


}
