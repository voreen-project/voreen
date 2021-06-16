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

#include <memory>

#include "voreen/core/io/volumeserializer.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/datastructures/volume/volumeurl.h"

#include "../ensembleanalysismodule.h"

#ifdef VRN_MODULE_HDF5
#include "modules/hdf5/io/hdf5volumereader.h"
#endif

namespace voreen {

/**
 * If a volume reader (or file format) does not support a disk representation, a Swap disk can be added.
 * If a RAM representation is requested on a volume with a swap representation, it will be loaded from disk
 * on demand, which typically takes longer than just loading a disk representation, but serves as
 * (experimental) workaround for missing disk representations.
 */
class VolumeRAMSwap : public VolumeDisk {
public:

    static bool TryAddVolumeSwap(VolumeBase* volume) {
        // If the memory manager is not initialized, we would not benefit from a swap disk..
        if(!VolumeMemoryManager::isInited()) {
            return false;
        }

        if(!volume || volume->hasRepresentation<VolumeDisk>()) {
            return false;
        }

        if(volume->getOrigin().getURL().empty()) {
            return false;
        }

        std::unique_ptr<VolumeRAMSwap> swap;
        try {
            swap.reset(new VolumeRAMSwap(volume));
        } catch(tgt::IOException&) {
            return false;
        }

        volume->addRepresentation(swap.release());

        LDEBUG("Added Swap Disk for " << volume->getOrigin().getURL());
        return true;
    }


    VolumeRAMSwap(VolumeBase* volume)
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
        // The Disk -> RAM conversion calls this function, hence we directly return the representation.
        return loadFromDisk();
    }

    VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
        return loadBrick(tgt::svec3(0, 0, firstZSlice),
                         tgt::svec3(getDimensions().xy(), lastZSlice-firstZSlice+1));
    }
    VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {

        tgtAssert(VolumeMemoryManager::isInited(), "Volume Memory Manager not initialized");
        boost::lock_guard<boost::recursive_mutex> guard(*VolumeMemoryManager::getRef().getMutex());

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
        LINFO("Reloading volume " << url_.getURL() << ". This may take a while...");

        tgtAssert(VolumeMemoryManager::isInited(), "Volume Memory Manager not initialized");
        auto& vmm = VolumeMemoryManager::getRef();
        boost::lock_guard<boost::recursive_mutex> guard(*vmm.getMutex());

        // Try to free memory before loading into RAM.
        if(!vmm.requestMainMemory(owner_)) {
            // We will potentially run out of memory...
            LWARNING("Memory requirement not met, try to increase available CPU memory in the settings");
        }

        std::unique_ptr<Volume> volume(dynamic_cast<Volume*>(EnsembleVolumeReader().getVolumeReader(url_.getPath())->read(url_)));
        tgtAssert(volume, "volume expected");
        tgtAssert(volume->getNumRepresentations() == 1, "Single RAM representation expected");

        VolumeRAMRepresentationLock lock(volume.get());
        VolumeRAM* ramRepresentation = volume->getWritableRepresentation<VolumeRAM>();
        volume->releaseAllRepresentations(); // Give away ownership.
        return ramRepresentation;
    }

private:

    VolumeBase* owner_;
    std::string hash_;
    VolumeURL url_;

    static const std::string loggerCat_;
};

const std::string VolumeRAMSwap::loggerCat_ = "VolumeRAMSwap";


/////////////////////////////////////////////////////////////////////
// EnsembleVolumeReader
/////////////////////////////////////////////////////////////////////

EnsembleVolumeReader::EnsembleVolumeReader(ProgressBar* progressBar)
    : volumeSerializerPopulator_(progressBar)
{
    extensions_ = volumeSerializerPopulator_.getSupportedReadExtensions();
}

VolumeReader* EnsembleVolumeReader::create(ProgressBar* progressBar) const {
    return new EnsembleVolumeReader(progressBar);
}

bool EnsembleVolumeReader::canRead(const std::string& path) {
    return getVolumeReader(path) != nullptr;
}

VolumeReader* EnsembleVolumeReader::getVolumeReader(const std::string& path) const {

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

VolumeList* EnsembleVolumeReader::read(const std::string& url) {
    VolumeReader* reader = getVolumeReader(url);
    if(!reader) {
        return nullptr;
    }
    return reader->read(url);
}

std::vector<VolumeURL> EnsembleVolumeReader::listVolumes(const std::string& url) const {
    VolumeReader* reader = getVolumeReader(url);
    if(!reader) {
        return std::vector<VolumeURL>();
    }
    return reader->listVolumes(url);
}

VolumeBase* EnsembleVolumeReader::read(const VolumeURL& url) {
    VolumeReader* reader = getVolumeReader(url.getPath());
    if (!reader) {
        LERRORC("voreen.EnsembleVolumeReader", "No suitable reader found for " << url.getURL());
        return nullptr;
    }

    VolumeBase* volume = nullptr;
    try {
        volume = reader->read(url);
    } catch (tgt::FileException& e) {
    }

    if (!volume) {
        LERRORC("voreen.EnsembleVolumeReader", "Could not read volume " << url.getURL());
        return nullptr;
    }

    // Add swap disk, if required.
    EnsembleAnalysisModule* instance = EnsembleAnalysisModule::getInstance();
    if(instance && instance->getForceDiskRepresentation()) {
        bool success = VolumeRAMSwap::TryAddVolumeSwap(volume);
        if(success) {

            // Enforce derived data calculation right before we add a swap disk.
            // The reason for this is that otherwise the volume has possibly to be loaded twice when added to an
            // ensemble at a later point, because we delete the RAM representation below.
            // The latter we need to do, because loading a volume with a volume reader will not request
            // memory from the memory manager.
            LDEBUGC("voreen.EnsembleVolumeReader", "Calculating derived data...");
            volume->getDerivedData<VolumeMinMax>();
            if(volume->getNumChannels() > 1) {
                volume->getDerivedData<VolumeMinMaxMagnitude>();
            }

            volume->removeRepresentation<VolumeRAM>();
        }
    }

    return volume;
}


}
