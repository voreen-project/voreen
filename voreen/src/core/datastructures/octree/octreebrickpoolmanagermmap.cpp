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

#include "voreen/core/datastructures/octree/octreebrickpoolmanagermmap.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/progressreporter.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "tgt/filesystem.h"

#include <vector>

namespace voreen {

static inline uint64_t numBricksInFile(uint64_t fileIndex) {
    return 1 << fileIndex;
}
static inline std::array<uint16_t, 4> packAddr(uint64_t addr) {
    return { static_cast<uint16_t>(addr), static_cast<uint16_t>(addr >> 16), static_cast<uint16_t>(addr >> 32), static_cast<uint16_t>(addr >> 48) };
}
static inline uint64_t unpackAddr(std::array<uint16_t, 4> packed) {
    return static_cast<uint64_t>(packed[0]) | static_cast<uint64_t>(packed[1]) << 16 | static_cast<uint64_t>(packed[2]) << 32 | static_cast<uint64_t>(packed[3]) << 48;
}

OctreeBrickPoolManagerMmapStorageIndex::OctreeBrickPoolManagerMmapStorageIndex(uint64_t brickAddr) {
    tgtAssert(brickAddr != 0, "Invalid brick addr");
    tgtAssert(brickAddr != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS, "Invalid brick addr");

    fileIndex_ = 0;
    while((brickAddr >> (fileIndex_+1)) != 0) {
        fileIndex_++;
    }
    brickIndexInFile_ = brickAddr - numBricksInFile(fileIndex_);
}


OctreeBrickPoolManagerMmap::OctreeBrickPoolManagerMmap(const std::string& brickPoolPath, const std::string& bufferFilePrefix)
    : OctreeBrickPoolManagerBase()
    , brickPoolPath_(brickPoolPath)
    , bufferFilePrefix_(bufferFilePrefix)
    , storage_()
    , freeListHead_(NO_BRICK_ADDRESS)
    , freeListMutex_()
    , nextNewBrickAddr_(1)
{
    brickPoolPath_ = tgt::FileSystem::absolutePath(brickPoolPath);
    bufferFilePrefix_ = (!bufferFilePrefix.empty() ? bufferFilePrefix : "brickbuffer_");
}

OctreeBrickPoolManagerMmap::~OctreeBrickPoolManagerMmap() {
}

OctreeBrickPoolManagerMmap* OctreeBrickPoolManagerMmap::create() const {
    return new OctreeBrickPoolManagerMmap("", "");
}

//-----------------------------------------------------------------------------------------------------------------------
//      DE-/INITIALIZE
//-----------------------------------------------------------------------------------------------------------------------
void OctreeBrickPoolManagerMmap::initialize(size_t brickMemorySizeInByte) {
    OctreeBrickPoolManagerBase::initialize(brickMemorySizeInByte);
}

void OctreeBrickPoolManagerMmap::deinitialize() {
    flushPoolToDisk();

    OctreeBrickPoolManagerBase::deinitialize();
}

//-----------------------------------------------------------------------------------------------------------------------
//      DE-/SERIALIZATION
//-----------------------------------------------------------------------------------------------------------------------
void OctreeBrickPoolManagerMmap::serialize(Serializer& s) const {
    OctreeBrickPoolManagerBase::serialize(s);

    s.serialize("brickPoolPath", tgt::FileSystem::cleanupPath(tgt::FileSystem::relativePath(brickPoolPath_,tgt::FileSystem::dirName(s.getDocumentPath())), false));
    s.serialize("bufferFilePrefix", bufferFilePrefix_);

    size_t numFiles;
    {
        boost::shared_lock<boost::shared_mutex> rlock(storageMutex_);
        numFiles = storage_.size();
    }

    std::vector<std::string> relBufferFiles;
    for(uint64_t i=0; i<numFiles; ++i) {
        relBufferFiles.push_back(tgt::FileSystem::cleanupPath(tgt::FileSystem::relativePath(storageFileName(i),tgt::FileSystem::dirName(s.getDocumentPath())),false));
    }
    s.serialize("bufferFiles", relBufferFiles);

    uint64_t nextNewBrickAddr = nextNewBrickAddr_;
    s.serialize("nextNewBrickAddr", nextNewBrickAddr);

    uint64_t freeListHead = freeListHead_;
    s.serialize("freeListHead", freeListHead);
}

void  OctreeBrickPoolManagerMmap::deserialize(Deserializer& s) {
    OctreeBrickPoolManagerBase::deserialize(s);

    s.deserialize("brickPoolPath", brickPoolPath_);
    brickPoolPath_ = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + brickPoolPath_, true);

    s.deserialize("bufferFilePrefix", bufferFilePrefix_);

    std::vector<std::string> relBufferFiles;
    s.deserialize("bufferFiles", relBufferFiles);
    {
        boost::unique_lock<boost::shared_mutex> rlock(storageMutex_);

        for(uint64_t i=0; i<relBufferFiles.size(); ++i) {
            size_t fileSize = numBricksInFile(i) * getBrickMemorySizeInByte();

            std::string filename = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + relBufferFiles[i],true);

            boost::iostreams::mapped_file_params openParams;
            openParams.path = filename;
            openParams.mode = std::ios::in | std::ios::out;
            openParams.length = fileSize;

            storage_.emplace_back(openParams);
        }
    }
    uint64_t nextNewBrickAddr;
    s.deserialize("nextNewBrickAddr", nextNewBrickAddr);
    nextNewBrickAddr_ = nextNewBrickAddr;

    uint64_t freeListHead;
    s.deserialize("freeListHead", freeListHead);
    freeListHead_ = freeListHead;
}

std::string OctreeBrickPoolManagerMmap::storageFileName(uint64_t fileIndex) const {
    //TODO revisit if naming scheme changes e.g. when saving with octreesave (it shouldn't).
    return tgt::FileSystem::cleanupPath(brickPoolPath_ + "/" + bufferFilePrefix_ + std::to_string(fileIndex));
}

//-----------------------------------------------------------------------------------------------------------------------
//      BRICK INTERACTION
//-----------------------------------------------------------------------------------------------------------------------
bool OctreeBrickPoolManagerMmap::isBrickInRAM(uint64_t virtualMemoryAddress) const {
    // We have not idea if a brick is in RAM or not, so we pessimistically
    // assume "no". This is important for timeouts in octree rendering/texture
    // compositon to work.
    return false;
}

const uint16_t* OctreeBrickPoolManagerMmap::getBrick(uint64_t virtualMemoryAddress, bool blocking) const {
    // For mmap, there is no difference between normal and writable bricks
    return getWritableBrick(virtualMemoryAddress, blocking);
}

uint16_t* OctreeBrickPoolManagerMmap::getWritableBrick(uint64_t virtualMemoryAddress, bool blocking) const {
    OctreeBrickPoolManagerMmapStorageIndex storageIndex(virtualMemoryAddress);

    {
        boost::upgrade_lock<boost::shared_mutex> rlock(storageMutex_);

        if(storageIndex.fileIndex_ >= storage_.size()) {
            boost::upgrade_to_unique_lock<boost::shared_mutex> wlock(rlock);

            // Need to allocate new storage block
            while(storageIndex.fileIndex_ >= storage_.size()) {
                uint64_t newFileIndex = storage_.size();

                size_t fileSize = numBricksInFile(newFileIndex) * getBrickMemorySizeInByte();

                boost::iostreams::mapped_file_params openParams;
                openParams.path = storageFileName(newFileIndex);
                openParams.mode = std::ios::in | std::ios::out;
                openParams.length = fileSize;
                openParams.new_file_size = fileSize; // Option: Do create the file (overwrite if it does exist)!

                boost::iostreams::mapped_file file(openParams);
                if(!file.is_open()) {
                    throw VoreenException("Unable to create brick file (possibly no space left on device?)");
                }

                storage_.emplace_back(file);
            }
        }

        tgtAssert(storageIndex.fileIndex_ < storage_.size(), "Somehow failed to allocate storage block");

        auto& file = storage_[storageIndex.fileIndex_];

        size_t brickByteOffset = storageIndex.brickIndexInFile_ * getBrickMemorySizeInByte();
        tgtAssert(brickByteOffset < file.size(), "Somehow failed to allocate storage block");

        char* brickPtr = file.begin() + brickByteOffset;

        return reinterpret_cast<uint16_t*>(brickPtr);
    }
}

void OctreeBrickPoolManagerMmap::releaseBrick(uint64_t virtualMemoryAddress, AccessMode mode) const {
    // Nothing to do here: Blocks are writting back to disk by OS.
}

uint64_t OctreeBrickPoolManagerMmap::allocateBrick() {

    // Fast path: Check if free list is empty. If this is the case we just
    // allocate a new brick. It doesn't matter if we miss an update to the free
    // list here: Some later allocation will get the brick from the list.
    if(freeListHead_.load() == NO_BRICK_ADDRESS)  {
        // No blocks in free list: allocate completely new block
        return nextNewBrickAddr_.fetch_add(1);
    }

    // Ok, we now know that we probably need to modify the free list.
    boost::unique_lock<boost::mutex> wlock(freeListMutex_);

    // Check again, now that we have exclusive access (there may have been
    // concurrent accesses in the meantime!)
    uint64_t currentListHeadAddr = freeListHead_.load();
    if(currentListHeadAddr == NO_BRICK_ADDRESS)  {
        // No blocks in free list: allocate completely new block
        return nextNewBrickAddr_.fetch_add(1);
    }
    // Now that we know that there is definitely a brick that the freeListHead_
    // points to, we can remove it from the free list:

    // Read addr block after head in list from header block itself
    const uint16_t* currentListHead = getBrick(currentListHeadAddr);
    std::array<uint16_t, 4> newListHeadAddrPacked;
    std::copy_n(currentListHead, newListHeadAddrPacked.size(), newListHeadAddrPacked.begin());

    // This next element in the free list is the new head, and we return the
    // old head (now removed from the list)
    freeListHead_.store(unpackAddr(newListHeadAddrPacked));
    return currentListHeadAddr;
    //LINFO("Reusing block from free list" << currentListHeadAddr);
}

void OctreeBrickPoolManagerMmap::deleteBrick(uint64_t virtualMemoryAddress) {
    tgtAssert(getBrickMemorySizeInByte() >= sizeof(uint64_t), "Bricks too small to host free list");
    //LINFO("Deleting block" << virtualMemoryAddress);

    uint16_t* newListHead = getWritableBrick(virtualMemoryAddress);

    boost::unique_lock<boost::mutex> wlock(freeListMutex_);

    // Save the current list header addr as data in the deleted block
    auto currentListHeadAddrPacked = packAddr(freeListHead_.load());
    std::copy(currentListHeadAddrPacked.begin(), currentListHeadAddrPacked.end(), newListHead);

    // Make the deleted block address the new free list header
    freeListHead_.store(virtualMemoryAddress);
}

void OctreeBrickPoolManagerMmap::flushPoolToDisk(ProgressReporter* progressReporter /*= 0*/) {
    // Bricks are automatically written back to disk.
    // It may be possible to call fsync or similar, but it is unclear if this is in line with the semantics of this function.
}

//-----------------------------------------------------------------------------------------------------------------------
//      GENERAL FUNCTIONS
//-----------------------------------------------------------------------------------------------------------------------
uint64_t OctreeBrickPoolManagerMmap::getBrickPoolMemoryUsed() const {
    return 0; // No way to query this. No memory is explicitly allocated.
}

std::string OctreeBrickPoolManagerMmap::getDescription() const {
    return "Brick Pool Manager that relies on the OS capabilities for brick caching.";
}
const std::string& OctreeBrickPoolManagerMmap::getBrickPoolPath() const {
    return brickPoolPath_;
}

uint64_t OctreeBrickPoolManagerMmap::getBrickPoolMemoryAllocated() const {
    return 0; // No way to query this. No memory is explicitly allocated.
}

} // namespace
