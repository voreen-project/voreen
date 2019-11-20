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

#include "voreen/core/datastructures/octree/octreebrickpoolmanagermmap.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/progressreporter.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "tgt/filesystem.h"

#include <vector>

namespace voreen {

static uint64_t numBricksInFile(uint64_t fileIndex) {
    return 1 << fileIndex;
}

OctreeBrickPoolManagerMmapStorageIndex::OctreeBrickPoolManagerMmapStorageIndex(uint64_t brickAddr) {
    tgtAssert(brickAddr != 0, "Invalid brick addr");

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
    // TODO
}

void OctreeBrickPoolManagerMmap::deinitialize() {
    // TODO
    flushPoolToDisk();

    OctreeBrickPoolManagerBase::deinitialize();
}

//-----------------------------------------------------------------------------------------------------------------------
//      DE-/SERIALIZATION
//-----------------------------------------------------------------------------------------------------------------------
void OctreeBrickPoolManagerMmap::serialize(Serializer& s) const {
    OctreeBrickPoolManagerBase::serialize(s);
    // TODO

    /*
    s.serialize("maxSingleBufferSizeBytes", maxBufferSizeBytes_);
    s.serialize("singleBufferSizeBytes",    singleBufferSizeBytes_);
    s.serialize("numBrickSlotsPerBuffer",   numBrickSlotsPerBuffer_);

    //save relative pathes
    std::vector<std::string> relBufferFiles;
    for(size_t i = 0; i < bufferFiles_.size(); i++)
        relBufferFiles.push_back(tgt::FileSystem::cleanupPath(tgt::FileSystem::relativePath(bufferFiles_[i],tgt::FileSystem::dirName(s.getDocumentPath())),false));
    s.serialize("bufferFiles", relBufferFiles);


    s.serialize("brickPoolPath", tgt::FileSystem::cleanupPath(tgt::FileSystem::relativePath(brickPoolPath_,tgt::FileSystem::dirName(s.getDocumentPath())),false));
    s.serialize("bufferFilePrefix", bufferFilePrefix_);

    s.serialize("nextVirtualMemoryAddress", nextVirtualMemoryAddress_);
    */
}

void  OctreeBrickPoolManagerMmap::deserialize(Deserializer& s) {
    OctreeBrickPoolManagerBase::deserialize(s);
    // TODO

    /*
    s.deserialize("maxSingleBufferSizeBytes", maxBufferSizeBytes_);
    s.deserialize("singleBufferSizeBytes",    singleBufferSizeBytes_);
    s.deserialize("numBrickSlotsPerBuffer", numBrickSlotsPerBuffer_);

    s.deserialize("bufferFiles", bufferFiles_);
    for(size_t i = 0; i < bufferFiles_.size(); i++)
        bufferFiles_[i] = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + bufferFiles_[i],true);

    s.deserialize("brickPoolPath", brickPoolPath_);
    brickPoolPath_ = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + brickPoolPath_,true);

    s.deserialize("bufferFilePrefix", bufferFilePrefix_);

    s.deserialize("nextVirtualMemoryAddress", nextVirtualMemoryAddress_);

    // check brick pool path
    if (!tgt::FileSystem::dirExists(brickPoolPath_))
        throw VoreenException("Brick pool path does not exist: " + brickPoolPath_);

    // make sure that buffer files are present
    if (bufferFiles_.empty())
        throw VoreenException("No brick buffer files");
    for (size_t i=0; i<bufferFiles_.size(); i++) {
        if (!tgt::FileSystem::fileExists(bufferFiles_.at(i)))
            throw VoreenException("Missing brick buffer file: " + bufferFiles_.at(i));
    }

    // check ram limit vs. buffer size
    if (2*singleBufferSizeBytes_ > ramLimitInBytes_)
        throw VoreenException("RAM memory limit is smaller than two times the size of a buffer. At least two buffer files have to fit in the RAM. "
            "[" + itos(ramLimitInBytes_) + " bytes < 2*" + itos(singleBufferSizeBytes_) + " bytes]");

    // init buffer vector
    maxNumBuffersInRAM_ = ramLimitInBytes_/singleBufferSizeBytes_;
    for (size_t i = 0; i < bufferFiles_.size(); i++)
        bufferVector_.push_back(new BufferEntry(numBrickSlotsPerBuffer_, 0, 0));
    */
}

//-----------------------------------------------------------------------------------------------------------------------
//      BRICK INTERACTION
//-----------------------------------------------------------------------------------------------------------------------
bool OctreeBrickPoolManagerMmap::isBrickInRAM(uint64_t virtualMemoryAddress) const {
    // Mmapped bricks are always in "RAM", as far as the program is concerned.
    return true;
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
                openParams.path = tgt::FileSystem::cleanupPath(brickPoolPath_ + "/" + bufferFilePrefix_ + std::to_string(newFileIndex));
                openParams.mode = std::ios::in | std::ios::out;
                openParams.length = fileSize;
                openParams.new_file_size = fileSize; // Option: Do create the file (overwrite if it does exist)!

                storage_.emplace_back(openParams);
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
    // Nothing to do here
}

uint64_t OctreeBrickPoolManagerMmap::allocateBrick() {
    return nextNewBrickAddr_.fetch_add(1); //TODO use free list
}

void OctreeBrickPoolManagerMmap::deleteBrick(uint64_t virtualMemoryAddress) {
    //TODO use free list
}

void OctreeBrickPoolManagerMmap::flushPoolToDisk(ProgressReporter* progressReporter /*= 0*/) {
    // TODO call fsync? See https://stackoverflow.com/questions/40972420/how-to-flush-memory-mapped-files-using-boosts-mapped-file-sink-class
}

//-----------------------------------------------------------------------------------------------------------------------
//      GENERAL FUNCTIONS
//-----------------------------------------------------------------------------------------------------------------------
uint64_t OctreeBrickPoolManagerMmap::getBrickPoolMemoryUsed() const {
    return 0; //TODO
}

std::string OctreeBrickPoolManagerMmap::getDescription() const {
    return "Brick Pool Manager that relies on the OS capabilities for brick caching.";
}

uint64_t OctreeBrickPoolManagerMmap::getBrickPoolMemoryAllocated() const {
    return 0; //TODO
}

} // namespace
