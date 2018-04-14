/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "octreebrickpoolmanagerdisksinglethreaded.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/serialization/serialization.h"

#include <time.h>
#include <fstream>

#include "tgt/assert.h"
#include "tgt/logmanager.h"
#include "tgt/tgt_math.h"
#include "tgt/filesystem.h"
#include "tgt/vector.h"

namespace voreen {

const std::string OctreeBrickPoolManagerDiskSingleThreaded::loggerCat_("voreen.OctreeBrickPoolManagerDisk");

OctreeBrickPoolManagerDiskSingleThreaded::OctreeBrickPoolManagerDiskSingleThreaded(const size_t maxBufferSize,
            const std::string& brickPoolPath, const std::string& bufferFilePrefix)
    : maxBufferSizeBytes_(maxBufferSize)
    , brickPoolPath_(brickPoolPath)
    , bufferFilePrefix_(bufferFilePrefix)
{
    brickPoolPath_ = tgt::FileSystem::absolutePath(brickPoolPath);
    bufferFilePrefix_ = (!bufferFilePrefix.empty() ? bufferFilePrefix : "brickbuffer");
}

OctreeBrickPoolManagerDiskSingleThreaded::~OctreeBrickPoolManagerDiskSingleThreaded() {
}

size_t OctreeBrickPoolManagerDiskSingleThreaded::getNumBrickBuffers() const {
    return bufferFiles_.size();
}

std::vector<std::string> OctreeBrickPoolManagerDiskSingleThreaded::getBufferFiles() const {
    return bufferFiles_;
}

void OctreeBrickPoolManagerDiskSingleThreaded::initialize(size_t brickMemorySizeInByte) {
    OctreeBrickPoolManagerBase::initialize(brickMemorySizeInByte);

    // round max buffer size down to next multiple of brick memory size
    if (maxBufferSizeBytes_ < getBrickMemorySizeInByte())
        throw VoreenException("Max brick buffer size is smaller than the memory size of a single brick "
                              "[" + itos(maxBufferSizeBytes_) + " bytes < " + itos(getBrickMemorySizeInByte()) + " bytes]");
    singleBufferSizeBytes_ = maxBufferSizeBytes_;
    if (!isMultipleOf(singleBufferSizeBytes_, getBrickMemorySizeInByte()))
        singleBufferSizeBytes_ = tgt::ifloor((float)maxBufferSizeBytes_ / (float)getBrickMemorySizeInByte()) * getBrickMemorySizeInByte();

    numBrickSlotsPerBuffer_ = singleBufferSizeBytes_ / brickMemorySizeInByte;

    if (!tgt::FileSystem::dirExists(brickPoolPath_))
        throw VoreenException("Brick pool path does not exist: " + brickPoolPath_);
}

void OctreeBrickPoolManagerDiskSingleThreaded::deinitialize() {
    bufferFiles_.clear();
    OctreeBrickPoolManagerBase::deinitialize();
}

char* OctreeBrickPoolManagerDiskSingleThreaded::allocateNextDiskBuffer() {

    LDEBUG("Allocating brick buffer " << bufferFiles_.size() << " (" << formatMemorySize(getBrickBufferSizeInBytes()) << ")");

    char* buffer = 0;
    try {
        buffer = new char[getBrickBufferSizeInBytes()];
    }
    catch(std::bad_alloc& e) {
        throw VoreenException("Failed to allocate brick buffer: " + std::string(e.what()));
    }

    std::stringstream path;
    path << brickPoolPath_ << "/" << bufferFilePrefix_ << itos(getNumBrickBuffers(), 10) << ".raw";

    std::ofstream outfile(path.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (outfile.fail())
        throw VoreenException("Could not open buffer file: " + path.str());
    outfile.write(buffer, getBrickBufferSizeInBytes());
    outfile.close();

    /* C-style I/O
    FILE* f = fopen(path.str().c_str(),"wb");
    setbuf(f,NULL);
    fwrite(buffer,1,getBrickBufferSizeInBytes(),f);
    fflush(f);
    fclose(f); */

    bufferFiles_.push_back(path.str());
    return buffer;
}

void OctreeBrickPoolManagerDiskSingleThreaded::saveBufferToDisk(size_t bufferIndex, const char* buffer) const {
    tgtAssert(bufferIndex < getNumBrickBuffers(), "invalid buffer index");
    tgtAssert(buffer, "null pointer passed");

    LDEBUG("Writing buffer " << bufferIndex << " to disk");

    std::ofstream outfile(bufferFiles_[bufferIndex].c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (outfile.fail())
        throw VoreenException("Could not open buffer file for writing: " + bufferFiles_[bufferIndex]);
    outfile.write(buffer, getBrickBufferSizeInBytes());
    outfile.close();

    /* C-style I/O
    FILE* f = fopen(bufferFiles_[bufferIndex].c_str(),"wb");
    setbuf(f,NULL);
    fwrite(buffer,1,getBrickBufferSizeInBytes(),f);
    fflush(f);
    fclose(f); */
}

char* OctreeBrickPoolManagerDiskSingleThreaded::loadDiskBuffer(size_t bufferIndex) const {
    tgtAssert(bufferIndex < getNumBrickBuffers(), "invalid buffer index");
    tgtAssert(getBufferFiles().size() == getNumBrickBuffers(), "number of buffer files does not match number of numbers");

    LDEBUG("Loading buffer " << bufferIndex << " from disk");

    const std::string bufferFile = getBufferFiles().at(bufferIndex);
    tgtAssert(!bufferFile.empty(), "buffer file path is empty");
    if (!tgt::FileSystem::fileExists(bufferFile))
        throw VoreenException("Buffer file does not exist: " + bufferFile);

    char* buffer = 0;
    try {
        buffer = new char[getBrickBufferSizeInBytes()];
    }
    catch(std::bad_alloc& e) {
        throw VoreenException("Failed to load disk buffer into RAM: " + std::string(e.what()));
    }

    std::ifstream infile(bufferFile.c_str(), std::ios::in | std::ios::binary);
    if(infile.fail())
        throw VoreenException("Could not open: " + bufferFile);
    infile.read(buffer, getBrickBufferSizeInBytes());
    infile.close();

    /* C-style I/O
    FILE* f = fopen(bufferFile.c_str(),"rb");
    setbuf(f,NULL);
    fread(buffer,1,getBrickBufferSizeInBytes(),f);
    fclose(f); */

    return buffer;
}

size_t OctreeBrickPoolManagerDiskSingleThreaded::getBrickBufferSizeInBytes() const {
    return singleBufferSizeBytes_;
}

size_t OctreeBrickPoolManagerDiskSingleThreaded::getNumBrickSlotsPerBuffer() const {
    return numBrickSlotsPerBuffer_;
}

void OctreeBrickPoolManagerDiskSingleThreaded::serialize(Serializer& s) const {
    OctreeBrickPoolManagerBase::serialize(s);

    s.serialize("maxSingleBufferSizeBytes", maxBufferSizeBytes_);
    s.serialize("singleBufferSizeBytes",    singleBufferSizeBytes_);

    s.serialize("numBrickSlotsPerBuffer", numBrickSlotsPerBuffer_);

    s.serialize("bufferFiles", bufferFiles_);
    s.serialize("brickPoolPath", brickPoolPath_);
    s.serialize("bufferFilePrefix", bufferFilePrefix_);
}

void OctreeBrickPoolManagerDiskSingleThreaded::deserialize(Deserializer& s) {
    OctreeBrickPoolManagerBase::deserialize(s);

    s.deserialize("maxSingleBufferSizeBytes", maxBufferSizeBytes_);
    s.deserialize("singleBufferSizeBytes",    singleBufferSizeBytes_);

    s.deserialize("numBrickSlotsPerBuffer", numBrickSlotsPerBuffer_);

    s.deserialize("bufferFiles", bufferFiles_);
    s.deserialize("brickPoolPath", brickPoolPath_);
    s.deserialize("bufferFilePrefix", bufferFilePrefix_);

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

}


//-------------------------------------------------------------------------------------------------
// OctreeBrickPoolManagerDiskLimitedRam

const std::string OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::loggerCat_("voreen.OctreeBrickPoolManagerDiskLimitedRam");

OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::OctreeBrickPoolManagerDiskSingleThreadedLimitedRam(const size_t maxSingleBufferSize, const size_t ramLimitInBytes,
                                         const std::string& brickPoolPath, const std::string& bufferFilePrefix, const AllocateStrategy strategy)
    : OctreeBrickPoolManagerDiskSingleThreaded(maxSingleBufferSize, brickPoolPath, bufferFilePrefix)
    , ramLimitInBytes_(ramLimitInBytes)
    , nextVirtualMemoryAddress_(0)
    , nextTimestemp_(0)
    , allocateStrategy_(strategy)
{
    tgtAssert(ramLimitInBytes_ > 0, "RAM limit must be larger than 0");
}

OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::~OctreeBrickPoolManagerDiskSingleThreadedLimitedRam() {
    brickBuffersInRAM_.clear();
}

OctreeBrickPoolManagerDiskSingleThreadedLimitedRam* OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::create() const {
    return new OctreeBrickPoolManagerDiskSingleThreadedLimitedRam(64 << 20, 512 << 20, "", "");
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::setRAMLimit(size_t ramLimitInBytes) {
    tgtAssert(ramLimitInBytes > 0, "RAM limit must be larger than 0");

    // to work properly, at least 2 buffers have to fit into the ram
    if (2*getBrickMemorySizeInByte() > ramLimitInBytes)
        throw VoreenException("RAM memory limit is smaller than two times the size of a buffer. At least two buffer files have to fit in the RAM. "
                              "[" + itos(ramLimitInBytes) + " bytes < 2*" + itos(getBrickMemorySizeInByte()) + " bytes]");

    if (ramLimitInBytes != ramLimitInBytes_) {
        ramLimitInBytes_ = ramLimitInBytes;

        // re-init RAM buffer vectors
        flushPoolToDisk();
        deletedBricks_.clear();
        bricksInUsage_.clear();
        brickBuffersInRAM_ = std::vector<BrickBufferInRAM>(ramLimitInBytes_ / getBrickBufferSizeInBytes(), BrickBufferInRAM());
        bricksOfBuffersInUsage_ = std::vector<size_t>(ramLimitInBytes_ / getBrickBufferSizeInBytes(), 0);
    }
}

size_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::getRAMLimit() const {
    return ramLimitInBytes_;
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::initialize(size_t brickMemorySizeInByte) {
    OctreeBrickPoolManagerDiskSingleThreaded::initialize(brickMemorySizeInByte);

    // to work properly, at least 2 buffers have to fit into the ram
    if (2*getBrickBufferSizeInBytes() > ramLimitInBytes_)
        throw VoreenException("RAM memory limit is smaller than two times the size of a buffer. At least two buffer files have to fit in the RAM. "
                              "[" + itos(ramLimitInBytes_) + " bytes < 2*" + itos(getBrickMemorySizeInByte()) + " bytes]");

    // set brick vector of buffers in the ram
    brickBuffersInRAM_ = std::vector<BrickBufferInRAM>(ramLimitInBytes_ / getBrickBufferSizeInBytes(), BrickBufferInRAM());
    bricksOfBuffersInUsage_ = std::vector<size_t>(ramLimitInBytes_ / getBrickBufferSizeInBytes(), 0);
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::deinitialize() {
    flushPoolToDisk();

    brickBuffersInRAM_.clear();
    deletedBricks_.clear();
    bricksInUsage_.clear();
    bricksOfBuffersInUsage_.clear();

    nextVirtualMemoryAddress_ = 0;

    OctreeBrickPoolManagerDiskSingleThreaded::deinitialize();
}

bool OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::isBrickInRAM(uint64_t virtualMemoryAddress) const {
    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / getBrickBufferSizeInBytes());
    for(size_t i = 0; i < brickBuffersInRAM_.size(); i++) {
        if(brickBuffersInRAM_[i].index_ == bufferID)
            return true;
    }
    return false;
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::flushPoolToDisk(ProgressReporter* progressReporter /*= 0*/) {
    try {
        for (size_t i=0; i<brickBuffersInRAM_.size(); i++) {
            if (brickBuffersInRAM_[i].writable_) {
                saveBufferToDisk(brickBuffersInRAM_[i].index_,brickBuffersInRAM_[i].pointer_);
                if (!brickBuffersInRAM_[i].inUsage_)
                    brickBuffersInRAM_[i].writable_ = false;
            }
        }
    }
    catch (VoreenException& e) {
        LERROR(e.what());
    }
}

size_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::isBufferInRAM(size_t bufferID) const {
    for(size_t i = 0; i < brickBuffersInRAM_.size(); i++) {
        if(brickBuffersInRAM_[i].index_ == bufferID)
            return i;
    }
    return (size_t)(-1);
}

uint64_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::allocateBrick() {
    switch(allocateStrategy_) {
    case IGNORE_DELETED_BRICKS:
        return ignoreFreeBricks();
        break;
    case USE_DELETED_BRICKS:
        return useFreeBricks();
        break;
    default:
        //shouldn't get here
        tgtAssert(false,"Shouldn't get here!!!");
        return 0;
        break;
    }
    //shouldn't get here
    tgtAssert(false,"Shouldn't get here!!!");
    return 0;
}

uint64_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::ignoreFreeBricks() {
    // TODO: checks
    size_t bufferID = static_cast<size_t>(nextVirtualMemoryAddress_ / getBrickBufferSizeInBytes());
    while (bufferID >= getBufferFiles().size()) {
    //load brick from disk
        BrickBufferInRAM* leastUsedBuffer = 0;
        for(size_t i = 0; i < brickBuffersInRAM_.size(); i++){
            if(!brickBuffersInRAM_[i].inUsage_  && (leastUsedBuffer == 0 || leastUsedBuffer->lastUsed_ > brickBuffersInRAM_[i].lastUsed_)) {
                leastUsedBuffer = const_cast<BrickBufferInRAM*>(&brickBuffersInRAM_[i]);
            }
        }
        //check, if all buffers are in use
        if(!leastUsedBuffer) {
            throw VoreenException("All buffers in RAM are in use!!!");
        }
        //save buffer if it is writable
        if(leastUsedBuffer->writable_)
            saveBufferToDisk(leastUsedBuffer->index_,leastUsedBuffer->pointer_);
        //LERROR("IGNORE kicked: "<< leastUsedBuffer->index_ << " loaded: " << getBufferFiles().size());
        //clean struct and return
        delete[] leastUsedBuffer->pointer_;
        leastUsedBuffer->index_ = getBufferFiles().size();
        leastUsedBuffer->pointer_ = allocateNextDiskBuffer(); //no bad_alloc, since memory has been freed
        leastUsedBuffer->lastUsed_ = increaseTimestemp();
        //LERROR("Loaded Buffer ID: " << getBufferFiles().size()-1); //TODO: DEBUG
    }
    //return value
    uint64_t returnValue = nextVirtualMemoryAddress_;
    nextVirtualMemoryAddress_ += static_cast<uint64_t>(getBrickMemorySizeInByte());
    return returnValue;
}

uint64_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::useFreeBricks() {
    //case1: actual buffer is not full -> use it
    //case2: we have free bricks -> use them
    //case3: we have to allocate a new buffer

    //case 1
    if (nextVirtualMemoryAddress_%getBrickBufferSizeInBytes() != 0) {
        uint64_t returnValue = nextVirtualMemoryAddress_;
        nextVirtualMemoryAddress_ += static_cast<uint64_t>(getBrickMemorySizeInByte());
        return returnValue;
    } else //case2
    if (!deletedBricks_.empty()) {
        uint64_t returnValue = deletedBricks_.back();
        deletedBricks_.pop_back();
        return returnValue;
    } else { //case3
        return ignoreFreeBricks();
    }
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::deleteBrick(uint64_t virtualMemoryAddress) {
    deletedBricks_.push_back(virtualMemoryAddress);
    releaseBrick(virtualMemoryAddress);
}

const uint16_t* OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::getBrick(uint64_t virtualMemoryAddress, bool blocking) const {
    //constraint SIZE_MAX equals no brick
    if(virtualMemoryAddress == (size_t)(-1)) return 0;
    //transform virtual memory address
    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / getBrickBufferSizeInBytes());
    size_t bufferOffset = virtualMemoryAddress % getBrickBufferSizeInBytes();
    //load/get buffer
    size_t bufferIndex = loadBrickBufferIntoRAM(bufferID);
    //buffer exists in ram
    if(bufferIndex < (size_t)(-1)) {
        setBrickUsed(virtualMemoryAddress, bufferIndex);
        return reinterpret_cast<uint16_t*>(brickBuffersInRAM_[bufferIndex].pointer_ + bufferOffset);
    } else //buffer could not be laoded
        return 0;
}

uint16_t* OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::getWritableBrick(uint64_t virtualMemoryAddress, bool blocking) const {
    //constraint SIZE_MAX equals no brick
    if(virtualMemoryAddress == (size_t)(-1)) return 0;
    //transform virtual memory address
    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / getBrickBufferSizeInBytes());
    size_t bufferOffset = virtualMemoryAddress % getBrickBufferSizeInBytes();
    //load/get buffer
    size_t bufferIndex = loadBrickBufferIntoRAM(bufferID);
    //buffer exists in ram
    if(bufferIndex < (size_t)(-1)) {
        brickBuffersInRAM_[bufferIndex].writable_ = true;
        setBrickUsed(virtualMemoryAddress, bufferIndex);
        return reinterpret_cast<uint16_t*>(brickBuffersInRAM_[bufferIndex].pointer_ + bufferOffset);
    } else //buffer could not be laoded
        return 0;
}

size_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::loadBrickBufferIntoRAM(size_t bufferID) const {
    // TODO: add checks
    if (bufferID >= getBufferFiles().size()) {
        tgtAssert(false,"Buffer has not been created");
        LERROR("loadBrickBufferIntoRAM(): Buffer has not been created");
        return (size_t)(-1);
    }
    else {
        size_t bufferIndex = isBufferInRAM(bufferID);
        if(bufferIndex < (size_t)(-1)) {
            //increase timestemp and return index
            brickBuffersInRAM_[bufferIndex].lastUsed_ = increaseTimestemp();
            return bufferIndex;
        } else {
            //load brick from disk
            BrickBufferInRAM* leastUsedBuffer = 0;
            for(size_t i = 0; i < brickBuffersInRAM_.size(); i++){
                if(!brickBuffersInRAM_[i].inUsage_  && (leastUsedBuffer == 0 || leastUsedBuffer->lastUsed_ > brickBuffersInRAM_[i].lastUsed_)) {
                    leastUsedBuffer = const_cast<BrickBufferInRAM*>(&brickBuffersInRAM_[i]);
                    bufferIndex = i;
                }
            }
            //check, if all buffers are in use
            if(!leastUsedBuffer) {
                throw VoreenException("All buffers in RAM are in use!!!");
            }
            //save buffer if it is writable
            if(leastUsedBuffer->writable_)
                saveBufferToDisk(leastUsedBuffer->index_,leastUsedBuffer->pointer_);
            //LERROR("LOAD kicked: "<< leastUsedBuffer->index_ << " loaded: " << bufferID);
            //clean struct and return
            delete[] leastUsedBuffer->pointer_;
            leastUsedBuffer->pointer_ = loadDiskBuffer(bufferID); //no bad_alloc, since memory has been freed
            leastUsedBuffer->index_ = bufferID;
            leastUsedBuffer->lastUsed_ = increaseTimestemp();
            leastUsedBuffer->writable_ = false;
            //LERROR("Loaded Buffer ID: " << bufferID); //TODO: DEBUG
            return bufferIndex;
        }
    }
    //shouldn't get here
    tgtAssert(false,"Shouldn't get here!!!")
    return (size_t)(-1);
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::releaseBrick(uint64_t virtualMemoryAddress, AccessMode mode) const {
    //remove address
    if(bricksInUsage_.erase(virtualMemoryAddress)){
        size_t bufferIndex = isBufferInRAM(static_cast<size_t>(virtualMemoryAddress / static_cast<uint64_t>(getBrickBufferSizeInBytes())));
        tgtAssert(bufferIndex < (size_t)(-1), "Used Buffer was not in RAM!!!");
        tgtAssert(bricksOfBuffersInUsage_[bufferIndex] > 0, "Buffer has less used bricks!!!")
        bricksOfBuffersInUsage_[bufferIndex] -= 1;
        if(bricksOfBuffersInUsage_[bufferIndex] == 0)
            brickBuffersInRAM_[bufferIndex].inUsage_ = false;
    }
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::setBrickUsed(uint64_t virtualMemoryAddress, size_t bufferIndex) const {
    //add address
    if(bricksInUsage_.insert(virtualMemoryAddress).second){
        bricksOfBuffersInUsage_[bufferIndex] += 1;
        if(bricksOfBuffersInUsage_[bufferIndex] == 1)
            brickBuffersInRAM_[bufferIndex].inUsage_ = true;
    }
}

size_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::increaseTimestemp() const {
    nextTimestemp_++;

    //handle overflow
    if(nextTimestemp_ == (size_t)(-1)){
        //find minimum
        size_t minV = -1;
        for(size_t i = 0; i < brickBuffersInRAM_.size(); ++i){
            if(brickBuffersInRAM_[i].lastUsed_ < minV)
                minV = brickBuffersInRAM_[i].lastUsed_;
        }
        //subtract min and find max
        size_t maxV = 0;
        for(size_t i = 0; i < brickBuffersInRAM_.size(); ++i){
            brickBuffersInRAM_[i].lastUsed_ -= minV;
            if(brickBuffersInRAM_[i].lastUsed_ > maxV)
                maxV = brickBuffersInRAM_[i].lastUsed_;
        }
        nextTimestemp_ = maxV+1;
    }
    return nextTimestemp_;
}

uint64_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::getBrickPoolMemoryAllocated() const {
    return static_cast<uint64_t>(getNumBrickBuffers())*static_cast<uint64_t>(getBrickBufferSizeInBytes());
}

uint64_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::getBrickPoolMemoryUsed() const {
    // TODO
    return 0;
}

size_t OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::getNumBrickBuffersInRAM() const {
    size_t numBuffers = 0;
    for (size_t i=0; i<brickBuffersInRAM_.size(); i++)
        if (brickBuffersInRAM_.at(i).pointer_)
            numBuffers++;
    return numBuffers;
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::serialize(Serializer& s) const {
    OctreeBrickPoolManagerDiskSingleThreaded::serialize(s);
    s.serialize("nextVirtualMemoryAddress", nextVirtualMemoryAddress_);
}

void OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::deserialize(Deserializer& s) {
    OctreeBrickPoolManagerDiskSingleThreaded::deserialize(s);

    s.deserialize("nextVirtualMemoryAddress", nextVirtualMemoryAddress_);

    //set brick vector of buffers in the ram
    brickBuffersInRAM_ = std::vector<BrickBufferInRAM>(ramLimitInBytes_/getBrickBufferSizeInBytes(), BrickBufferInRAM());
    bricksOfBuffersInUsage_ = std::vector<size_t>(ramLimitInBytes_/getBrickBufferSizeInBytes(),0);
}

std::string OctreeBrickPoolManagerDiskSingleThreadedLimitedRam::getDescription() const {
    std::string desc;
    desc += "Single Buffer Size: " + formatMemorySize(getBrickBufferSizeInBytes()) + ", ";
    desc += "Max RAM Usage: " + formatMemorySize(ramLimitInBytes_) + ", ";
    desc += "Memory Allocated: " + formatMemorySize(getBrickPoolMemoryAllocated()) + ", ";
    desc += "Brick Buffers in RAM: " + itos(getNumBrickBuffersInRAM()) + ", ";
    desc += "RAM Used: " + formatMemorySize(getNumBrickBuffersInRAM()*getBrickBufferSizeInBytes());
    return desc;
}

} // namespace
