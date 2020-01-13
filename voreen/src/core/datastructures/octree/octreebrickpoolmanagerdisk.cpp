/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/progressreporter.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "tgt/filesystem.h"

#include <vector>

namespace voreen {

void OctreeBrickPoolManagerDisk::BrickEntry::increaseInUse(size_t channel) {
    tgtAssert(channel < 4, "channel not between 0 and 3!");
    tgtAssert(!isInUse(channel,6), "Overflow. More than 7 threads are using the brick!");
    entry_ += (uint16_t)(1) << (4+3*channel);
}

void OctreeBrickPoolManagerDisk::BrickEntry::increaseAllInUse() {
    tgtAssert(!isInUse(0,6), "Overflow. More than 7 threads are using the brick in 1.channel!");
    tgtAssert(!isInUse(1,6), "Overflow. More than 7 threads are using the brick in 2.channel!");
    tgtAssert(!isInUse(2,6), "Overflow. More than 7 threads are using the brick in 3.channel!");
    tgtAssert(!isInUse(3,6), "Overflow. More than 7 threads are using the brick in 4.channel!");
    entry_ += 9360; //0010010010010000
}

void OctreeBrickPoolManagerDisk::BrickEntry::decreaseInUse(size_t channel) {
    tgtAssert(channel < 4, "channel not between 0 and 3!");
    tgtAssert(isInUse(channel,0), "Underflow. No thread is using the brick!");
    entry_ -= (uint16_t)(1) << (4+3*channel);
}

void OctreeBrickPoolManagerDisk::BrickEntry::decreaseAllInUse() {
    tgtAssert(isInUse(0,0), "Underflow. No thread is using the brick in 1.channel!");
    tgtAssert(isInUse(1,0), "Underflow. No thread is using the brick in 2.channel!");
    tgtAssert(isInUse(2,0), "Underflow. No thread is using the brick in 3.channel!");
    tgtAssert(isInUse(3,0), "Underflow. No thread is using the brick in 4.channel!");
    entry_ -= 9360; //0010010010010000
}

bool OctreeBrickPoolManagerDisk::BrickEntry::isInUse(size_t channel, uint16_t threshold) {
    tgtAssert(channel < 4, "channel not between 0 and 3!");
    uint16_t mask = (uint16_t)(112) << (3*channel); //0000000001110000
    return (((entry_ & mask) >> (4+3*channel)) > threshold);
}

bool OctreeBrickPoolManagerDisk::BrickEntry::isInUse(uint16_t threshold) {
    return (isInUse(0,threshold) || isInUse(1,threshold) || isInUse(2,threshold) || isInUse(3,threshold));

}

void OctreeBrickPoolManagerDisk::BrickEntry::setBeingWritten(bool written, size_t channel) {
    tgtAssert(channel < 4, "channel not between 0 and 3!");
    uint16_t mask = (uint16_t)(1) << channel; //0000000000000001
    if(written)
        entry_ |= mask;
    else
        entry_ &= ~mask;
}

void OctreeBrickPoolManagerDisk::BrickEntry::setBeingWritten(bool written) {
    uint16_t mask = 15; //0000000000001111
    if(written)
        entry_ |= mask;
    else
        entry_ &= ~mask;
}

bool OctreeBrickPoolManagerDisk::BrickEntry::isBeingWritten(size_t channel) {
    tgtAssert(channel < 4, "channel not between 0 and 3!");
    uint16_t mask = (uint16_t)(1) << channel; //0000000000000001
    return ((entry_ & mask) > 0);
}

bool OctreeBrickPoolManagerDisk::BrickEntry::isBeingWritten() {
    uint16_t mask = 15; //0000000000001111
    return ((entry_ & mask) > 0);
}

OctreeBrickPoolManagerDisk::OctreeBrickPoolManagerDisk(const size_t maxBufferSize, const size_t maxRamUsed,
        const std::string& brickPoolPath, const std::string& bufferFilePrefix)
    : OctreeBrickPoolManagerBase()
    , maxBufferSizeBytes_(maxBufferSize)
    , brickPoolPath_(brickPoolPath)
    , bufferFilePrefix_(bufferFilePrefix)
    , ramLimitInBytes_(maxRamUsed)
    , maxNumBuffersInRAM_(0)
    , nextVirtualMemoryAddress_(0)
    , numBuffersInRAM_(0)
{
    brickPoolPath_ = tgt::FileSystem::absolutePath(brickPoolPath);
    bufferFilePrefix_ = (!bufferFilePrefix.empty() ? bufferFilePrefix : "brickbuffer_");
}

OctreeBrickPoolManagerDisk::~OctreeBrickPoolManagerDisk() {
}

OctreeBrickPoolManagerDisk* OctreeBrickPoolManagerDisk::create() const {
    return new OctreeBrickPoolManagerDisk(64<<20, 512<<20, "");
}

//-----------------------------------------------------------------------------------------------------------------------
//      DE-/INITIALIZE
//-----------------------------------------------------------------------------------------------------------------------
void OctreeBrickPoolManagerDisk::initialize(size_t brickMemorySizeInByte) {
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

    // check ram limit vs. buffer size: at least 2 buffers have to fit into the ram
    if (2*singleBufferSizeBytes_ > ramLimitInBytes_)
        throw VoreenException("RAM memory limit is smaller than two times the size of a buffer. At least two buffer files have to fit in the RAM. "
            "[" + itos(ramLimitInBytes_) + " bytes < 2*" + itos(singleBufferSizeBytes_) + " bytes]");

    // define max value
    maxNumBuffersInRAM_ = ramLimitInBytes_/singleBufferSizeBytes_;
}

void OctreeBrickPoolManagerDisk::deinitialize() {
    flushPoolToDisk();

    for (size_t i = 0; i < bufferVector_.size(); i++) {
        delete bufferVector_[i];
        bufferVector_[i] = 0;
    }
    bufferVector_.clear();
    brickPoolManagerQueue_.clear();
    deletedBricks_.clear();
    numBuffersInRAM_ = 0;

    nextVirtualMemoryAddress_ = 0;

    OctreeBrickPoolManagerBase::deinitialize();
}

void OctreeBrickPoolManagerDisk::setRAMLimit(size_t ramLimitInBytes) {
    tgtAssert(ramLimitInBytes > 0, "RAM limit must be larger than 0");

    // check ram limit vs. buffer size: at least 2 buffers have to fit into the ram
    if (2*singleBufferSizeBytes_ > ramLimitInBytes)
        throw VoreenException("RAM memory limit is smaller than two times the size of a buffer. At least two buffer files have to fit in the RAM. "
            "[" + itos(ramLimitInBytes) + " bytes < 2*" + itos(singleBufferSizeBytes_) + " bytes]");

    if (ramLimitInBytes != ramLimitInBytes_) {
        // save non-persistent data to disk
        flushPoolToDisk();

        // delete RAM buffer vectors
        for (size_t i = 0; i < bufferVector_.size(); i++) {
            delete bufferVector_[i];
            bufferVector_[i] = 0;
        }
        bufferVector_.clear();
        brickPoolManagerQueue_.clear();
        deletedBricks_.clear();
        numBuffersInRAM_ = 0;

        // re-init RAM buffer vectors
        ramLimitInBytes_ = ramLimitInBytes;
        maxNumBuffersInRAM_ = ramLimitInBytes_/singleBufferSizeBytes_;

        for (size_t i = 0; i < bufferFiles_.size(); i++)
            bufferVector_.push_back(new BufferEntry(numBrickSlotsPerBuffer_, 0, 0));
    }
}

//-----------------------------------------------------------------------------------------------------------------------
//      DE-/SERIALIZATION
//-----------------------------------------------------------------------------------------------------------------------
void OctreeBrickPoolManagerDisk::serialize(Serializer& s) const {
    OctreeBrickPoolManagerBase::serialize(s);

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
}

void  OctreeBrickPoolManagerDisk::deserialize(Deserializer& s) {
    OctreeBrickPoolManagerBase::deserialize(s);

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
}

//-----------------------------------------------------------------------------------------------------------------------
//      BRICK INTERACTION
//-----------------------------------------------------------------------------------------------------------------------
bool OctreeBrickPoolManagerDisk::isBrickInRAM(uint64_t virtualMemoryAddress) const {
    boost::unique_lock<boost::mutex> lock(mutex_);
    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / singleBufferSizeBytes_);
    tgtAssert(bufferID < bufferVector_.size(), "bufferID not in vector");
    return bufferVector_[bufferID]->isInRAM_;
}

const uint16_t* OctreeBrickPoolManagerDisk::getBrick(uint64_t virtualMemoryAddress, ChannelSelection channels, bool blocking) const {
    boost::unique_lock<boost::mutex> lock(mutex_);
    //constraint SIZE_MAX equals no brick
    if(virtualMemoryAddress == std::numeric_limits<uint64_t>::max()) return 0;

    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / singleBufferSizeBytes_);
    size_t bufferOffset = static_cast<size_t>(virtualMemoryAddress % singleBufferSizeBytes_);
    size_t index = bufferOffset/getBrickMemorySizeInByte();

    tgtAssert(bufferID < bufferVector_.size(), "bufferID not in vector");
    if(bufferVector_[bufferID]->isInRAM_) {
        //increase use counter to prevent buffer removel
        tgtAssert(bufferVector_[bufferID]->inUse_ != 255,"Overflow in use buffer!");
        bufferVector_[bufferID]->inUse_++;
        bufferVector_[bufferID]->bricksInUse_[index].increaseInUse(channels);
        //wait for other threads to release brick
        if(blocking) {
            while(bufferVector_[bufferID]->bricksInUse_[index].isBeingWritten(channels))
                cond_.wait(lock);
        } else {
            if(bufferVector_[bufferID]->bricksInUse_[index].isBeingWritten(channels))
                throw BrickIsBeingWrittenException();
        }
        brickPoolManagerQueue_.pushToFront(bufferVector_[bufferID]->node_);
        return reinterpret_cast<uint16_t*>(bufferVector_[bufferID]->data_ + bufferOffset);
    } else {
        BufferEntry* entry = loadBufferFromDisk(bufferID,blocking,lock); //inUse_(1)
        entry->bricksInUse_[index].increaseInUse(channels);
        return reinterpret_cast<uint16_t*>(entry->data_ + bufferOffset);
    }
}

uint16_t* OctreeBrickPoolManagerDisk::getWritableBrick(uint64_t virtualMemoryAddress, ChannelSelection channels, bool blocking) const {
    boost::unique_lock<boost::mutex> lock(mutex_);
    //constraint SIZE_MAX equals no brick
    if(virtualMemoryAddress == (uint64_t)(-1)) return 0;

    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / singleBufferSizeBytes_);
    size_t bufferOffset = static_cast<size_t>(virtualMemoryAddress % singleBufferSizeBytes_);
    size_t index = bufferOffset/getBrickMemorySizeInByte();

    tgtAssert(bufferID < bufferVector_.size(), "bufferID not in vector");
    if(bufferVector_[bufferID]->isInRAM_) {
        //increase use counter to prevent buffer removel
        tgtAssert(bufferVector_[bufferID]->inUse_ != 255,"Overflow in use buffer!");
        bufferVector_[bufferID]->inUse_++;
        bufferVector_[bufferID]->bricksInUse_[index].increaseInUse(channels);
        if(blocking) {
            while(bufferVector_[bufferID]->bricksInUse_[index].isBeingWritten(channels) ||
                bufferVector_[bufferID]->bricksInUse_[index].isInUse(channels, 1)) //1, since it has been increased
                cond_.wait(lock);
        } else {
            if(bufferVector_[bufferID]->bricksInUse_[index].isBeingWritten(channels)) {
                throw BrickIsBeingWrittenException();
            }
            if(bufferVector_[bufferID]->bricksInUse_[index].isInUse(channels, 1)) {
                throw BrickIsInUseException();
            }
        }
        bufferVector_[bufferID]->bricksInUse_[index].setBeingWritten(true, channels);
        brickPoolManagerQueue_.pushToFront(bufferVector_[bufferID]->node_);
        bufferVector_[bufferID]->mustBeSavedToDisk_ = true;
        return reinterpret_cast<uint16_t*>(bufferVector_[bufferID]->data_ + bufferOffset);
    } else {
        BufferEntry* entry = loadBufferFromDisk(bufferID,blocking,lock); //inUse_(1)
        entry->bricksInUse_[index].increaseInUse(channels);
        entry->bricksInUse_[index].setBeingWritten(true, channels);
        entry->mustBeSavedToDisk_ = true;
        return reinterpret_cast<uint16_t*>(entry->data_ + bufferOffset);
    }
}

void OctreeBrickPoolManagerDisk::releaseBrick(uint64_t virtualMemoryAddress, ChannelSelection channels, AccessMode mode) const {
    boost::unique_lock<boost::mutex> lock(mutex_);

    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / singleBufferSizeBytes_);
    size_t bufferOffset = static_cast<size_t>(virtualMemoryAddress % singleBufferSizeBytes_);
    size_t index = bufferOffset/getBrickMemorySizeInByte();

    tgtAssert(bufferID < bufferVector_.size(), "bufferID not in vector");
    if(bufferVector_[bufferID]->isInRAM_) {
        if(mode == WRITE)
            bufferVector_[bufferID]->bricksInUse_[index].setBeingWritten(false,channels);
        tgtAssert(bufferVector_[bufferID]->bricksInUse_[index].isInUse(channels, 0), "releaseBrick called on brick without being in use!");
        bufferVector_[bufferID]->bricksInUse_[index].decreaseInUse(channels);
        tgtAssert(bufferVector_[bufferID]->inUse_ > 0, "releaseBrick called on brick without being in use!");
        bufferVector_[bufferID]->inUse_--;
        tgtAssert(!(!(bufferVector_[bufferID]->bricksInUse_[index].isInUse(channels,0)) && (bufferVector_[bufferID]->bricksInUse_[index].isBeingWritten(channels))), "brick has no handle but is still being written");
        cond_.notify_all();
    } else {
        tgtAssert(false, "buffer not in RAM!")
    }
}

uint64_t OctreeBrickPoolManagerDisk::allocateBrick() {
    boost::unique_lock<boost::mutex> lock(mutex_);
    //case1: actual buffer is not full -> use it
    //case2: we have free bricks -> use them
    //case3: we have to allocate a new buffer

    //case 1
    if (nextVirtualMemoryAddress_%singleBufferSizeBytes_ != 0) {
        uint64_t returnValue = nextVirtualMemoryAddress_;
        nextVirtualMemoryAddress_ += static_cast<uint64_t>(getBrickMemorySizeInByte());
        return returnValue;
    } else //case2
    if (!deletedBricks_.empty()) {
        uint64_t returnValue = deletedBricks_.back();
        deletedBricks_.pop_back();
        return returnValue;
    } else { //case3

/*        std::ofstream outfile(path.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if(outfile.fail())
            throw VoreenException("Could not open: " + path.str());
        outfile << itos(0,(int)singleBufferSizeBytes_);
        outfile.close();*/


        std::stringstream path;
        path << brickPoolPath_ << "/" << bufferFilePrefix_ << itos(bufferFiles_.size(), 10);
        bufferFiles_.push_back(path.str());

        bufferVector_.push_back(new BufferEntry(numBrickSlotsPerBuffer_, 0, 0));

        if(numBuffersInRAM_ == maxNumBuffersInRAM_) {
            //find LRU buffer
            size_t removeBuffer = brickPoolManagerQueue_.last_->previous_->data_;
            tgtAssert(bufferVector_.size() > removeBuffer, "buffer is not in ram!")
            if(bufferVector_[removeBuffer]->inUse_ > 0) {
                tgtAssert(false,"All bricks are in use!");
                LERROR("All bricks are in use!");
                throw VoreenException("All bricks are in use!");
            }
            //safe old buffer
            if(bufferVector_[removeBuffer]->mustBeSavedToDisk_)
                saveBufferToDisk(removeBuffer);
            //clean up
            if(removeBuffer != brickPoolManagerQueue_.removeLast()) {
                tgtAssert(false, "something went wrong!");
                LERROR("something went wrong!");
            }
            delete[] bufferVector_[removeBuffer]->data_;
            bufferVector_[removeBuffer]->data_ = 0;
            bufferVector_[removeBuffer]->node_ = 0;
            bufferVector_[removeBuffer]->isInRAM_ = false;
            //LERROR("kicked: " << removeBuffer << " loaded: " << bufferID);
            numBuffersInRAM_--;
        }

        char* buffer = 0;
        try {
            buffer = new char[singleBufferSizeBytes_];
        } catch(std::bad_alloc& e) {
            tgtAssert(false,e.what());
            LERROR(e.what());
            throw VoreenException(e.what());
        }

        size_t bufferID = bufferVector_.size()-1;
        BrickPoolManagerQueueNode<size_t>* node = brickPoolManagerQueue_.insertToFront(bufferID);
        bufferVector_.back()->data_ = buffer;
        bufferVector_[bufferID]->isInRAM_ = true;
        bufferVector_[bufferID]->inUse_ = 0;
        bufferVector_[bufferID]->mustBeSavedToDisk_ = true;
        bufferVector_[bufferID]->node_ = node;
        numBuffersInRAM_++;


        uint64_t returnValue = nextVirtualMemoryAddress_;
        nextVirtualMemoryAddress_ += static_cast<uint64_t>(getBrickMemorySizeInByte());
        return returnValue;
    }
}

void OctreeBrickPoolManagerDisk::deleteBrick(uint64_t virtualMemoryAddress) {
    boost::unique_lock<boost::mutex> lock(mutex_);

    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / singleBufferSizeBytes_);
    size_t bufferOffset = static_cast<size_t>(virtualMemoryAddress % singleBufferSizeBytes_);
    size_t index = bufferOffset/getBrickMemorySizeInByte();

    //wait for brick to be not used
    tgtAssert(bufferID < bufferVector_.size(), "bufferID not in vector");
    if(bufferVector_[bufferID]->isInRAM_) {
        while(bufferVector_[bufferID]->bricksInUse_[index].isInUse((uint16_t) 0)) {
                cond_.wait(lock);
        }
    }

    deletedBricks_.push_back(virtualMemoryAddress);
}

//-----------------------------------------------------------------------------------------------------------------------
//      DISK INTERACTION
//-----------------------------------------------------------------------------------------------------------------------
OctreeBrickPoolManagerDisk::BufferEntry* OctreeBrickPoolManagerDisk::loadBufferFromDisk(size_t bufferID, bool blocking, boost::unique_lock<boost::mutex> &lock) const {
    tgtAssert(bufferID < bufferVector_.size(), "bufferID not in vector");
    //check, if bufferID is valid
    if (bufferID >= bufferFiles_.size()) {
        tgtAssert(false,"Buffer has not been created");
        LERROR("loadBrickFromDisk(): Buffer has not been created");
        throw VoreenException("Buffer has not been created");
    }
    else {
        //kick out old buffer, if RAM is to full
        if(numBuffersInRAM_ == maxNumBuffersInRAM_) {
            //find LRU buffer
            size_t removeBuffer = brickPoolManagerQueue_.last_->previous_->data_;
            tgtAssert(bufferVector_.size() > removeBuffer, "buffer is not in ram!")
            if(blocking) {
                while (bufferVector_[removeBuffer]->inUse_ > 0) {
                    cond_.wait(lock);
                    removeBuffer = brickPoolManagerQueue_.last_->previous_->data_;
                }
            } else {
                if(bufferVector_[removeBuffer]->inUse_ > 0) {
                    throw AllBuffersInUseException();
                }
            }
            //safe old buffer
            if(bufferVector_[removeBuffer]->mustBeSavedToDisk_)
                saveBufferToDisk(removeBuffer);
            //clean up
            if(removeBuffer != brickPoolManagerQueue_.removeLast()) {
                tgtAssert(false, "something went wrong!");
                LERROR("something went wrong!");
            }
            delete[] bufferVector_[removeBuffer]->data_;
            bufferVector_[removeBuffer]->data_ = 0;
            bufferVector_[removeBuffer]->node_ = 0;
            bufferVector_[removeBuffer]->isInRAM_ = false;
            //LERROR("kicked: " << removeBuffer << " loaded: " << bufferID);
            numBuffersInRAM_--;
        }

        const std::string bufferFile = bufferFiles_.at(bufferID);
        tgtAssert(!bufferFile.empty(), "buffer file path is empty");
        if (!tgt::FileSystem::fileExists(bufferFile)) {
            tgtAssert(false,"Buffer file does not exists!");
            LERROR("Buffer file does not exists!");
            throw VoreenException("Buffer file does not exists!");
        }
        char* buffer = 0;
        try {
            buffer = new char[singleBufferSizeBytes_];
        } catch(std::bad_alloc& e) {
            tgtAssert(false,e.what());
            LERROR(e.what());
            throw VoreenException(e.what());
        }
        std::ifstream infile(bufferFile.c_str(), std::ios::in | std::ios::binary);
        if(infile.fail()) {
            delete buffer;
            tgtAssert(false,"Could not open buffer file!");
            LERROR("Could not open buffer file!");
            throw VoreenException("Could not open buffer file!");
        }
        infile.read(buffer,singleBufferSizeBytes_);
        tgtAssert(!infile.bad(), "reading from disk went wrong");
        infile.close();

        BrickPoolManagerQueueNode<size_t>* node = brickPoolManagerQueue_.insertToFront(bufferID);
        bufferVector_[bufferID]->data_ = buffer;
        bufferVector_[bufferID]->isInRAM_ = true;
        bufferVector_[bufferID]->inUse_ = 1;
        bufferVector_[bufferID]->mustBeSavedToDisk_ = false;
        bufferVector_[bufferID]->node_ = node;
        numBuffersInRAM_++;

        return bufferVector_[bufferID];
    }
    tgtAssert(false, "Should not get here");
    throw VoreenException("Should not get here");
}

void OctreeBrickPoolManagerDisk::saveBufferToDisk(const size_t bufferID) const {
    tgtAssert(bufferID < bufferVector_.size(), "bufferID not in vector");
    if (bufferID >= bufferFiles_.size()) {
        tgtAssert(false,"Buffer has not been created");
        LERROR("saveBrickToDisk(): Buffer has not been created");
        return ;
    }
    else {
        const std::string bufferFile = bufferFiles_.at(bufferID);
        tgtAssert(!bufferFile.empty(), "buffer file path is empty");
        //may not exist, since allocate does not write the file
        /*if (!tgt::FileSystem::fileExists(bufferFile)) {
            tgtAssert(false,"Buffer file does not exists!");
            LERROR("Buffer file does not exists!");
            return ;
        }*/

        std::ofstream outfile(bufferFile.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if(outfile.fail()) {
            tgtAssert(false,"Could not open buffer file!");
            LERROR("Could not open buffer file!");
            return ;
        }
        outfile.write(bufferVector_[bufferID]->data_,singleBufferSizeBytes_);
        tgtAssert(!outfile.bad(), "writing brick to disk went wrong");
        outfile.close();
        bufferVector_[bufferID]->mustBeSavedToDisk_ = false;
    }
}

void OctreeBrickPoolManagerDisk::flushPoolToDisk(ProgressReporter* progressReporter /*= 0*/) {
    // determine number of buffers to be saved to disk
    size_t numToBeSavedToDisk = 0;
    for (size_t i = 0; i < bufferVector_.size(); i++) {
        if (bufferVector_[i]->mustBeSavedToDisk_)
            numToBeSavedToDisk++;
    }

    // save buffers to disk
    size_t numSavedToDisk = 0;
    for (size_t i = 0; i < bufferVector_.size(); i++) {
        if (bufferVector_[i]->mustBeSavedToDisk_) {
            tgtAssert(bufferVector_[i]->isInRAM_,"buffer not in ram!");
            saveBufferToDisk(i);
            numSavedToDisk++;
            if (progressReporter)
                progressReporter->setProgress((float)(numSavedToDisk) / (float)numToBeSavedToDisk);
            //LINFO("Saved " << numSavedToDisk << "/" << numToBeSavedToDisk);
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------
//      GENERAL FUNCTIONS
//-----------------------------------------------------------------------------------------------------------------------
uint64_t OctreeBrickPoolManagerDisk::getBrickPoolMemoryUsed() const {
    boost::unique_lock<boost::mutex> lock(mutex_);
    return static_cast<size_t>(numBuffersInRAM_ * singleBufferSizeBytes_);
}

std::string OctreeBrickPoolManagerDisk::getDescription() const {
    // TODO
    return "not available";
}

uint64_t OctreeBrickPoolManagerDisk::getBrickPoolMemoryAllocated() const {
    return (uint64_t)bufferFiles_.size() * (uint64_t)singleBufferSizeBytes_;
}

} // namespace
