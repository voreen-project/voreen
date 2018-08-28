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

#ifndef VRN_OCTREEBRICKPOOLMANAGERDISK_H
#define VRN_OCTREEBRICKPOOLMANAGERDISK_H

#include "voreen/core/datastructures/octree/octreebrickpoolmanager.h"
#include "voreen/core/datastructures/octree/brickpoolmanagerqueue.h"

#include <map>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

#include "tgt/assert.h"

namespace voreen {

/**
 * Interface class for disk manager exceptions
 */
class VRN_CORE_API OctreeBrickPoolManagerDiskException : public VoreenException {
public:
    OctreeBrickPoolManagerDiskException(
        const std::string& what = "Something went wrong in the OctreeBrickPoolManagerDisk class!")
        : VoreenException(what)
    {}
};

class VRN_CORE_API BrickIsBeingWrittenException : public OctreeBrickPoolManagerDiskException {
public:
    BrickIsBeingWrittenException(
        const std::string& what = "Requested brick is currently written!")
        : OctreeBrickPoolManagerDiskException(what)
    {}
};

class VRN_CORE_API BrickIsInUseException : public OctreeBrickPoolManagerDiskException {
public:
    BrickIsInUseException(
        const std::string& what = "Requested brick is in use and cannot been written!")
        : OctreeBrickPoolManagerDiskException(what)
    {}
};

class VRN_CORE_API AllBuffersInUseException : public OctreeBrickPoolManagerDiskException {
public:
    AllBuffersInUseException(
        const std::string& what = "All buffers are in use!")
        : OctreeBrickPoolManagerDiskException(what)
    {}
};

/**
 * Class used to load/save brick buffers from/to the disk.
 */
class VRN_CORE_API OctreeBrickPoolManagerDisk : public OctreeBrickPoolManagerBase {
    friend class IncrementalVolumeOctree;

    enum ChannelSelection {
        CHANNEL_1 = 0,
        CHANNEL_2 = 1,
        CHANNEL_3 = 2,
        CHANNEL_4 = 3,
        ALL_CHANNELS = 4
    };

    /**
     * Pair used to count the number of handles using the brick and if the brick is been written.
     * @first count of handles (initialized with 0)
     * @second determines, if the brick is been written (initilaized with false)
     */
    struct BrickEntry {
        uint16_t entry_;

        BrickEntry() : entry_(0) {}

        void increaseInUse(size_t channel);
        void increaseAllInUse();
        void decreaseInUse(size_t channel);
        void decreaseAllInUse();
        bool isInUse(size_t channel, uint16_t threshold = 0);
        bool isInUse(uint16_t threshold = 0);

        void setBeingWritten(bool written, size_t channel);
        void setBeingWritten(bool written);
        bool isBeingWritten(size_t channel);
        bool isBeingWritten();

        inline void increaseInUse(ChannelSelection channels) { (channels == ALL_CHANNELS ? increaseAllInUse() : increaseInUse((size_t)channels)); }
        inline void decreaseInUse(ChannelSelection channels) { (channels == ALL_CHANNELS ? decreaseAllInUse() : decreaseInUse((size_t)channels)); }
        inline bool isInUse(ChannelSelection channels, uint16_t threshold) { return (channels == ALL_CHANNELS ? isInUse(threshold) : isInUse((size_t)channels, threshold)); }
        inline void setBeingWritten(bool written, ChannelSelection channels) { return (channels == ALL_CHANNELS ? setBeingWritten(written) : setBeingWritten(written,(size_t)channels)); }
        inline bool isBeingWritten(ChannelSelection channels) { return (channels == ALL_CHANNELS ? isBeingWritten() : isBeingWritten((size_t)channels)); }
    };

    /**
     * Struct used by the bufferMap to take track of the buffers in RAM
     */
    struct BufferEntry {
        bool isInRAM_;                               //<
        bool mustBeSavedToDisk_;                     //< flag, if the buffer must be saved to disk
        char* data_;                                 //< pointer to the buffer data
        uint8_t inUse_;                              //< counter of handels using the buffer
        BrickEntry* bricksInUse_;                    //< array to all bricks in the buffer (counting handels)
        BrickPoolManagerQueueNode<size_t>* node_;    //< pointer to the queue for least resently used update

        BufferEntry(size_t numberOfBricks, char* data, BrickPoolManagerQueueNode<size_t>* node)
            : isInRAM_(false)
            , mustBeSavedToDisk_(false)
            , data_(data)
            , inUse_(0)
            , bricksInUse_(new BrickEntry[numberOfBricks])
            , node_(node)
        {}
        ~BufferEntry() {
            tgtAssert(inUse_ == 0, "buffer still in use");
            delete[] bricksInUse_;
            delete data_;
        }
    };

public:
    /** Constructor */
    OctreeBrickPoolManagerDisk(const size_t maxSingleBufferSize, const size_t ramLimit,
                               const std::string& brickPoolPath, const std::string& bufferFilePrefix = "");
    /** Destructor */
    ~OctreeBrickPoolManagerDisk();
    /// @see VoreenSerializableObject
    std::string getClassName() const {return "OctreeBrickPoolManagerDisk"; }
    /// @see VoreenSerializableObject
    OctreeBrickPoolManagerDisk* create() const;

    /// Sets the maximum amount of RAM in bytes the brick pool manager is allowed to use.
    void setRAMLimit(size_t ramLimitInBytes);

    //brick interaction
    bool isBrickInRAM(uint64_t virtualMemoryAddress) const;

    const uint16_t* getBrick(uint64_t virtualMemoryAddress, ChannelSelection channels, bool blocking = true) const;
    uint16_t* getWritableBrick(uint64_t virtualMemoryAddress, ChannelSelection channels, bool blocking = true) const;

    void releaseBrick(uint64_t virtualMemoryAddress, ChannelSelection channels, AccessMode mode = READ) const;

    const uint16_t* getBrick(uint64_t virtualMemoryAddress, bool blocking = true) const {
        return getBrick(virtualMemoryAddress, ALL_CHANNELS, blocking);
    }

    uint16_t* getWritableBrick(uint64_t virtualMemoryAddress, bool blocking = true) const {
        return getWritableBrick(virtualMemoryAddress, ALL_CHANNELS, blocking);
    }

    void releaseBrick(uint64_t virtualMemoryAddress, AccessMode mode = READ) const {
        releaseBrick(virtualMemoryAddress, ALL_CHANNELS, mode);
    }

    virtual void flushPoolToDisk(ProgressReporter* progressReporter = 0);

    uint64_t allocateBrick();
    void deleteBrick(uint64_t virtualMemoryAddress);

    /// @see VoreenSerializableObject
    virtual void serialize(Serializer& s) const;
    /// @see VoreenSerializableObject
    virtual void deserialize(Deserializer& s);

    /// general functions
    virtual uint64_t getBrickPoolMemoryUsed() const;
    virtual uint64_t getBrickPoolMemoryAllocated() const;
    virtual std::string getDescription() const;

protected:
    virtual void initialize(size_t brickMemorySizeInByte);
    virtual void deinitialize();
private:
    /// Disk Interaction
    /**
     * Loads a single buffer from the disk.
     * @note This function is not protected by a mutex.
     */
    BufferEntry* loadBufferFromDisk(size_t bufferID, bool blocking, boost::unique_lock<boost::mutex> &lock) const;

    /**
     * Saves a single buffer to the disk.
     * @note This function is not protected by a mutex.
     */
    void saveBufferToDisk(const size_t bufferID) const;

    //--------------------
    //  members
    //--------------------

    /// general defines
    std::string brickPoolPath_;                 //< directory where the buffer files are stored
    std::string bufferFilePrefix_;              //< filename prefix of the buffer files (may be empty)

    size_t maxBufferSizeBytes_;                 //< maximum buffer size in byte (as passed to the constructor)
    size_t singleBufferSizeBytes_;              //< actual size of a single buffer in bytes (next smaller multiple of brick memory size)
    size_t numBrickSlotsPerBuffer_;             //< number of brick slots per buffer

    size_t ramLimitInBytes_;                    //< maximal number of bytes used in the ram
    size_t maxNumBuffersInRAM_;                 //< maximal ram usage maped to map size
    uint64_t nextVirtualMemoryAddress_;         //< virtual memory address of next allocated brick

    std::vector<std::string> bufferFiles_;      //< disk files storing the brick buffers

    ///brick ram management
    mutable size_t numBuffersInRAM_;                               //<
    mutable std::vector<BufferEntry*> bufferVector_;               //<
    mutable BrickPoolManagerQueue<size_t> brickPoolManagerQueue_;  //< least resently used queue
    std::vector<uint64_t> deletedBricks_;                          //< bricks, which had been deleted

    ///multi threaded
    mutable boost::mutex mutex_;                    //< mutex to handle multi-threaded access
    mutable boost::condition_variable cond_;        //< condidion to handle thread sleeping
};

} // namespace

#endif
