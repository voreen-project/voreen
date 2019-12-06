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

#ifndef VRN_OCTREEBRICKPOOLMANAGER_H
#define VRN_OCTREEBRICKPOOLMANAGER_H

#include "voreen/core/utils/exception.h"
#include "voreen/core/voreenobject.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <vector>
#include <string>

#include <boost/thread/mutex.hpp>
#include "tgt/vector.h"

namespace voreen {

class ProgressReporter;

class OctreeBrickPoolManagerBase;

// Brick smart pointer that automatically releases its content on destruction. Read/write variant.
//
// For obvious reasons, this struct must not outlive the OctreeBrickPoolManager passen on creation.
struct BrickPoolBrick {
    BrickPoolBrick() = delete;
    BrickPoolBrick(const BrickPoolBrick&) = delete;
    BrickPoolBrick& operator=(const BrickPoolBrick&) = delete;

    BrickPoolBrick(uint64_t addr, const tgt::svec3& brickDataSize, const OctreeBrickPoolManagerBase& pool);
    ~BrickPoolBrick();
    float getVoxelNormalized(const tgt::svec3& pos) const;
    VolumeAtomic<uint16_t>& data();
    const VolumeAtomic<uint16_t>& data() const;

private:
    uint64_t addr_;
    VolumeAtomic<uint16_t> data_;
    const OctreeBrickPoolManagerBase& pool_;
};

// Brick smart pointer that automatically releases its content on destruction. Read-only variant.
//
// For obvious reasons, this struct must not outlive the OctreeBrickPoolManager passen on creation.
struct BrickPoolBrickConst {
    BrickPoolBrickConst() = delete;
    BrickPoolBrickConst(const BrickPoolBrickConst&) = delete;
    BrickPoolBrickConst& operator=(const BrickPoolBrickConst&) = delete;

    BrickPoolBrickConst(uint64_t addr, const tgt::svec3& brickDataSize, const OctreeBrickPoolManagerBase& pool);
    ~BrickPoolBrickConst();

    float getVoxelNormalized(const tgt::svec3& pos) const;
    const VolumeAtomic<uint16_t>& data() const;

private:
    uint64_t addr_;
    const VolumeAtomic<uint16_t> data_;
    const OctreeBrickPoolManagerBase& pool_;
};

/**
 * Base class for brick pool managers that provide a virtual memory space for octree bricks.
 * The brick pool may be either completely stored in the RAM, or the bricks are streamed
 * from the disk.
 *
 * @see OctreeBrickPoolManagerRAM
 * @see OctreeBrickPoolManagerDisk
 */
class VRN_CORE_API OctreeBrickPoolManagerBase : public VoreenSerializableObject {

public:

    static const uint64_t NO_BRICK_ADDRESS = std::numeric_limits<uint64_t>::max();

    enum AccessMode{
        READ,
        WRITE
    };

    OctreeBrickPoolManagerBase();
    virtual ~OctreeBrickPoolManagerBase();

    /**
     * Allocates a new brick and returns a pointer to the virtual memory.
     * To get the real brick pointer, getBrick() or getWritableBrick() has to be called
     *
     * @return pointer to the virtual adress of the allocated brick.
     */
    virtual uint64_t allocateBrick() = 0;

    /**
     * Free the memory of the brick.
     */
    virtual void deleteBrick(uint64_t virtualMemoryAddress) = 0;

    /**
     * Returns a read-only pointer to the brick stored at the passed virtual memory address.
     * If the memory address is unallocated, the null pointer is returned.
     *  @note The brick will be marked as in use.
     */
    virtual const uint16_t* getBrick(uint64_t virtualMemoryAddress, bool blocking = true) const = 0;

    /**
     * Returns a writable pointer to the brick stored at the passed virtual memory address.
     * If the memory address is unallocated, the null pointer is returned.
     * @note The brick will be marked as in use.
     */
    virtual uint16_t* getWritableBrick(uint64_t virtualMemoryAddress, bool blocking = true) const = 0;

    /**
     * Marks the brick as no longer in use.
     */
    virtual void releaseBrick(uint64_t virtualMemoryAddress, AccessMode mode = READ) const = 0;

    /**
     * Returns whether the brick stored at the passed virtual memory address
     * is currently loaded into RAM. Use this, if you want to prevent disk
     * access when acquiring a brick,
     */
    virtual bool isBrickInRAM(uint64_t virtualMemoryAddress) const = 0;

    /**
     * Flushes all unwritten bricks to the disk, thereby making sure that the
     * entire brick pool is persistently stored.
     */
    virtual void flushPoolToDisk(ProgressReporter* progressReporter = 0) = 0;

    /// Returns the memory size of one brick in byte.
    size_t getBrickMemorySizeInByte() const;

    /// Returns the amount of memory in bytes that has been allocated for the entire brick pool.
    virtual uint64_t getBrickPoolMemoryAllocated() const = 0;

    /// Returns the amount of memory in bytes that is actually used for storing the bricks.
    virtual uint64_t getBrickPoolMemoryUsed() const = 0;

    /// Returns a string containing the major properties of the brick pool manager.
    virtual std::string getDescription() const;

    /// Returns true, if the pool manager has been initialized (by the octree).
    bool isInitialized() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

    /**
     * Initializes and checks the brick pool. Is to be called by the Octree when the brick pool manager is assigned to it.
     *
     * @param brickMemorySizeInByte memory size of a brick in byte
     *
     * @throw VoreenException If initialization failed.
     */
    virtual void initialize(size_t brickMemorySizeInByte);

    /// Frees allocated resources. Is to be called by the Octree.
    virtual void deinitialize();

protected:
    static const std::string loggerCat_;

private:
    size_t brickMemorySizeInByte_;  ///< number of bytes of one brick

    bool initialized_;              ///< true, if the manager has been initialized (by the octree)
};

//-------------------------------------------------------------------------------------------------

/**
 * Basic brick pool manager that stores the entire brick in RAM.
 */
class VRN_CORE_API OctreeBrickPoolManagerRAM : public OctreeBrickPoolManagerBase {

public:
    /**
     * @param maxSingleBuff Maximum size of a *single* brick buffer in byte. A brick buffer can potentially
     *  store multiple bricks. The actual buffer size is determined by rounding down to the nearest
     *  multiple of the brick memory size.
     */
    OctreeBrickPoolManagerRAM(size_t maxSingleBufferSize = 1<<26 /* 64 MB */);
    virtual ~OctreeBrickPoolManagerRAM();
    OctreeBrickPoolManagerRAM* create() const;

    std::string getClassName() const { return "OctreeBrickPoolManagerRAM"; }

    virtual uint64_t allocateBrick();
    virtual void deleteBrick(uint64_t virtualMemoryAddress);

    virtual const uint16_t* getBrick(uint64_t virtualMemoryAddress, bool blocking = true) const;
    virtual uint16_t* getWritableBrick(uint64_t virtualMemoryAddress, bool blocking = true) const;

    virtual void releaseBrick(uint64_t virtualMemoryAddress, AccessMode mode = READ) const;

    virtual bool isBrickInRAM(uint64_t virtualMemoryAddress) const;

    virtual void flushPoolToDisk(ProgressReporter* progressReporter = 0);

    /// Returns the byte size of a single brick buffer (a brick buffer can store multiple bricks).
    size_t getBrickBufferSizeInBytes() const;

    /// Returns the number of allocated brick buffers (a brick buffer can store multiple bricks).
    size_t getNumBrickBuffers() const;

    /// Returns the amount of memory in bytes that has been allocated for the entire brick pool.
    virtual uint64_t getBrickPoolMemoryAllocated() const;

    /// Returns the amount of memory in bytes that is actually used for storing the bricks.
    virtual uint64_t getBrickPoolMemoryUsed() const;

    virtual std::string getDescription() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

protected:
    virtual void initialize(size_t brickByteSize);
    virtual void deinitialize();

    static const std::string loggerCat_;

private:
    size_t maxSingleBufferSizeBytes_;     ///< maximum single buffer size in byte (as passed to the constructor)
    size_t singleBufferSizeBytes_;        ///< actual size of a single buffer in bytes (next smaller multiple of brick memory size)
    uint64_t nextVirtualMemoryAddress_;   ///< virtual memory address of next allocated brick

    std::vector<char*> brickBuffers_;     ///< pointer to the buffers storing the bricks

    mutable boost::mutex mutex_;          ///< mutex for handling multi-threaded brick access
};


} // namespace

#endif
