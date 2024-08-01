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

#ifndef VRN_OCTREEBRICKPOOLMANAGERDISKSINGLETHREADED_H
#define VRN_OCTREEBRICKPOOLMANAGERDISKSINGLETHREADED_H

#include "voreen/core/datastructures/octree/octreebrickpoolmanager.h"

#include <iostream>
#include <fstream>
#include <set>
namespace voreen {

/**
 * Class used to load/save brick buffers from/to the disk.
 */
class VRN_CORE_API OctreeBrickPoolManagerDiskSingleThreaded : public OctreeBrickPoolManagerBase {

public:
    /**
     * @param maxSingleBufferSize size of a single brick buffer in bytes.
     *  Note that a brick buffer can potentially store multiple bricks.
     * @param brickPoolPath location of the brick buffer files on the disk. When \p brickPoolExists is true, all files at the brickPoolPath
     *  matching the \p bufferFilePrefix are considered buffer files, ordered alphabetically.
     */
    OctreeBrickPoolManagerDiskSingleThreaded(const size_t maxSingleBufferSize,
        const std::string& brickPoolPath, const std::string& bufferFilePrefix = "");
    virtual ~OctreeBrickPoolManagerDiskSingleThreaded();

    /// Returns the byte size of a single brick buffer (a brick buffer can store multiple bricks).
    size_t getBrickBufferSizeInBytes() const;

    /// Returns the number of allocated brick buffers (a brick buffer can store multiple bricks).
    virtual size_t getNumBrickBuffers() const;

    /// Returns the number of bricks per buffers (a brick buffer can store multiple bricks).
    virtual size_t getNumBrickSlotsPerBuffer() const;

    /// Returns the paths to the disk files holding the brick buffers.
    std::vector<std::string> getBufferFiles() const;

    ///sets the brickPath and the prefix
    virtual void initialize(size_t brickMemorySizeInByte);
    virtual void deinitialize();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:
    ///allocates the next disk buffer
    virtual char* allocateNextDiskBuffer();

    ///saves a buffer to the disk
    virtual void saveBufferToDisk(size_t bufferIndex, const char* buffer) const;

    ///loads a buffer from the disk
    virtual char* loadDiskBuffer(size_t bufferIndex) const;

    std::vector<std::string> bufferFiles_;      //< disk files storing the brick buffers
    std::string brickPoolPath_;                 //< directory where the buffer files are stored
    std::string bufferFilePrefix_;              //< filename prefix of the buffer files (may be empty)

    static const std::string loggerCat_;

private:
    size_t maxBufferSizeBytes_;                 //< maximum buffer size in byte (as passed to the constructor)
    size_t singleBufferSizeBytes_;              //< actual size of a single buffer in bytes (next smaller multiple of brick memory size)
    size_t numBrickSlotsPerBuffer_;             //< number of brick slots per buffer
};

//-------------------------------------------------------------------------------------------------

class VRN_CORE_API OctreeBrickPoolManagerDiskSingleThreadedLimitedRam : public OctreeBrickPoolManagerDiskSingleThreaded {
public:
    /**
     * Different strategies for handling freeBricks
     */
    enum AllocateStrategy {
        IGNORE_DELETED_BRICKS,
        USE_DELETED_BRICKS
    };

    OctreeBrickPoolManagerDiskSingleThreadedLimitedRam(const size_t maxSingleBufferSize, const size_t ramLimitInBytes,
                                            const std::string& brickPoolPath, const std::string& bufferFilePrefix = "",
                                            const AllocateStrategy strategy = IGNORE_DELETED_BRICKS);
    virtual ~OctreeBrickPoolManagerDiskSingleThreadedLimitedRam();
    OctreeBrickPoolManagerDiskSingleThreadedLimitedRam* create() const;

    std::string getClassName() const { return "OctreeBrickPoolManagerDiskSingleThreaded"; }

    /// Sets the maximum amount of RAM in bytes the brick pool manager is allowed to use.
    void setRAMLimit(size_t ramLimitInBytes);

    /// Returns the maximum amount of RAM in bytes the brick pool manager is allowed to use.
    size_t getRAMLimit() const;

    //----------------------------
    //  Interface functions
    //----------------------------
    virtual uint64_t allocateBrick();
    virtual void deleteBrick(uint64_t virtualMemoryAddress);
    virtual const uint16_t* getBrick(uint64_t virtualMemoryAddress, bool blocking = true) const;
    virtual uint16_t* getWritableBrick(uint64_t virtualMemoryAddress, bool blocking = true) const;
    virtual void releaseBrick(uint64_t virtualMemoryAddress, AccessMode mode = READ) const;
    virtual bool isBrickInRAM(uint64_t virtualMemoryAddress) const;
    virtual void flushPoolToDisk(ProgressReporter* progressReporter = 0);
    virtual uint64_t getBrickPoolMemoryUsed() const;
    virtual uint64_t getBrickPoolMemoryAllocated() const;

    size_t getNumBrickBuffersInRAM() const;

    virtual std::string getDescription() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:
    /// Calculates the maximum buffers keeped in ram
    virtual void initialize(size_t brickMemorySizeInByte);
    virtual void deinitialize();

        /**
         * Struct for the buffers keeped in ram
         */
        struct BrickBufferInRAM{
        char* pointer_;         //< pointer to the ram buffer
        size_t index_;          //< index of the buffer
        bool inUsage_;          //< determines if a brick of the buffer is in use atm
        size_t lastUsed_;       //< last time the buffer has been used
        bool writable_;         //< is buffer writable i.e. it has to be saved after usage

        BrickBufferInRAM() : pointer_(0), index_(-1), inUsage_(false), lastUsed_(0), writable_(false) {}
        BrickBufferInRAM(char* pointer, size_t index, bool inUsage, size_t lastUsed, bool writable)
            : pointer_(pointer), index_(index), inUsage_(inUsage), lastUsed_(lastUsed), writable_(writable) {}
        ~BrickBufferInRAM() {delete[] pointer_;}
        };

    static const std::string loggerCat_;

private:
    ///allocateBrick uses this, if allocateStategy is IGNORE_FREE_BRICKS
    uint64_t ignoreFreeBricks();
    ///allocateBrick uses this, if allocateStategy is USE_FREE_BRICKS
    uint64_t useFreeBricks();
    ///gets the index of the BrickBufferInRAM according to the virtual memory address (used in get(w)Brick)
    size_t loadBrickBufferIntoRAM(size_t bufferIndex) const;
    ///checks, if the buffer is in ram and returs its index or -1
    size_t isBufferInRAM(size_t bufferID) const;
    ///marks a brick as used
    void setBrickUsed(uint64_t virtualMemoryAddress, size_t bufferIndex) const;
    ///this function handels overflows of the timestemp
    size_t increaseTimestemp() const;

    //----------------------------
    //  Member
    //----------------------------
    mutable std::vector<BrickBufferInRAM> brickBuffersInRAM_;  //< pointer to the buffers storing the bricks
    //mutable std::map<size_t,size_t> bufferRAMMap_;           //< map to check, which buffers are in RAM bufferID->bufferIndex

    std::vector<uint64_t> deletedBricks_;               //< bricks, which had been deleted
    mutable std::set<uint64_t> bricksInUsage_;            //< bricks, which are in usage
    mutable std::vector<size_t> bricksOfBuffersInUsage_;        //< vector to check, if a buffer is in usage

    size_t ramLimitInBytes_;                            //< maximal number of bytes used in the ram
    uint64_t nextVirtualMemoryAddress_;                 //< virtual memory address of next allocated brick
    mutable size_t nextTimestemp_;                      //< next timestemp used to mark least used buffer

    AllocateStrategy allocateStrategy_;                 //< strategy for using free bricks

};

} // namespace

#endif
