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

#ifndef VRN_OCTREEBRICKPOOLMANAGERMMAP_H
#define VRN_OCTREEBRICKPOOLMANAGERMMAP_H

#include "voreen/core/datastructures/octree/octreebrickpoolmanager.h"
#include "voreen/core/datastructures/octree/brickpoolmanagerqueue.h"

#include <map>
#include <atomic>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/thread.hpp>

#include "tgt/assert.h"

namespace voreen {

struct OctreeBrickPoolManagerMmapStorageIndex {
    OctreeBrickPoolManagerMmapStorageIndex(uint64_t brickAddr);
    uint64_t fileIndex_;
    uint64_t brickIndexInFile_;
};

/**
 * Class used to load/save brick buffers from/to the disk.
 */
class VRN_CORE_API OctreeBrickPoolManagerMmap : public OctreeBrickPoolManagerBase {

public:
    /** Constructor */
    OctreeBrickPoolManagerMmap(const std::string& brickPoolPath, const std::string& bufferFilePrefix);
    /** Destructor */
    ~OctreeBrickPoolManagerMmap();
    /// @see VoreenSerializableObject
    std::string getClassName() const {return "OctreeBrickPoolManagerMmap"; }
    /// @see VoreenSerializableObject
    OctreeBrickPoolManagerMmap* create() const;

    //brick interaction
    bool isBrickInRAM(uint64_t virtualMemoryAddress) const;
    const uint16_t* getBrick(uint64_t virtualMemoryAddress, bool blocking = true) const;

    uint16_t* getWritableBrick(uint64_t virtualMemoryAddress, bool blocking = true) const;

    virtual void releaseBrick(uint64_t virtualMemoryAddress, AccessMode mode = READ) const;

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

    std::string storageFileName(uint64_t fileIndex) const;

protected:
    virtual void initialize(size_t brickMemorySizeInByte);
    virtual void deinitialize();
private:

    std::string brickPoolPath_;                 //< directory where the buffer files are stored
    std::string bufferFilePrefix_;              //< filename prefix of the buffer files (may be empty)

    std::atomic<uint64_t> freeListHead_;      // Next reusable previously deleted block. This block contains the address to the next block in the list.
    std::atomic<uint64_t> nextNewBrickAddr_; // Addr (index) of the next (entirely new) brick that can be allocated

    mutable boost::shared_mutex storageMutex_; // R/RW for storage_, i.e,. either: read (shared) or modify (exclusive)
    mutable std::vector<boost::iostreams::mapped_file> storage_;
};

} // namespace

#endif
