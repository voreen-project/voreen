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

#include "voreen/core/datastructures/octree/octreebrickpoolmanager.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/serialization/serialization.h"

#include <time.h>
#include <limits>
#include <fstream>
#include <stdint.h>
#include <boost/thread/locks.hpp>

#include "tgt/assert.h"
#include "tgt/logmanager.h"
#include "tgt/tgt_math.h"
#include "tgt/filesystem.h"
#include "tgt/vector.h"

using tgt::svec3;
using tgt::vec3;
//-----------------------------------------------------------------------------

namespace voreen {
namespace {
static float brickToNorm(uint16_t val) {
    const float NORM_TO_BRICK_FACTOR = 1.0f/0xffff;
    return static_cast<float>(val) * NORM_TO_BRICK_FACTOR;
}
}

/// BrickPoolBrick -------------------------------------------------------------
BrickPoolBrick::BrickPoolBrick(BrickPoolBrick&& other)
    : addr_(other.addr_)
    , data_(other.data_.voxel(), other.data_.getDimensions(), false) // data is not owned!
    , pool_(other.pool_)
{
    other.addr_ = OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS;
}
BrickPoolBrick& BrickPoolBrick::operator=(BrickPoolBrick&& other) {
    if(&other != this) {
        this->~BrickPoolBrick();
        new(this) BrickPoolBrick(std::move(other));
    }
    return *this;
}

BrickPoolBrick::BrickPoolBrick(uint64_t addr, const tgt::svec3& brickDataSize, const OctreeBrickPoolManagerBase& pool)
    : addr_(addr)
    , data_(pool.getWritableBrick(addr_), brickDataSize, false) // data is not owned!
    , pool_(pool)
{
}
BrickPoolBrick::~BrickPoolBrick() {
    if(addr_ != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
        pool_.releaseBrick(addr_, OctreeBrickPoolManagerBase::WRITE);
    }
}
float BrickPoolBrick::getVoxelNormalized(const tgt::svec3& pos) const {
    return brickToNorm(data_.voxel(pos));
}
VolumeAtomic<uint16_t>& BrickPoolBrick::data() {
    return data_;
}
const VolumeAtomic<uint16_t>& BrickPoolBrick::data() const {
    return data_;
}

/// BrickPoolBrickConst --------------------------------------------------------
BrickPoolBrickConst::BrickPoolBrickConst(BrickPoolBrickConst&& other)
    : addr_(other.addr_)
    , data_(const_cast<uint16_t*>(other.data_.voxel()) /*const_cast is fine, as the VolumeAtomic is const*/, other.data_.getDimensions(), false) // data is not owned!
    , pool_(other.pool_)
{
    other.addr_ = OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS;
}
BrickPoolBrickConst& BrickPoolBrickConst::operator=(BrickPoolBrickConst&& other) {
    if(&other != this) {
        this->~BrickPoolBrickConst();
        new(this) BrickPoolBrickConst(std::move(other));
    }
    return *this;
}

BrickPoolBrickConst::BrickPoolBrickConst(uint64_t addr, const tgt::svec3& brickDataSize, const OctreeBrickPoolManagerBase& pool)
    : addr_(addr)
    , data_(const_cast<uint16_t*>(pool.getBrick(addr_)) /*const_cast is fine, as the VolumeAtomic is const*/, brickDataSize, false) // data is not owned!
    , pool_(pool)
{
}
BrickPoolBrickConst::~BrickPoolBrickConst() {
    if(addr_ != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
        pool_.releaseBrick(addr_, OctreeBrickPoolManagerBase::READ);
    }
}
float BrickPoolBrickConst::getVoxelNormalized(const tgt::svec3& pos) const {
    return brickToNorm(data_.voxel(pos));
}
const VolumeAtomic<uint16_t>& BrickPoolBrickConst::data() const {
    return data_;
}

/// BrickPoolManagerBase -------------------------------------------------------
const std::string OctreeBrickPoolManagerBase::loggerCat_("voreen.OctreeBrickPoolManagerBase");

OctreeBrickPoolManagerBase::OctreeBrickPoolManagerBase()
    : brickMemorySizeInByte_(0)
    , initialized_(false)
{}

OctreeBrickPoolManagerBase::~OctreeBrickPoolManagerBase() {
    if (isInitialized())
        LWARNING("~OctreeBrickPoolManagerBase(): not deinitialized");
}

size_t OctreeBrickPoolManagerBase::getBrickMemorySizeInByte() const {
    return brickMemorySizeInByte_;
}

bool OctreeBrickPoolManagerBase::isInitialized() const {
    return initialized_;
}

void OctreeBrickPoolManagerBase::initialize(size_t brickMemorySizeInByte) {
    if (isInitialized())
        deinitialize();

    brickMemorySizeInByte_ = brickMemorySizeInByte;

    initialized_ = true;
}

void OctreeBrickPoolManagerBase::deinitialize() {
    if (!isInitialized()) {
        LWARNING("deinitialize(): not initialized");
    }
    initialized_ = false;
}

std::string OctreeBrickPoolManagerBase::getDescription() const {
    return "<not available>";
}

void OctreeBrickPoolManagerBase::serialize(Serializer& s) const {
    s.serialize("brickMemorySizeInByte", brickMemorySizeInByte_);
    s.serialize("initialized", initialized_);
}

void OctreeBrickPoolManagerBase::deserialize(Deserializer& s) {
    s.deserialize("brickMemorySizeInByte", brickMemorySizeInByte_);
    s.deserialize("initialized", initialized_);
}

//-------------------------------------------------------------------------------------------------
// OctreeBrickPoolManagerRAM

const std::string OctreeBrickPoolManagerRAM::loggerCat_("voreen.OctreeBrickPoolManagerRAM");

OctreeBrickPoolManagerRAM::OctreeBrickPoolManagerRAM(size_t maxSingleBufferSize)
    : OctreeBrickPoolManagerBase()
    , maxSingleBufferSizeBytes_(maxSingleBufferSize)
    , nextVirtualMemoryAddress_(0)
{
    tgtAssert(maxSingleBufferSize > 0, "max buffer size must be greater zero");
}

OctreeBrickPoolManagerRAM::~OctreeBrickPoolManagerRAM() {
    if (isInitialized())
        deinitialize();
}

OctreeBrickPoolManagerRAM* OctreeBrickPoolManagerRAM::create() const {
    return new OctreeBrickPoolManagerRAM();
}

size_t OctreeBrickPoolManagerRAM::getNumBrickBuffers() const {
    return brickBuffers_.size();
}

uint64_t OctreeBrickPoolManagerRAM::allocateBrick() {
    tgtAssert(isInitialized(), "no initialized");

    boost::lock_guard<boost::mutex> lock(mutex_);

    uint64_t bufferID = nextVirtualMemoryAddress_ / getBrickBufferSizeInBytes();
    uint64_t bufferOffset = nextVirtualMemoryAddress_ % getBrickBufferSizeInBytes();
    while (bufferID >= brickBuffers_.size()){
        char* test = new char[getBrickBufferSizeInBytes()];
        brickBuffers_.push_back(test);
    }
    uint64_t returnValue = nextVirtualMemoryAddress_;
    nextVirtualMemoryAddress_ += getBrickMemorySizeInByte();
    return returnValue;
}

void OctreeBrickPoolManagerRAM::deleteBrick(uint64_t virtualMemoryAddress) {
    //TODO
}

const uint16_t* OctreeBrickPoolManagerRAM::getBrick(uint64_t virtualMemoryAddress, bool /*blocking*/) const {
    // ignore blocking parameter, since we only operate in RAM anyway

    boost::lock_guard<boost::mutex> lock(mutex_);

    tgtAssert(isInitialized(), "no initialized");
    tgtAssert(virtualMemoryAddress == 0 || (virtualMemoryAddress == std::numeric_limits<uint64_t>::max()) || (isMultipleOf(virtualMemoryAddress, (uint64_t)getBrickMemorySizeInByte())),
            "brick virtual memory address is not a multiple of the brick memory size");

    size_t bufferID = static_cast<size_t>(virtualMemoryAddress / getBrickBufferSizeInBytes());
    uint64_t bufferOffset = virtualMemoryAddress % getBrickBufferSizeInBytes();
    if (bufferID >= brickBuffers_.size())
        return 0;
    else
        return reinterpret_cast<uint16_t*>(brickBuffers_[bufferID] + bufferOffset);
}

uint16_t* OctreeBrickPoolManagerRAM::getWritableBrick(uint64_t virtualMemoryAddress, bool blocking) const {
    //boost::lock_guard<boost::mutex> lock(mutex_);

    return const_cast<uint16_t*>(getBrick(virtualMemoryAddress));
}

void OctreeBrickPoolManagerRAM::releaseBrick(uint64_t /*virtualMemoryAddress*/, AccessMode /*mode*/) const {
    //TODO
}

bool OctreeBrickPoolManagerRAM::isBrickInRAM(uint64_t /*virtualMemoryAddress*/) const {
    return true;
}

void OctreeBrickPoolManagerRAM::flushPoolToDisk(ProgressReporter* progressReporter /*= 0*/) {
    // TODO
}

size_t OctreeBrickPoolManagerRAM::getBrickBufferSizeInBytes() const {
    return singleBufferSizeBytes_;
}

uint64_t OctreeBrickPoolManagerRAM::getBrickPoolMemoryAllocated() const {
    return static_cast<uint64_t>(getNumBrickBuffers()) * static_cast<uint64_t>(getBrickBufferSizeInBytes());
}

uint64_t OctreeBrickPoolManagerRAM::getBrickPoolMemoryUsed() const {
    return nextVirtualMemoryAddress_;
}

void OctreeBrickPoolManagerRAM::serialize(Serializer& s) const {

    boost::lock_guard<boost::mutex> lock(mutex_);

    OctreeBrickPoolManagerBase::serialize(s);

    s.serialize("maxSingleBufferSizeBytes", maxSingleBufferSizeBytes_);
    s.serialize("singleBufferSizeBytes",    singleBufferSizeBytes_);
    s.serialize("nextVirtualMemoryAddress", nextVirtualMemoryAddress_);
    s.serialize("numBuffers", brickBuffers_.size());

    // determine output path for brick buffers
    const std::string octreeFile = s.getDocumentPath();
    const std::string octreePath = tgt::FileSystem::dirName(octreeFile);
    if (octreeFile.empty() || !tgt::FileSystem::dirExists(octreePath))
        throw SerializationException("Octree path does not exist: " + octreePath);

    // create brick sub-directory in octree path
    std::string brickBufferPath = tgt::FileSystem::cleanupPath(octreePath + "/" + getClassName());
    if (!tgt::FileSystem::dirExists(brickBufferPath)) {
        if (!tgt::FileSystem::createDirectory(brickBufferPath))
            throw SerializationException("Failed to create brick buffer directory: " + brickBufferPath);
    }
    tgtAssert(tgt::FileSystem::dirExists(brickBufferPath), "brick buffer path not created");

    // serialize brick buffers
    for (size_t i=0; i<brickBuffers_.size(); i++) {
        const std::string bufferFile = brickBufferPath + "/rambuffer_" + itos(i, 10, '0') + ".raw";
        std::fstream fileStream(bufferFile.c_str(), std::ios_base::out | std::ios_base::binary);
        if (fileStream.fail())
            throw SerializationException("Failed to open file '" + bufferFile + "' for writing");

        try {
            fileStream.write(brickBuffers_.at(i), singleBufferSizeBytes_);
        }
        catch (std::exception& e) {
            fileStream.close();
            throw SerializationException("Failed to write RAM buffer to file '" + bufferFile + "': " + std::string(e.what()));
        }
        fileStream.close();
    }

}

void OctreeBrickPoolManagerRAM::deserialize(Deserializer& s) {

    boost::lock_guard<boost::mutex> lock(mutex_);

    OctreeBrickPoolManagerBase::deserialize(s);

    s.deserialize("maxSingleBufferSizeBytes", maxSingleBufferSizeBytes_);
    s.deserialize("singleBufferSizeBytes",    singleBufferSizeBytes_);
    s.deserialize("nextVirtualMemoryAddress", nextVirtualMemoryAddress_);
    int numBuffers = 0;
    s.deserialize("numBuffers", numBuffers);
    if (numBuffers <= 0 || numBuffers > 1<<20)
        throw SerializationException("OctreeBrickPoolManagerRAM: invalid buffer count: " + itos(numBuffers));

    // determine brick buffer path
    const std::string octreeFile = s.getDocumentPath();
    const std::string octreePath = tgt::FileSystem::dirName(octreeFile);
    const std::string brickBufferPath = tgt::FileSystem::cleanupPath(octreePath + "/" + getClassName());
    if (brickBufferPath.empty() || !tgt::FileSystem::dirExists(brickBufferPath))
        throw SerializationException("Brick buffer path does not exist: " + brickBufferPath);

    // load brick buffers from file
    for (int i=0; i<numBuffers; i++) {
        // open brick buffer file
        const std::string bufferFile = tgt::FileSystem::cleanupPath(brickBufferPath + "/rambuffer_" + itos(i, 10, '0') + ".raw");
        std::fstream fileStream(bufferFile.c_str(), std::ios_base::in | std::ios_base::binary);
        if (fileStream.fail())
            throw SerializationException("Failed to open brick buffer file '" + bufferFile + "' for reading");

        // allocate brick buffer
        char* buffer = 0;
        try {
            buffer = new char[singleBufferSizeBytes_];
        }
        catch (std::bad_alloc&) {
            throw SerializationException("Bad allocation when creating brick buffer");
        }
        tgtAssert(buffer, "no buffer");
        brickBuffers_.push_back(buffer);

        // read brick buffer from file
        try {
            fileStream.read(buffer, singleBufferSizeBytes_);
        }
        catch (std::exception& e) {
            fileStream.close();
            throw SerializationException("Failed to read brick buffer from file '" + bufferFile + "': " + std::string(e.what()));
        }
        fileStream.close();
    }

}

void OctreeBrickPoolManagerRAM::initialize(size_t brickMemorySizeInByte) {

    boost::lock_guard<boost::mutex> lock(mutex_);

    OctreeBrickPoolManagerBase::initialize(brickMemorySizeInByte);

    if (maxSingleBufferSizeBytes_ < getBrickMemorySizeInByte())
        throw VoreenException("Max brick buffer size is smaller than the memory size of a single brick "
                              "[" + itos(maxSingleBufferSizeBytes_) + " bytes < " + itos(getBrickMemorySizeInByte()) + " bytes]");

    // round max buffer size down to next multiple of brick memory size
    singleBufferSizeBytes_ = maxSingleBufferSizeBytes_;
    if (!isMultipleOf(singleBufferSizeBytes_, getBrickMemorySizeInByte()))
        singleBufferSizeBytes_ = tgt::ifloor((float)maxSingleBufferSizeBytes_ / (float)getBrickMemorySizeInByte()) * getBrickMemorySizeInByte();
}

void OctreeBrickPoolManagerRAM::deinitialize() {

    boost::lock_guard<boost::mutex> lock(mutex_);

    for (size_t i=0; i<brickBuffers_.size(); i++)
        delete[] brickBuffers_.at(i);
    brickBuffers_.clear();
    nextVirtualMemoryAddress_ = 0;

    OctreeBrickPoolManagerBase::deinitialize();
}

std::string OctreeBrickPoolManagerRAM::getDescription() const {
    std::string desc;
    desc += "Single Buffer Size: " + formatMemorySize(singleBufferSizeBytes_) + ", ";
    desc += "Num Buffers Allocated: " + itos(brickBuffers_.size()) + ", ";
    desc += "Memory Allocated: " + formatMemorySize(getBrickPoolMemoryAllocated());
    return desc;
}

} // namespace
