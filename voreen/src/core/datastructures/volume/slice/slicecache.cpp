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

#include "voreen/core/datastructures/volume/slice/slicecache.h"
#include "voreen/core/utils/hashing.h"

namespace voreen {

SliceCache::SliceCache(Processor* owner, size_t cacheSize)
    : owner_(owner)
    , cacheSize_(cacheSize)
    , currentBackgroundThread_(0)
    , tmpCache_()
{
    tgtAssert(owner_, "null pointer passed");
}

SliceCache::~SliceCache() {
    clear();
}

size_t SliceCache::getCacheSize() const {
    return cacheSize_;
}

void SliceCache::setCacheSize(size_t cacheSize) {
    cacheSize_ = cacheSize;
    while (slices_.size() > cacheSize_) {
        delete slices_.back().slice_;
        slices_.pop_back();
    }
}

void SliceCache::clear() {
    boost::lock_guard<boost::mutex> lock(cacheAccessMutex_);

    if (currentBackgroundThread_ && currentBackgroundThread_->isRunning())
        currentBackgroundThread_->interruptAndJoin();
    delete currentBackgroundThread_;
    currentBackgroundThread_ = 0;

    for (std::list<CacheEntry>::iterator it = slices_.begin(); it != slices_.end(); ++it) {
        tgtAssert(it->slice_, "cache entry does not store slice");
        delete it->slice_;
    }
    slices_.clear();

    delete tmpCache_.slice_;
    tmpCache_.slice_ = 0;
}

SliceTexture* SliceCache::getVolumeSlice(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray /*= 0*/,
    size_t levelOfDetail /*= 0*/, clock_t timeLimit /*= 0*/, bool* complete /*= 0*/, bool asynchronous /*= false*/) const
{
    tgtAssert(volume, "null pointer passed");

    std::string hash = getHash(volume, alignment, sliceIndex, shiftArray, levelOfDetail);

    // check if slice is present in cache, otherwise create it and add it to cache (if complete)
    SliceTexture* slice = findSliceInCache(hash, true);
    bool sliceComplete = true;

    if (!slice) {

        if (!asynchronous) { // create slice synchronously
            slice = SliceHelper::getVolumeSlice(volume, alignment, sliceIndex, shiftArray, levelOfDetail, timeLimit, &sliceComplete);
            if (slice) {
                addSliceToCache(slice, hash, sliceComplete);
                cleanupCache();
                tgtAssert(slices_.size() <= cacheSize_, "invalid cache size");
            }
        }
        else { // create slice in background thread
            bool sliceCurrentlyCreated = currentBackgroundThread_ && currentBackgroundThread_->hash_ == hash;
            if (!sliceCurrentlyCreated) {

                // interrupt/delete current thread
                if (currentBackgroundThread_) {
                    if (currentBackgroundThread_->isRunning()) {
                        currentBackgroundThread_->interruptAndJoin();
                    }
                    delete currentBackgroundThread_;
                    currentBackgroundThread_ = 0;
                }

                cleanupCache();

                // start new thread
                currentBackgroundThread_ = new SliceCreationBackgroundThread(owner_, this, volume, alignment, sliceIndex, shiftArray, levelOfDetail, hash);
                currentBackgroundThread_->run();
            }

            sliceComplete = false;
            slice = 0;
        }
    }

    if (complete)
        *complete = sliceComplete;

    return slice;
}

SliceTexture* SliceCache::getVolumeSlice(const VolumeBase* volume, tgt::plane pl, float samplingRate, bool asynchronous) const {
    tgtAssert(volume, "null pointer passed");

    std::string hash = getHash(volume, pl, samplingRate);

    // check if slice is present in cache, otherwise create it and add it to cache
    SliceTexture* slice = findSliceInCache(hash, true);

    if (!slice) {
        slice = SliceHelper::getVolumeSlice(volume, pl, samplingRate);
        if (slice) {
            addSliceToCache(slice, hash);
            cleanupCache();
            tgtAssert(slices_.size() <= cacheSize_, "invalid cache size");
        }
    }

    return slice;
}

// private

void SliceCache::addSliceToCache(SliceTexture* slice, const std::string& hash, bool finishedSlice) const {
    boost::lock_guard<boost::mutex> lock(cacheAccessMutex_);

    tgtAssert(slice, "null pointer passed");
    tgtAssert(!hash.empty(), "hash is empty");

    // store unfinished or "not" cached slices to be able to delete them.
    if(!finishedSlice || (cacheSize_ == 0)) {
        delete tmpCache_.slice_;
        tmpCache_.slice_ = slice;
        tmpCache_.hash_ = hash;
    } else { //add slice to cache list
        CacheEntry entry;
        entry.hash_ = hash;
        entry.slice_ = slice;
        slices_.push_front(entry);
    }
}


void SliceCache::cleanupCache() const {
    boost::lock_guard<boost::mutex> lock(cacheAccessMutex_);

    while (slices_.size() > cacheSize_) {
        delete slices_.back().slice_;
        slices_.pop_back();
    }
    tgtAssert(slices_.size() <= cacheSize_, "invalid cache size");
}

SliceTexture* SliceCache::findSliceInCache(const std::string& hash, bool updateUsage) const {
    boost::lock_guard<boost::mutex> lock(cacheAccessMutex_);

    tgtAssert(!hash.empty(), "hash is empty");
    for (std::list<CacheEntry>::iterator it = slices_.begin(); it != slices_.end(); ++it) {
        if (it->hash_ == hash) {
            CacheEntry entry = *it;
            if (updateUsage) { // move found entry to front
                slices_.erase(it);
                slices_.push_front(entry);
            }
            tgtAssert(entry.slice_, "cache entry does not contain slice");
            return entry.slice_;
        }
    }
    return 0;
}

bool SliceCache::hasSliceInCache(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray, size_t levelOfDetail /*= 0*/) const {
    return (findSliceInCache(getHash(volume, alignment, sliceIndex, shiftArray, levelOfDetail), false) != 0);
}

bool SliceCache::hasSliceInCache(const VolumeBase* volume, tgt::plane pl, float samplingRate) const {
    return (findSliceInCache(getHash(volume, pl, samplingRate), false) != 0);
}

std::string SliceCache::getHash(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray, size_t levelOfDetail) const {
    std::ostringstream configStr;
    configStr << volume->getHash();
    configStr << alignment;
    configStr << sliceIndex;
    if(shiftArray) {
        for(size_t c = 0; c < volume->getNumChannels(); c++)
            configStr << shiftArray[c];
    }
    configStr << levelOfDetail;
    return VoreenHash::getHash(configStr.str());
}

std::string SliceCache::getHash(const VolumeBase* volume, tgt::plane pl, float samplingRate) const {
    std::ostringstream configStr;
    configStr << volume->getHash();
    configStr << pl;
    configStr << samplingRate;
    return VoreenHash::getHash(configStr.str());

}

//--------------------------------------------------------------------------
//  Background Thread
//--------------------------------------------------------------------------
SliceCreationBackgroundThread::SliceCreationBackgroundThread(Processor* processor, const SliceCache* sliceCache,
           const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray,
            size_t levelOfDetail, const std::string& hash)
    : BackgroundThread()
    , sliceCache_(sliceCache)
    , volume_(volume)
    , alignment_(alignment)
    , sliceIndex_(sliceIndex)
    , shiftArray_(0)
    , levelOfDetail_(levelOfDetail)
    , hash_(hash)
    {
    if(shiftArray) {
        shiftArray_ = new int[volume->getNumChannels()];
        memcpy(&shiftArray_,&shiftArray,sizeof(int)*volume->getNumChannels());
    }
    //tgtAssert(processor_, "null pointer passed as processor");
    tgtAssert(volume_, "null pointer passed as volume");
    tgtAssert(!hash_.empty(), "passed hash is empty");
}

void SliceCreationBackgroundThread::threadMain() {
    SliceTexture* slice = SliceHelper::getVolumeSlice(volume_, alignment_, sliceIndex_, shiftArray_, levelOfDetail_);
    delete[] shiftArray_;
    if (slice) {
        sliceCache_->addSliceToCache(slice, hash_);
    }
}

} // namespace voreen
