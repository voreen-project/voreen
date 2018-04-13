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

#ifndef VRN_SLICECACHE_H
#define VRN_SLICECACHE_H

#include "voreen/core/datastructures/volume/slice/slicehelper.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"

#include "voreen/core/utils/backgroundthread.h"
#include <boost/thread/mutex.hpp>

namespace voreen {

class Processor;
class SliceCreationBackgroundThread;

/**
 * Fixed-size LRU cache for 2D volume slices. On each getVolumeSlice() call the SliceCache checks whether
 * a slice with the specified parameters is present in the cache and returns it, if this is the case.
 * If not, the slice is extracted from the underlying volume via the SliceHelper and stored in the cache.
 *
 * Uncomplete slices (see VolumeOctree) will be stored in the tmpCache.
 */
class VRN_CORE_API SliceCache {
    friend class SliceCreationBackgroundThread;
public:
    SliceCache(Processor* owner, size_t cacheSize);
    ~SliceCache();

    size_t getCacheSize() const;
    void setCacheSize(size_t cacheSize);
    void clear();

    SliceTexture* getVolumeSlice(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray,
        size_t levelOfDetail = 0, clock_t timeLimit = 0, bool* complete = 0, bool asynchronous = false) const;
    SliceTexture* getVolumeSlice(const VolumeBase* volume, tgt::plane pl, float samplingRate, bool asynchronous = false) const;

    bool hasSliceInCache(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray, size_t levelOfDetail = 0) const;
    bool hasSliceInCache(const VolumeBase* volume, tgt::plane pl, float samplingRate) const;

private:
    struct CacheEntry {
        std::string hash_;  //< consisting of the input volume hash and a hash of the slice parameters
        SliceTexture* slice_;
        CacheEntry() : hash_(""), slice_(0) {}
    };

    void addSliceToCache(SliceTexture* slice, const std::string& hash, bool finishedSlice = true) const;
    void cleanupCache() const;
    SliceTexture* findSliceInCache(const std::string& hash, bool updateUsage) const;

    std::string getHash(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray, size_t levelOfDetail) const;
    std::string getHash(const VolumeBase* volume, tgt::plane pl, float samplingRate) const;

    Processor* owner_;

    mutable std::list<CacheEntry> slices_;
    size_t cacheSize_;

    // Used to store the last unfinished slice. Needed to delete it, since it is not stored in slices_.
    // Also stores the last requested slice, if cacheSize is zero. Needed for destruction.
    mutable CacheEntry tmpCache_;


    mutable SliceCreationBackgroundThread* currentBackgroundThread_;
    mutable boost::mutex cacheAccessMutex_; ///< mutex for synchronizing cache accesses
};


class SliceCreationBackgroundThread : public BackgroundThread {
    friend class SliceCache;
public:
    SliceCreationBackgroundThread(Processor* processor, const SliceCache* sliceCache,
           const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray,
            size_t levelOfDetail, const std::string& hash);

protected:
    virtual void threadMain();
    const SliceCache* sliceCache_;
    const VolumeBase* volume_;
    SliceAlignment alignment_;
    size_t sliceIndex_;
    int* shiftArray_;
    size_t levelOfDetail_;
    std::string hash_;
};


} // namespace voreen

#endif
