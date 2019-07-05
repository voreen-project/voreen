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

#ifndef VRN_VOLUMEMASK_H
#define VRN_VOLUMEMASK_H


#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/io/progressreporter.h"
#include "modules/bigdataimageprocessing/datastructures/lz4slicevolume.h"
#include "tgt/vector.h"
#include "tgt/filesystem.h"
#include "../util/tasktimelogger.h"
#include "surface.h"

#include <boost/iostreams/device/mapped_file.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <set>
#include <list>
#include <algorithm>

namespace voreen {

struct ScrapeIterationDescriptor {
    enum Dimension {
        DimX = 0,
        DimY = 1,
        DimZ = 2
    };
    enum Direction {
        Forwards = 1,
        Backwards = -1
    };

    ScrapeIterationDescriptor(Dimension dim, Direction dir)
        : dim(dim)
        , dir(dir)
        , deletedLastIteration(-1)
    {
    }
    tgt::ivec3 getNeightbor(tgt::svec3 pos) {
        tgt::ivec3 neighbor = pos;
        neighbor[dim] += dir;
        return neighbor;
    }
    std::string to_string() const {
        std::ostringstream out;
        out <<  " [" << dim << ", " << dir << "]";
        return out.str();
    }
    Dimension dim;
    Direction dir;
    size_t deletedLastIteration;
};


enum class VolumeMaskValue {
    FIXED_OBJECT = 0b11,
    OBJECT = 0b10,
    OUTSIDE_VOLUME = 0b01,
    BACKGROUND = 0b00,
};

#define VOLUME_MASK_STORAGE_BLOCK_SIZE_X 32
#define VOLUME_MASK_STORAGE_BLOCK_SIZE_Y 32
#define VOLUME_MASK_STORAGE_BLOCK_SIZE_Z 32
class VolumeMaskStorage {

public:
    VolumeMaskStorage(std::string filename, tgt::svec3 dimensions);

    ~VolumeMaskStorage();

    void set(const tgt::svec3& pos, VolumeMaskValue val);

    VolumeMaskValue get(const tgt::svec3& pos) const;

    static const size_t VOXELS_PER_BYTE = 4;

    tgt::svec3 dimensions_;
private:
    std::pair<size_t, char> indexAndSubindexFor(const tgt::svec3& pos) const {
        tgt::svec3 in_block_pos(pos.x % VOLUME_MASK_STORAGE_BLOCK_SIZE_X, pos.y % VOLUME_MASK_STORAGE_BLOCK_SIZE_Y, pos.z % VOLUME_MASK_STORAGE_BLOCK_SIZE_Z);
        tgt::svec3 block_pos(pos.x / VOLUME_MASK_STORAGE_BLOCK_SIZE_X, pos.y / VOLUME_MASK_STORAGE_BLOCK_SIZE_Y, pos.z / VOLUME_MASK_STORAGE_BLOCK_SIZE_Z);
        size_t in_block_linear = in_block_pos.x + VOLUME_MASK_STORAGE_BLOCK_SIZE_X*(in_block_pos.y + VOLUME_MASK_STORAGE_BLOCK_SIZE_Y * in_block_pos.z);
        size_t block_linear = block_pos.x + dimInBlocks_.x*(block_pos.y + dimInBlocks_.y * block_pos.z);
        size_t voxelpos = block_linear*VOLUME_MASK_STORAGE_BLOCK_SIZE_X*VOLUME_MASK_STORAGE_BLOCK_SIZE_Y*VOLUME_MASK_STORAGE_BLOCK_SIZE_Z + in_block_linear;
        return std::make_pair(voxelpos/VOXELS_PER_BYTE, voxelpos % VOXELS_PER_BYTE);
    }

    size_t data_index(size_t linearIndex) const;
    tgt::svec3 dimInBlocks_;

    const std::string filename_;
    boost::iostreams::mapped_file file_;
};

/**
 * Helper class that stores a binary volume in ram and provides methods to thin it.
 */
class VolumeMask {
public:

    enum ThinningAlgorithm {
        MA,
        CHEN,
        LEE,
        LEE_NO_LINE_PRESERVATION,
        IMPROVED,
        IMPROVED_NO_LINE_PRESERVATION,
    };

    VolumeMask(const LZ4SliceVolume<uint8_t>& vol, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress);
    VolumeMask(VolumeMask&& other);
    ~VolumeMask();

    void set(const tgt::svec3& pos, VolumeMaskValue val); //TODO should not be public, as it may invalidate surface information. Instead: use something like "fix" which can only act on foreground voxels
    VolumeMaskValue get(const tgt::ivec3& pos, VolumeMaskValue outsideVolumeValue = VolumeMaskValue::OBJECT) const;

    // Creates a VolumeRAM using the provided format and the contents of this mask.
    // Background voxels can be expected to be zero, while foreground voxels will
    // be > 0.
    voreen::VolumeRAM* toVolumeRAM(std::string format) const;

    // Remove layers of voxels from the surface of all objects using the provided algorithm
    template<ThinningAlgorithm T>
    void skeletonize(size_t maxIterations, ProgressReporter& progress);

    // Methods that thin a volume using Ma and Sonkas method
    void thinnerByMa(size_t& deleted);
    bool deletableMa(const tgt::svec3& pos) const;

    // Methods that thin a volume using the method of Chen et al.
    void thinnerByChen(size_t& deleted);
    bool deletableChen(const tgt::svec3& pos) const;

    bool isSurfaceVoxel(const tgt::svec3& pos) const;
    bool isLine(const tgt::ivec3& prev, const tgt::ivec3& start, int minLength) const;
    bool isEndVoxel(const tgt::svec3& pos, int lineVoxelChainMinLenght = 0) const;
    bool isEulerInvariantVoxel(const tgt::svec3& pos) const;
    bool isSimple(const tgt::svec3& pos) const;
    bool isSingleBorderVoxel(const tgt::svec3& pos) const;

    // Remove simple points from the volume. Use with care: This may not preserve
    // topology in all cases.
    void removeSimplePoints(size_t& deleted);

    template<class L>
    void thinnerByLee(ScrapeIterationDescriptor& scrapeDescriptor, size_t& numberOfDeletedVoxels, ProgressReporter& progress);

    // Improved method for volume thinning. Based on the method by Chen et al.
    // This method generally provides the "best" (i.e., least erroneous branches
    // and most medially positioned skeleton) results.
    template<class L>
    void scrape(ScrapeIterationDescriptor& scrapeDescriptor, size_t& numberOfDeletedVoxels, ProgressReporter& progress);

    const tgt::svec3& getDimensions() const;
    size_t numOriginalForegroundVoxels_;

private:
    uint64_t toLinearPos(const tgt::svec3& pos);
    tgt::svec3 fromLinearPos(uint64_t pos);

    VolumeMaskStorage data_;
    tgt::vec3 spacing_;

    StoredSurface surfaceFile_;

    static const std::string loggerCat_;
};

struct NoFixedForeground {
    void advance() {
    }
    bool isForeground(size_t x, size_t y) {
        return false;
    }
};

// --------------------------------------------------------------------------------------------------
// Implementation -----------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------



template<class L>
void VolumeMask::scrape(ScrapeIterationDescriptor& scrapeDescriptor, size_t& numberOfDeletedVoxels, ProgressReporter& progress) {
    SurfaceReader surfaceReader(surfaceFile_);
    SurfaceBuilder builder;

    SurfaceSlices<4> surface;

    std::vector<tgt::svec3> to_delete_prev_prev;
    std::vector<tgt::svec3> to_delete_prev;
    std::vector<tgt::svec3> to_delete_current;

    size_t deleted_this_it = 0;
    auto deleteAll = [this, &deleted_this_it] (std::vector<tgt::svec3>& voxels, SurfaceSlices<4>& surface) {
        for(const auto& pos : voxels) {
            if(L::deletableScraping(*this, pos)) {
                set(pos, VolumeMaskValue::BACKGROUND);
                ++deleted_this_it;

                // Add neighbors to active surface (as their status might have changed now)
                {
                    int z = -1;
                    for(int y = -1; y != 2; ++y) {
                        for(int x = -1; x != 2; ++x) {
                            tgt::ivec3 p = tgt::ivec3(pos) + tgt::ivec3(x, y, z);
                            if(get(p, VolumeMaskValue::BACKGROUND) == VolumeMaskValue::OBJECT) {
                                surface.m<3>().push_back(toLinearPos(tgt::svec3(p)));
                            }
                        }
                    }
                }

                {
                    int z = 0;
                    for(int y = -1; y != 2; ++y) {
                        for(int x = -1; x != 2; ++x) {
                            if(x != 0 || y != 0) {
                                tgt::ivec3 p = tgt::ivec3(pos) + tgt::ivec3(x, y, z);
                                if(get(p, VolumeMaskValue::BACKGROUND) == VolumeMaskValue::OBJECT) {
                                    surface.m<2>().push_back(toLinearPos(tgt::svec3(p)));
                                }
                            }
                        }
                    }
                }

                {
                    int z = 1;
                    for(int y = -1; y != 2; ++y) {
                        for(int x = -1; x != 2; ++x) {
                            tgt::ivec3 p = tgt::ivec3(pos) + tgt::ivec3(x, y, z);
                            if(get(p, VolumeMaskValue::BACKGROUND) == VolumeMaskValue::OBJECT) {
                                surface.m<1>().push_back(toLinearPos(tgt::svec3(p)));
                            }
                        }
                    }
                }

            } else {
                //Voxel is potentially deletable, but was not, so it will (or may at least) become inactive for now
                // => we do not: "surface.m<2>().push_back(toLinearPos(pos));"
            }
        }
    };

    size_t z = 0;
    size_t currentVoxel = 0;
    size_t updateInterval = 10000;
    for(uint64_t linearPos = -1; surfaceReader.read(linearPos); ++currentVoxel) {
        tgtAssert(linearPos < tgt::hmul(getDimensions()), "Invalid linear pos read from file");
        if((currentVoxel % updateInterval) == 0) {
            progress.setProgress(static_cast<float>(numberOfDeletedVoxels)/numOriginalForegroundVoxels_);
        }

        tgt::svec3 pos = fromLinearPos(linearPos);

        if(get(pos, VolumeMaskValue::BACKGROUND) != VolumeMaskValue::OBJECT) {
            // Voxel might have been added to surface before deletion (i.e. after deletion of a voxel in a previous slice)
            // In this case we can skip the voxel.
            continue;
        }

        while(pos.z != z) {
            tgtAssert(pos.z > z, "pos too small");

            // 1
            deleteAll(to_delete_prev_prev, surface);
            to_delete_prev_prev.clear();

            std::swap(to_delete_prev_prev, to_delete_prev);
            std::swap(to_delete_prev, to_delete_current);
            // 2
            surface.advance(builder);

            ++z;
        }

        // 3
        if(get(scrapeDescriptor.getNeightbor(pos), VolumeMaskValue::OBJECT) == VolumeMaskValue::BACKGROUND) {
            if (isEulerInvariantVoxel(pos)) {
                to_delete_current.push_back(pos);
            }
        } else {
            // Active surface voxel is not valid for current scraping direction
            surface.m<0>().push_back(linearPos);
        }
    }
    for(int i=0; i<5; ++i) {
        deleteAll(to_delete_prev_prev, surface);
        to_delete_prev_prev.clear();

        std::swap(to_delete_prev_prev, to_delete_prev);
        std::swap(to_delete_prev, to_delete_current);

        surface.advance(builder);
    }
    tgtAssert(to_delete_prev_prev.empty() && to_delete_prev.empty() && to_delete_current.empty(), "Deletion unfinished");
    tgtAssert(surface.m<0>().empty() && surface.m<1>().empty() && surface.m<2>().empty() && surface.m<3>().empty(), "Writing surface back unfinished");

    surfaceFile_ = std::move(builder).finalize();


    progress.setProgress(static_cast<float>(numberOfDeletedVoxels)/numOriginalForegroundVoxels_);
    numberOfDeletedVoxels += deleted_this_it;
    scrapeDescriptor.deletedLastIteration = deleted_this_it;
}

struct LinePreservingScraping {
    static bool deletableScraping(const VolumeMask& m, const tgt::svec3& pos) {
        return m.get(pos, VolumeMaskValue::BACKGROUND) != VolumeMaskValue::BACKGROUND
            && !m.isEndVoxel(pos, 0 /* Parameter! One might want to consider other values*/)
            && m.isEulerInvariantVoxel(pos)
            && m.isSimple(pos)
        ;
    }
};

struct NoLinePreservingScraping {
    static bool deletableScraping(const VolumeMask& m, const tgt::svec3& pos) {
        return m.get(pos, VolumeMaskValue::BACKGROUND) != VolumeMaskValue::BACKGROUND
            && m.isEulerInvariantVoxel(pos)
            && m.isSimple(pos)
        ;
    }
};


struct ScrapeIterationQueueElement {
    ScrapeIterationQueueElement(ScrapeIterationDescriptor descriptor, float initialProgress)
        : descriptor(descriptor)
        , accumulatedProgress(initialProgress)
    {
    }

    ScrapeIterationDescriptor descriptor;
    float accumulatedProgress;
};

inline bool operator <(const ScrapeIterationQueueElement& lhs, const ScrapeIterationQueueElement& rhs) {
    return lhs.accumulatedProgress < rhs.accumulatedProgress;
}

class ScrapeIterationCoordinator {
public:
    ScrapeIterationCoordinator(tgt::vec3 volumeSpacing)
        : volumeSpacing_(volumeSpacing)
        , nextIterations_()
    {
        nextIterations_.push_back(ScrapeIterationQueueElement(ScrapeIterationDescriptor(ScrapeIterationDescriptor::DimX, ScrapeIterationDescriptor::Forwards), 0.0f));
        nextIterations_.push_back(ScrapeIterationQueueElement(ScrapeIterationDescriptor(ScrapeIterationDescriptor::DimX, ScrapeIterationDescriptor::Backwards), 0.0f));
        nextIterations_.push_back(ScrapeIterationQueueElement(ScrapeIterationDescriptor(ScrapeIterationDescriptor::DimY, ScrapeIterationDescriptor::Forwards), 0.0f));
        nextIterations_.push_back(ScrapeIterationQueueElement(ScrapeIterationDescriptor(ScrapeIterationDescriptor::DimY, ScrapeIterationDescriptor::Backwards), 0.0f));
        nextIterations_.push_back(ScrapeIterationQueueElement(ScrapeIterationDescriptor(ScrapeIterationDescriptor::DimZ, ScrapeIterationDescriptor::Forwards), 0.0f));
        nextIterations_.push_back(ScrapeIterationQueueElement(ScrapeIterationDescriptor(ScrapeIterationDescriptor::DimZ, ScrapeIterationDescriptor::Backwards), 0.0f));
    }

    ScrapeIterationDescriptor& next() {
        auto nextIt = std::min_element(nextIterations_.begin(), nextIterations_.end()); //Returns iterator to the FIRST smallest element
        nextIterations_.push_back(*nextIt);
        nextIterations_.erase(nextIt);

        auto& next = nextIterations_.back();
        next.accumulatedProgress += volumeSpacing_[next.descriptor.dim];
        return next.descriptor;
    }

    bool finishedOnAllAxes() const {
        for(const auto& elm: nextIterations_) {
            if(elm.descriptor.deletedLastIteration != 0) {
                return false;
            }
        }
        return true;
    }

private:
    tgt::vec3 volumeSpacing_;
    std::list<ScrapeIterationQueueElement> nextIterations_;
};

template<VolumeMask::ThinningAlgorithm T>
void VolumeMask::skeletonize(size_t maxIterations, ProgressReporter& progress) {
    TaskTimeLogger _("Skeletonize VolumeMask", tgt::Info);
    size_t iteration = 0;
    size_t deleted = 0;
    size_t deleted_prev_it = -1;
    ScrapeIterationCoordinator scrapeCoordinator(spacing_);
    bool exitConditionReached = false;
    while(!exitConditionReached && maxIterations > iteration) {
        deleted_prev_it = deleted;
        std::string it_info = "";
        auto& next = scrapeCoordinator.next();
        switch(T) {
            case VolumeMask::MA:
                thinnerByMa(deleted);
                break;
            case VolumeMask::CHEN:
                thinnerByChen(deleted);
                break;
            case VolumeMask::LEE:
                it_info = next.to_string();
                thinnerByLee<LinePreservingScraping>(next, deleted, progress);
                break;
            case VolumeMask::LEE_NO_LINE_PRESERVATION:
                it_info = next.to_string();
                thinnerByLee<NoLinePreservingScraping>(next, deleted, progress);
                break;
            case VolumeMask::IMPROVED:
                it_info = next.to_string();
                scrape<LinePreservingScraping>(next, deleted, progress);
                break;
            case VolumeMask::IMPROVED_NO_LINE_PRESERVATION:
                it_info = next.to_string();
                scrape<NoLinePreservingScraping>(next, deleted, progress);
                break;
        }
        size_t deleted_this_it = deleted - deleted_prev_it;
        LINFO("Thinning iteration " + std::to_string(iteration) + it_info + " deleted " + std::to_string(deleted_this_it) + " voxels.");
        ++iteration;

        switch(T) {
            case VolumeMask::MA:
            case VolumeMask::CHEN:
                if(deleted == deleted_prev_it) {
                    exitConditionReached = true;
                }
                break;
            case VolumeMask::LEE:
            case VolumeMask::LEE_NO_LINE_PRESERVATION:
            case VolumeMask::IMPROVED:
            case VolumeMask::IMPROVED_NO_LINE_PRESERVATION:
                if(scrapeCoordinator.finishedOnAllAxes()) {
                    exitConditionReached = true;
                }
                break;
        }
    }
    progress.setProgress(1.0f);
    LINFO("Thinning iteration finished and deleted " + std::to_string(deleted) + " voxels.");
}

template<class L>
void VolumeMask::thinnerByLee(ScrapeIterationDescriptor& scrapeDescriptor, size_t& numberOfDeletedVoxels, ProgressReporter& progress) {
    std::vector<tgt::svec3> marked;
    tgt::ivec3 dim = getDimensions();
    for(int z = 0; z<dim.z; ++z) {
        for(int y = 0; y<dim.y; ++y) {
            for(int x = 0; x<dim.x; ++x) {
                tgt::svec3 pos = tgt::svec3(x,y,z);

                if(get(pos) == VolumeMaskValue::OBJECT
                        && L::deletableScraping(*this, pos)
                        && get(scrapeDescriptor.getNeightbor(pos), VolumeMaskValue::OBJECT) == VolumeMaskValue::BACKGROUND) {
                    marked.push_back(pos);
                }
            }
        }
    }

    size_t deletedThisIt = 0;
    for(auto& pos : marked) {
        if(isSimple(pos)) {
            set(pos, VolumeMaskValue::BACKGROUND);
            ++deletedThisIt;
        }
    }
    numberOfDeletedVoxels += deletedThisIt;
    scrapeDescriptor.deletedLastIteration = deletedThisIt;
}


} // namespace voreen

#endif // VRN_VOLUMETHINNING_H
