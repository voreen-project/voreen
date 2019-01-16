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

#include "volumemask.h"

#include "deletetemplates.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "tgt/filesystem.h"

#include <csignal>

namespace {

template <typename T>
int8_t sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
int8_t readRelative(const voreen::VolumeMask* vol, const tgt::svec3& p, int dx, int dy, int dz) {
    return sgn(vol->get(tgt::ivec3(p) + tgt::ivec3(dx, dy, dz)));
}
int8_t readRelative(const voreen::VolumeMask* vol, const tgt::ivec3& p, int dx, int dy, int dz) {
    return sgn(vol->get(p + tgt::ivec3(dx, dy, dz)));
}
bool fitsTemplate(const voreen::VolumeMask* vol, const tgt::svec3& pos, size_t templateNum) {
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                if(readRelative(vol, pos, dx, dy, dz)*deleteTemplates[templateNum][dz+1][dy+1][dx+1] < 0) { //does not fit template
                    return false;
                }
            }
        }
    }
    return fitsExpandedTemplate(vol, pos, templateNum);
    //return true;
}
bool isObjectPoint(const voreen::VolumeMask* vol, const tgt::svec3& p, int dx, int dy, int dz) {
    return vol->get(tgt::ivec3(p.x+dx,p.y+dy,p.z+dz)) != voreen::VolumeMaskValue::BACKGROUND;
}
bool isObjectPoint(const voreen::VolumeMask* vol, const tgt::ivec3& p, int dx, int dy, int dz) {
    return vol->get(tgt::ivec3(p.x+dx,p.y+dy,p.z+dz)) != voreen::VolumeMaskValue::BACKGROUND;
}
bool isTailPoint(const voreen::VolumeMask* vol, const tgt::svec3& pos) {
    size_t numNeighbours = 0;
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                if(isObjectPoint(vol, pos, dx, dy, dz)) {
                    ++numNeighbours;
                }
            }
        }
    }
    numNeighbours -= 1; //the center voxel does not count
    if(numNeighbours == 1) {
        return true;
    }
    if(numNeighbours != 2) {
        return false;
    }
    return (isObjectPoint(vol, pos, 0, 0,-1) && (isObjectPoint(vol, pos, 1, 0, 0) ^ isObjectPoint(vol, pos, 0, 1, 0)))
        || (isObjectPoint(vol, pos, 0,-1, 0) && (isObjectPoint(vol, pos, 1, 0, 0) ^ isObjectPoint(vol, pos, 0, 0, 1)))
        || (isObjectPoint(vol, pos,-1, 0, 0) && (isObjectPoint(vol, pos, 0, 1, 0) ^ isObjectPoint(vol, pos, 0, 0, 1)));
}

bool isNearBorderObjectPoint(const voreen::VolumeMask* vol, const tgt::svec3& pos) {
    if(vol->get(pos)==voreen::VolumeMaskValue::BACKGROUND) {
        return false;
    }
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                if(!isObjectPoint(vol, pos, dx, dy, dz)) {
                    return true;
                }
            }
        }
    }
    return false;
}

template<int8_t Z, int8_t Y, int8_t X>
void minimize(uint8_t (&neighborhood)[3][3][3], uint8_t& val) {
    if(   X < 0 || X >= 3
       || Y < 0 || Y >= 3
       || Z < 0 || Z >= 3) {
        return;
    }
    val = std::min(val, neighborhood[Z][Y][X]);
}
template<int8_t Z, int8_t Y, int8_t X>
uint8_t min(uint8_t (&neighborhood)[3][3][3]) {
    uint8_t id = 255;
    minimize<Z-1,Y-1,X-1>(neighborhood, id);
    minimize<Z-1,Y-1,X  >(neighborhood, id);
    minimize<Z-1,Y-1,X+1>(neighborhood, id);
    minimize<Z-1,Y  ,X-1>(neighborhood, id);
    minimize<Z-1,Y  ,X  >(neighborhood, id);
    minimize<Z-1,Y  ,X+1>(neighborhood, id);
    minimize<Z-1,Y+1,X-1>(neighborhood, id);
    minimize<Z-1,Y+1,X  >(neighborhood, id);
    minimize<Z-1,Y+1,X+1>(neighborhood, id);
    minimize<Z  ,Y-1,X-1>(neighborhood, id);
    minimize<Z  ,Y-1,X  >(neighborhood, id);
    minimize<Z  ,Y-1,X+1>(neighborhood, id);
    minimize<Z  ,Y  ,X-1>(neighborhood, id);
    minimize<Z  ,Y  ,X  >(neighborhood, id);
    minimize<Z  ,Y  ,X+1>(neighborhood, id);
    minimize<Z  ,Y+1,X-1>(neighborhood, id);
    minimize<Z  ,Y+1,X  >(neighborhood, id);
    minimize<Z  ,Y+1,X+1>(neighborhood, id);
    minimize<Z+1,Y-1,X-1>(neighborhood, id);
    minimize<Z+1,Y-1,X  >(neighborhood, id);
    minimize<Z+1,Y-1,X+1>(neighborhood, id);
    minimize<Z+1,Y  ,X-1>(neighborhood, id);
    minimize<Z+1,Y  ,X  >(neighborhood, id);
    minimize<Z+1,Y  ,X+1>(neighborhood, id);
    minimize<Z+1,Y+1,X-1>(neighborhood, id);
    minimize<Z+1,Y+1,X  >(neighborhood, id);
    minimize<Z+1,Y+1,X+1>(neighborhood, id);
    return id;
}

template<int8_t Z, int8_t Y, int8_t X>
void linkNeighbors(uint8_t (&neighborhood)[3][3][3]) {
    if(X == 1 && Y == 1 && Z == 1){
        return;
    }
    if(neighborhood[Z][Y][X] == 255) {
        return;
    }
    uint8_t minID = min<Z, Y, X>(neighborhood);
    neighborhood[Z][Y][X] = minID;
}

int8_t D_EULER_CHARACTERISTIC [128] = {
     1, -1, -1,  1, -3, -1, -1,  1,
    -1,  1,  1, -1,  3,  1,  1, -1,
    -3, -1,  3,  1,  1, -1,  3,  1,
    -1,  1,  1, -1,  3,  1,  1, -1,
    -3,  3, -1,  1,  1,  3, -1,  1,
    -1,  1,  1, -1,  3,  1,  1, -1,
     1,  3,  3,  1,  5,  3,  3,  1,
    -1,  1,  1, -1,  3,  1,  1, -1,
    -7, -1, -1,  1, -3, -1, -1,  1,
    -1,  1,  1, -1,  3,  1,  1, -1,
    -3, -1,  3,  1,  1, -1,  3,  1,
    -1,  1,  1, -1,  3,  1,  1, -1,
    -3,  3, -1,  1,  1,  3, -1,  1,
    -1,  1,  1, -1,  3,  1,  1, -1,
     1,  3,  3,  1,  5,  3,  3,  1,
    -1,  1,  1, -1,  3,  1,  1, -1
};
int8_t OCTANT_VOXEL_INDEX [8][7] {
    {  0,  1,  3,  4,  9, 10, 12},
    {  2,  1,  5,  4, 11, 10, 13},
    {  6,  7,  3,  4, 14, 15, 12},
    {  8,  7,  5,  4, 16, 15, 13},
    { 17, 18, 20, 21,  9, 10, 12},
    { 19, 18, 22, 21, 11, 10, 13},
    { 23, 24, 20, 21, 14, 15, 12},
    { 25, 24, 22, 21, 16, 15, 13}
};


} // anonymous namespace

namespace voreen {

const std::string VolumeMask::loggerCat_("voreen.vesseltopology.volumemask");

VolumeMask::VolumeMask(VolumeMask&& other)
    : data_(std::move(other.data_))
    , spacing_(other.spacing_)
    , surfaceFile_(std::move(other.surfaceFile_))
{
}

VolumeMask::VolumeMask(const LZ4SliceVolume<uint8_t>& vol, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress)
    : data_(VoreenApplication::app()->getUniqueTmpFilePath(), vol.getDimensions())
    , surfaceFile_("", 0)
    , numOriginalForegroundVoxels_(0)
    , spacing_(vol.getMetaData().getSpacing())
{
    TaskTimeLogger _("Create VolumeMask", tgt::Info);
    SubtaskProgressReporterCollection<2> subtaskReporters(progress);
    tgt::svec3 dimensions = vol.getDimensions();

    tgtAssert(!sampleMask || dimensions == sampleMask->getDimensions(), "Sample mask dimension mismatch");

    auto calc_num_chunks = [] (size_t dim, size_t size) {
        if((dim % size) == 0) {
            return dim/size;
        } else {
            return dim/size+1;
        }
    };

    auto calc_chunk_size = [&] (size_t chunk, size_t num_chunks, size_t max_chunk_size, size_t dimensions) {
        size_t size = chunk != num_chunks - 1 || (dimensions % max_chunk_size == 0) ? max_chunk_size : dimensions % max_chunk_size; //TODO check
        tgtAssert(size > 0, "Invalid chunk size");
        return size;
    };

    const size_t max_chunk_size_x = VOLUME_MASK_STORAGE_BLOCK_SIZE_X;
    const size_t num_chunks_x = calc_num_chunks(dimensions.x, max_chunk_size_x);
    const size_t max_chunk_size_y = VOLUME_MASK_STORAGE_BLOCK_SIZE_Y;
    const size_t num_chunks_y = calc_num_chunks(dimensions.y, max_chunk_size_y);
    const size_t max_chunk_size_z = VOLUME_MASK_STORAGE_BLOCK_SIZE_Z;
    const size_t num_chunks_z = calc_num_chunks(dimensions.z, max_chunk_size_z);

    // Initialize volumetric data (data_) in chunks_z (more friendly to disk layout)
    for(size_t chunk_z = 0; chunk_z < num_chunks_z; ++chunk_z) {
        subtaskReporters.get<0>().setProgress(static_cast<float>(chunk_z)/num_chunks_z);

        const size_t chunk_size_z = calc_chunk_size(chunk_z, num_chunks_z, max_chunk_size_z, dimensions.z);

        size_t slab_begin = chunk_z*max_chunk_size_z;
        size_t slab_end = slab_begin + chunk_size_z;
        auto volSlab = vol.loadSlab(slab_begin, slab_end);
        boost::optional<VolumeAtomic<uint8_t>> sampleMaskSlab = sampleMask ? boost::optional<VolumeAtomic<uint8_t>>(sampleMask->loadSlab(slab_begin, slab_end)) : boost::none;

        for(size_t chunk_y = 0; chunk_y < num_chunks_y; ++chunk_y) {
            const size_t chunk_size_y = calc_chunk_size(chunk_y, num_chunks_y, max_chunk_size_y, dimensions.y);

            for(size_t chunk_x = 0; chunk_x < num_chunks_x; ++chunk_x) {
                const size_t chunk_size_x = calc_chunk_size(chunk_x, num_chunks_x, max_chunk_size_x, dimensions.x);

                for(size_t dz = 0; dz<chunk_size_z; ++dz) {
                    size_t z = chunk_z*max_chunk_size_z + dz;

                    for(size_t dy = 0; dy<chunk_size_y; ++dy) {
                        size_t y = chunk_y*max_chunk_size_y + dy;

                        for(size_t dx = 0; dx<chunk_size_x; ++dx) {
                            size_t x = chunk_x*max_chunk_size_x + dx;

                            const tgt::svec3 p(x,y,z);
                            const tgt::svec3 slabPos(x,y, dz);
                            VolumeMaskValue val;
                            if(sampleMaskSlab && sampleMaskSlab->voxel(slabPos) == 0) {
                                val = VolumeMaskValue::OUTSIDE_VOLUME;
                            } else if(volSlab.getVoxelNormalized(slabPos) > 0) {
                                val = VolumeMaskValue::OBJECT;
                                ++numOriginalForegroundVoxels_;
                            } else {
                                val = VolumeMaskValue::BACKGROUND;
                            }
                            // NOTE! This is a shortcut that only works under the assumption that newly created
                            // files are filled with zeros, which I'm PRETTY sure of, but not 100% certain.
                            if(val != VolumeMaskValue::BACKGROUND) {
                                data_.set(p, val);
                            } else {
                                tgtAssert(data_.get(p) == VolumeMaskValue::BACKGROUND, "File was not zeroed before filling");
                            }
                        }
                    }
                }
            }
        }
    }

    // Initialize surface
    SurfaceBuilder builder;
    for(size_t z = 0; z<dimensions.z; ++z) {
        subtaskReporters.get<1>().setProgress(static_cast<float>(z)/dimensions.z);
        for(size_t y = 0; y<dimensions.y; ++y) {
            for(size_t x = 0; x<dimensions.x; ++x) {
                const tgt::svec3 p(x,y,z);
                if(get(p, VolumeMaskValue::BACKGROUND) == VolumeMaskValue::OBJECT && isSurfaceVoxel(p)) {
                    builder.push(toLinearPos(p));
                }
            }
        }
    }
    surfaceFile_ = std::move(builder).finalize();
    progress.setProgress(1.0f);
}


VolumeMask::~VolumeMask() {
    tgt::FileSystem::deleteFile(surfaceFile_.filename_);
}

void VolumeMask::set(const tgt::svec3& pos, VolumeMaskValue val) {
    tgtAssert(val == VolumeMaskValue::FIXED_OBJECT || data_.get(pos) != VolumeMaskValue::FIXED_OBJECT, "tried to change fixed object voxel");
    data_.set(pos, val);
}
VolumeMaskValue VolumeMask::get(const tgt::ivec3& pos, VolumeMaskValue outsideVolumeValue) const {
    auto dim = getDimensions();
    if(     pos.x < 0 || pos.x >= static_cast<int>(dim.x)
         || pos.y < 0 || pos.y >= static_cast<int>(dim.y)
         || pos.z < 0 || pos.z >= static_cast<int>(dim.z)) {
        return outsideVolumeValue;
    } else {
        VolumeMaskValue val = data_.get(tgt::svec3(pos));
        if(val == VolumeMaskValue::OUTSIDE_VOLUME) {
            // Value might be masked
            return outsideVolumeValue;
        } else {
            return val;
        }
    }
}
voreen::VolumeRAM* VolumeMask::toVolumeRAM(std::string format) const {
    voreen::VolumeFactory vf;
    auto dim = getDimensions();
    voreen::VolumeRAM* output = vf.create(format, dim);
    for(size_t z = 0; z<dim.z; ++z) {
        for(size_t y = 0; y<dim.y; ++y) {
            for(size_t x = 0; x<dim.x; ++x) {
                const tgt::svec3 p(x,y,z);
                VolumeMaskValue val = get(p, VolumeMaskValue::BACKGROUND);
                tgtAssert( val==VolumeMaskValue::OBJECT || val==VolumeMaskValue::FIXED_OBJECT || val==VolumeMaskValue::BACKGROUND, "Invalid value in volume mask: ");
                if(val == VolumeMaskValue::BACKGROUND) {
                    output->setVoxelNormalized(0, p);
                } else /*val == OBJECT || FIXED_OBJECT*/ {
                    output->setVoxelNormalized(1, p);
                }
            }
        }
    }
    return output;
}
bool VolumeMask::deletableMa(const tgt::svec3& pos) const {
    if(isTailPoint(this, pos)) {
        return false;
    }
    for(size_t templateNum = 0; templateNum < sizeof(deleteTemplates)/sizeof(deleteTemplates[0]); ++templateNum) {
        if(fitsTemplate(this, pos, templateNum)) {
            //LINFO(std::to_string(pos.x) + "," + std::to_string(pos.y) + "," + std::to_string(pos.z) + " fits template " + templateNames[templateNum]);
            return true;
        }
    }
    return false;
}

void VolumeMask::thinnerByMa(size_t& numberOfDeletedVoxels) {
    numberOfDeletedVoxels = 0;
    std::list<tgt::svec3> marked;
    auto dim = getDimensions();
    for(int z = dim.z -1 ; z>=0; --z) {
        for(int y = dim.y -1 ; y>=0; --y) {
            for(int x = dim.x -1 ; x>=0; --x) {
                tgt::svec3 pos = tgt::svec3(x,y,z);
                if(get(pos, VolumeMaskValue::BACKGROUND) == VolumeMaskValue::OBJECT && isNearBorderObjectPoint(this, pos)) {
                //if(isNearBorderObjectPoint(this, pos) && isEulerInvariantVoxel(pos)) {
                    marked.push_back(pos);
                }
            }
        }
    }
    size_t deletedThisIt;
    do {
        deletedThisIt = 0;
        for(auto pos = marked.begin(); pos != marked.end();) {
            if(deletableMa(*pos)) {
                set(*pos, VolumeMaskValue::BACKGROUND);
                ++deletedThisIt;
                pos = marked.erase(pos);
            } else {
                ++pos;
            }
        }
        numberOfDeletedVoxels += deletedThisIt;
    } while(deletedThisIt > 0);
}

bool VolumeMask::isSurfaceVoxel(const tgt::svec3& pos) const {
    if(!isObjectPoint(this, pos, 0, 0, 0)) return false;
    if(!isObjectPoint(this, pos, 1, 0, 0)) return true;
    if(!isObjectPoint(this, pos,-1, 0, 0)) return true;
    if(!isObjectPoint(this, pos, 0, 1, 0)) return true;
    if(!isObjectPoint(this, pos, 0,-1, 0)) return true;
    if(!isObjectPoint(this, pos, 0, 0, 1)) return true;
    if(!isObjectPoint(this, pos, 0, 0,-1)) return true;
    return false;
}

bool VolumeMask::isLine(const tgt::ivec3& prev, const tgt::ivec3& start, int minLength) const {
    if(minLength == 0) {
        return true;
    }
    tgt::ivec3 neighbor;
    size_t numObjects = 0;
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                tgt::ivec3 p = start + tgt::ivec3(dx, dy, dz);
                if(get(p) != VolumeMaskValue::BACKGROUND && p != prev && (dx != 0 || dy != 0 || dz != 0)) {
                    ++numObjects;
                    neighbor = tgt::ivec3(p);
                }
            }
        }
    }
    return (numObjects == 1) && isLine(start, neighbor, minLength - 1);
}
bool VolumeMask::isEndVoxel(const tgt::svec3& pos, int lineVoxelChainMinLength) const {
    return isLine(pos, pos, lineVoxelChainMinLength + 1);
}
bool VolumeMask::isEulerInvariantVoxel(const tgt::svec3& pos) const {
    bool neighborhood[26];
    int i = 0;
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                if(dx != 0 || dy != 0 ||  dz!= 0) {
                    neighborhood[i++] = get(tgt::ivec3(pos.x + dx, pos.y + dy, pos.z + dz)) != VolumeMaskValue::BACKGROUND;
                }
            }
        }
    }
    int eulerVal = 0;
    for(int octant = 0; octant < 8; ++octant) {
        int index = 0;
        for(int voxel = 0; voxel < 7; ++voxel) {
            index <<= 1;
            index |= neighborhood[OCTANT_VOXEL_INDEX[octant][voxel]];
        }
        eulerVal += D_EULER_CHARACTERISTIC[index];
    }
    return eulerVal == 0;
}
bool VolumeMask::isSimple(const tgt::svec3& pos) const {
    if(get(pos)==VolumeMaskValue::BACKGROUND) {
        return false;
    }
    uint8_t neighborhood[3][3][3];
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                neighborhood[dz+1][dy+1][dx+1] = 255;
            }
        }
    }
    uint8_t id = 0;
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                if((dx != 0 || dy != 0 ||  dz!= 0) && isObjectPoint(this, pos, dx, dy, dz)) {
                    neighborhood[dz+1][dy+1][dx+1] = id++;
                } else {
                    neighborhood[dz+1][dy+1][dx+1] = 255;
                }
            }
        }
    }
    //if(id == 1 /* exactly one neighbor */) {
    //    return false;
    //}
    tgtAssert(id <= 28, "impossible");
    uint8_t maxID = -2;
    uint8_t prevMaxID = -1;
    while(maxID != prevMaxID) {
        linkNeighbors<0, 0, 0>(neighborhood);
        linkNeighbors<0, 0, 1>(neighborhood);
        linkNeighbors<0, 0, 2>(neighborhood);
        linkNeighbors<0, 1, 0>(neighborhood);
        linkNeighbors<0, 1, 1>(neighborhood);
        linkNeighbors<0, 1, 2>(neighborhood);
        linkNeighbors<0, 2, 0>(neighborhood);
        linkNeighbors<0, 2, 1>(neighborhood);
        linkNeighbors<0, 2, 2>(neighborhood);
        linkNeighbors<1, 0, 0>(neighborhood);
        linkNeighbors<1, 0, 1>(neighborhood);
        linkNeighbors<1, 0, 2>(neighborhood);
        linkNeighbors<1, 1, 0>(neighborhood);
        //linkNeighbors<1, 1, 1>(neighborhood);
        linkNeighbors<1, 1, 2>(neighborhood);
        linkNeighbors<1, 2, 0>(neighborhood);
        linkNeighbors<1, 2, 1>(neighborhood);
        linkNeighbors<1, 2, 2>(neighborhood);
        linkNeighbors<2, 0, 0>(neighborhood);
        linkNeighbors<2, 0, 1>(neighborhood);
        linkNeighbors<2, 0, 2>(neighborhood);
        linkNeighbors<2, 1, 0>(neighborhood);
        linkNeighbors<2, 1, 1>(neighborhood);
        linkNeighbors<2, 1, 2>(neighborhood);
        linkNeighbors<2, 2, 0>(neighborhood);
        linkNeighbors<2, 2, 1>(neighborhood);
        linkNeighbors<2, 2, 2>(neighborhood);

        prevMaxID = maxID;
        maxID = min<1, 1, 1>(neighborhood);
        for(int dz = -1; dz<=1; ++dz) {
            for(int dy = -1; dy<=1; ++dy) {
                for(int dx = -1; dx<=1; ++dx) {
                    uint8_t id = neighborhood[dz+1][dy+1][dx+1];
                    if(id != 255 && id > maxID) {
                        maxID = id;
                    }
                }
            }
        }
    }
    if(maxID == 0) {
        return true;
    } else {
        return false;
    }
}

bool VolumeMask::isSingleBorderVoxel(const tgt::svec3& pos) const {
    int numNeighbours = 0;
    for(int dz = -1; dz<=1; ++dz) {
        for(int dy = -1; dy<=1; ++dy) {
            for(int dx = -1; dx<=1; ++dx) {
                if(get(tgt::ivec3(pos.x+dx, pos.y+dy, pos.z+dz), VolumeMaskValue::BACKGROUND) != VolumeMaskValue::BACKGROUND) {
                    ++numNeighbours;
                }
            }
        }
    }
    return (numNeighbours - 1 /* minus one for pos itself */) < 1;
}

bool VolumeMask::deletableChen(const tgt::svec3& pos) const {
    return isObjectPoint(this, pos, 0, 0, 0)
        //&& isSurfaceVoxel(pos)
        && !isEndVoxel(pos, 0) // Parameter! One might want to consider other values
        && isEulerInvariantVoxel(pos)
        && isSimple(pos);
}

void VolumeMask::thinnerByChen(size_t& numberOfDeletedVoxels) {
    numberOfDeletedVoxels = 0;
    std::list<tgt::svec3> marked;
    auto dim = getDimensions();
    for(int z = dim.z -1 ; z>=0; --z) {
        for(int y = dim.y -1 ; y>=0; --y) {
            for(int x = dim.x -1 ; x>=0; --x) {
                tgt::svec3 pos = tgt::svec3(x,y,z);
                //if(isNearBorderObjectPoint(this, pos)) {
                if(get(pos, VolumeMaskValue::BACKGROUND) == VolumeMaskValue::OBJECT && isSurfaceVoxel(pos)) {
                    marked.push_back(pos);
                }
            }
        }
    }

    size_t deletedThisIt;
    deletedThisIt = 0;
    for(auto pos = marked.begin(); pos != marked.end();) {
        if(deletableChen(*pos)) {
            set(*pos, VolumeMaskValue::BACKGROUND);
            ++deletedThisIt;
            pos = marked.erase(pos);
        } else {
            ++pos;
        }
    }
    numberOfDeletedVoxels += deletedThisIt;
}
const tgt::svec3& VolumeMask::getDimensions() const {
    return data_.dimensions_;
}

uint64_t VolumeMask::toLinearPos(const tgt::svec3& pos) {
    return pos.x + data_.dimensions_.x*(pos.y + data_.dimensions_.y*pos.z);
}
tgt::svec3 VolumeMask::fromLinearPos(uint64_t pos) {
    tgt::svec3 p;
    p.x = pos % data_.dimensions_.x;
    pos /= data_.dimensions_.x;
    p.y = pos % data_.dimensions_.y;
    pos /= data_.dimensions_.y;
    p.z = pos;
    return p;
}

/// VolumeMaskStorage ----------------------------------------------------------
VolumeMaskStorage::VolumeMaskStorage(std::string filename, tgt::svec3 dimensions)
    : file_()
    , dimensions_(dimensions)
    , filename_(filename)
    , dimInBlocks_(
            (dimensions.x / VOLUME_MASK_STORAGE_BLOCK_SIZE_X) + (dimensions.x % VOLUME_MASK_STORAGE_BLOCK_SIZE_X == 0 ? 0 : 1),
            (dimensions.y / VOLUME_MASK_STORAGE_BLOCK_SIZE_Y) + (dimensions.y % VOLUME_MASK_STORAGE_BLOCK_SIZE_Y == 0 ? 0 : 1),
            (dimensions.z / VOLUME_MASK_STORAGE_BLOCK_SIZE_Z) + (dimensions.z % VOLUME_MASK_STORAGE_BLOCK_SIZE_Z == 0 ? 0 : 1)
            )
{
    size_t numVoxelsPhyical = tgt::hmul(dimInBlocks_)*VOLUME_MASK_STORAGE_BLOCK_SIZE_X*VOLUME_MASK_STORAGE_BLOCK_SIZE_Y*VOLUME_MASK_STORAGE_BLOCK_SIZE_Z;
    size_t fileSize = numVoxelsPhyical/VOXELS_PER_BYTE + (((numVoxelsPhyical%VOXELS_PER_BYTE) == 0) ? 0 : 1);

    boost::iostreams::mapped_file_params openParams;
    openParams.path = filename_;
    openParams.mode = std::ios::in | std::ios::out;
    openParams.length = fileSize;
    openParams.new_file_size = fileSize; // Option: Do create the file (overwrite if it does not exist)!

    file_.open(openParams);
}
VolumeMaskStorage::~VolumeMaskStorage() {
    file_.close();
    tgt::FileSystem::deleteFile(filename_);
}
void VolumeMaskStorage::set(const tgt::svec3& pos, VolumeMaskValue val) {
    tgtAssert(file_.is_open(), "File is not open");
    tgtAssert(pos.x < dimensions_.x && pos.y < dimensions_.y && pos.z < dimensions_.z, "Invalid pos");

    size_t index;
    uint8_t subindex;
    std::tie(index, subindex) = indexAndSubindexFor(pos);

    char currentByte = file_.data()[index];
    char currentByteExceptNewValue = currentByte & (~(0b11 << (2*subindex)));
    char newByte = currentByteExceptNewValue | (static_cast<uint8_t>(val) << (2*subindex));
    file_.data()[index] = newByte;
}
VolumeMaskValue VolumeMaskStorage::get(const tgt::svec3& pos) const {
    tgtAssert(file_.is_open(), "File is not open");
    tgtAssert(pos.x < dimensions_.x && pos.y < dimensions_.y && pos.z < dimensions_.z, "Invalid pos");

    size_t index;
    uint8_t subindex;
    std::tie(index, subindex) = indexAndSubindexFor(pos);
    return static_cast<VolumeMaskValue>((file_.const_data()[index] >> (2*subindex)) & 0b11);
}


} // namespace voreen
