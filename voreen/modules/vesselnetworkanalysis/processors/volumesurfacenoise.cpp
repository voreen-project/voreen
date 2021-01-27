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

#include "volumesurfacenoise.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include <unordered_map>
#include <vector>
#include <random>

namespace voreen {

const std::string VolumeSurfaceNoise::loggerCat_("voreen.vesselnetworkanalysis.volumesurfacenoise");

VolumeSurfaceNoise::VolumeSurfaceNoise()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output", false, Processor::VALID /* do not update if something connects to this port*/)
    , enabledProp_("enabledProp","Enabled",true)
    , noiseAmount_("noiseAmount", "Noise Amount", 0, 0, 1)
    , binarizationThreshold_("binarizationThreshold", "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , usePredeterminedSeed_("usePredeterminedSeed", "Use Predetermined Seed", false)
    , predeterminedSeed_("predeterminedSeed", "Seed", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , preserveTopology_("preserveTopology", "Preserve Topology", true)
    , randomDevice_()
{
    addPort(inport_);
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(noiseAmount_);
    addProperty(binarizationThreshold_);
    addProperty(usePredeterminedSeed_);
    addProperty(predeterminedSeed_);
    addProperty(preserveTopology_);
}

VolumeSurfaceNoise::~VolumeSurfaceNoise() {
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

struct VoxelSet {
    std::vector<tgt::ivec3> voxels;
    std::unordered_map<tgt::ivec3, size_t> voxelIndices;

    VoxelSet()
        : voxels()
        , voxelIndices()
    {
    }

    void insert(tgt::ivec3 p) {
        if(voxelIndices.count(p) == 0) {
            voxels.push_back(p);
            voxelIndices.insert({p, voxels.size()-1});
        }
    }
    void remove(tgt::ivec3 p) {
        auto it = voxelIndices.find(p);
        if(it != voxelIndices.end()) {
            size_t index = it->second;
            if(index != voxels.size()-1) {
                std::swap(voxels[index], voxels.back());
            }
            voxels.pop_back();
            voxelIndices.erase(it);
        }
    }

    tgt::ivec3 drawRandom(std::mt19937& randomEngine) {
        size_t random_index = std::uniform_int_distribution<size_t>(0, voxels.size()-1)(randomEngine);
        tgtAssert(random_index < voxels.size(), "Invalid index");
        return voxels[random_index];
    }

    size_t size() const {
        tgtAssert(voxels.size() == voxelIndices.size(), "Size mismatch");
        return voxels.size();
    }
};

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

void VolumeSurfaceNoise::process() {
    if(!enabledProp_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    if(!inport_.getData()) {
        return;
    }
    const VolumeBase& inVol = *inport_.getData();
    const tgt::ivec3 dim = inVol.getDimensions();

    const VolumeRAM* inputPtr = inVol.getRepresentation<VolumeRAM>();
    const VolumeRAM& input = *inputPtr;

    float binarizationThresholdNormalized;
    if(inVol.hasMetaData(VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING)) {
        // If the input volume does not have a RealWorldMapping we need to convert the binarizationThreshold to a normalized value.
        binarizationThresholdNormalized = inVol.getRealWorldMapping().realWorldToNormalized(binarizationThreshold_.get());
    } else {
        // If the input volume does not have a RealWorldMapping we expect RW values to be normalized.
        binarizationThresholdNormalized = binarizationThreshold_.get();
    }

    voreen::VolumeRAM* outputPtr = VolumeFactory().create("uint8", dim);
    voreen::VolumeRAM& output = *outputPtr;

    for(int z = 0; z < dim.z; ++z) {
        for(int y = 0; y < dim.y; ++y) {
            for(int x = 0; x < dim.x; ++x) {
                if(input.getVoxelNormalized(x, y, z) > binarizationThresholdNormalized) {
                    outputPtr->setVoxelNormalized(1.0, x, y, z);
                } else {
                    outputPtr->setVoxelNormalized(0.0, x, y, z);
                }
            }
        }
    }

    const bool FOREGROUND = true;
    const bool BACKGROUND = false;

    auto setVoxel = [&] (tgt::ivec3 p, bool newVal) {
        if(newVal == FOREGROUND){
            output.setVoxelNormalized(1.0, p.x, p.y, p.z);
        } else {
            output.setVoxelNormalized(0.0, p.x, p.y, p.z);
        }
    };
    auto getVoxel = [&] (tgt::ivec3 p, bool outsideVolumeValue) {
        if (0 > p.x || p.x >= dim.x ||
            0 > p.y || p.y >= dim.y ||
            0 > p.z || p.z >= dim.z) {
            return outsideVolumeValue;
        }
        if(output.getVoxelNormalized(p.x, p.y, p.z) > binarizationThresholdNormalized){
            return FOREGROUND;
        } else {
            return BACKGROUND;
        }
    };

    auto isSurfaceVoxel = [&] (tgt::ivec3 p) {
        if (0 > p.x || p.x >= dim.x ||
            0 > p.y || p.y >= dim.y ||
            0 > p.z || p.z >= dim.z) {
            return false;
        }
        if(getVoxel(p, FOREGROUND) == FOREGROUND) {
            return getVoxel(tgt::ivec3(p.x+1, p.y  , p.z  ), FOREGROUND) == BACKGROUND
                || getVoxel(tgt::ivec3(p.x-1, p.y  , p.z  ), FOREGROUND) == BACKGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y+1, p.z  ), FOREGROUND) == BACKGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y-1, p.z  ), FOREGROUND) == BACKGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y  , p.z+1), FOREGROUND) == BACKGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y  , p.z-1), FOREGROUND) == BACKGROUND;
        } else {
            return getVoxel(tgt::ivec3(p.x+1, p.y  , p.z  ), BACKGROUND) == FOREGROUND
                || getVoxel(tgt::ivec3(p.x-1, p.y  , p.z  ), BACKGROUND) == FOREGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y+1, p.z  ), BACKGROUND) == FOREGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y-1, p.z  ), BACKGROUND) == FOREGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y  , p.z+1), BACKGROUND) == FOREGROUND
                || getVoxel(tgt::ivec3(p.x  , p.y  , p.z-1), BACKGROUND) == FOREGROUND;
        }
    };

    auto isEulerInvariantVoxel = [&] (tgt::ivec3 pos) {
        bool neighborhood[26];
        int i = 0;
        for(int dz = -1; dz<=1; ++dz) {
            for(int dy = -1; dy<=1; ++dy) {
                for(int dx = -1; dx<=1; ++dx) {
                    if(dx != 0 || dy != 0 ||  dz!= 0) {
                        neighborhood[i++] = (getVoxel(tgt::ivec3(pos.x + dx, pos.y + dy, pos.z + dz), FOREGROUND) == FOREGROUND);
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
    };

    auto isSimple = [&] (const tgt::svec3& pos) {
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
                    if((dx != 0 || dy != 0 ||  dz!= 0)) {
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
    };

    VoxelSet surface;

    auto tryAddToSurface = [&] (tgt::ivec3 p) {
        if(isSurfaceVoxel(p)) {
            surface.insert(p);
        }
    };

    for(int z = 0; z < dim.z; ++z) {
        for(int y = 0; y < dim.y; ++y) {
            for(int x = 0; x < dim.x; ++x) {
                tryAddToSurface(tgt::ivec3(x, y, z));
            }
        }
    }

    std::mt19937 randomEngine(randomDevice_());
    if(usePredeterminedSeed_.get()) {
        randomEngine.seed(predeterminedSeed_.get());
    }

    const float noiseAmount = noiseAmount_.get();
    const size_t numSurfaceVoxelsToAdd = std::min(surface.size(), static_cast<size_t>(surface.size() * noiseAmount));

    bool preserveTopology = preserveTopology_.get();
    for(size_t i=0; i<numSurfaceVoxelsToAdd && surface.size() > 0;) {
        auto p = surface.drawRandom(randomEngine);

        surface.remove(p);
        if(preserveTopology && (!isSimple(p) || !isEulerInvariantVoxel(p))) {
            continue;
        }
        if(getVoxel(p, FOREGROUND) == FOREGROUND) {
            setVoxel(p, BACKGROUND);
        } else {
            setVoxel(p, FOREGROUND);
        }

        tryAddToSurface(tgt::ivec3(p.x  , p.y  , p.z  ));
        tryAddToSurface(tgt::ivec3(p.x+1, p.y  , p.z  ));
        tryAddToSurface(tgt::ivec3(p.x-1, p.y  , p.z  ));
        tryAddToSurface(tgt::ivec3(p.x  , p.y+1, p.z  ));
        tryAddToSurface(tgt::ivec3(p.x  , p.y-1, p.z  ));
        tryAddToSurface(tgt::ivec3(p.x  , p.y  , p.z+1));
        tryAddToSurface(tgt::ivec3(p.x  , p.y  , p.z-1));

        ++i;
    }

    Volume* outvol = new Volume(outputPtr, &inVol);

    outport_.setData(outvol);
}
VoreenSerializableObject* VolumeSurfaceNoise::create() const {
    return new VolumeSurfaceNoise();
}

} // namespace voreen
