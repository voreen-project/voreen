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

#include "volumesurfacenoise.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include <boost/functional/hash.hpp>

#include <unordered_map>
#include <vector>
#include <random>

//Hash implementation for tgt::ivec3
namespace std {
    template <> struct hash<tgt::ivec3>
    {
        size_t operator()(const tgt::ivec3 & p) const
        {
            std::array<int, 3> array{ {p.x, p.y, p.z} };
            return boost::hash_value(array);
        }
    };
}

namespace voreen {

const std::string VolumeSurfaceNoise::loggerCat_("voreen.vesselnetworkanalysisextra.volumesurfacenoise");

VolumeSurfaceNoise::VolumeSurfaceNoise()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output", false, Processor::VALID /* do not update if something connects to this port*/)
    , enabledProp_("enabledProp","Enabled",true)
    , noiseAmount_("noiseAmount", "Noise Amount", 0, 0, 1)
    , binarizationThreshold_("binarizationThreshold", "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , usePredeterminedSeed_("usePredeterminedSeed", "Use Predetermined Seed", false)
    , predeterminedSeed_("predeterminedSeed", "Seed", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , randomDevice_()
{
    addPort(inport_);
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(noiseAmount_);
    addProperty(binarizationThreshold_);
    addProperty(usePredeterminedSeed_);
    addProperty(predeterminedSeed_);
}

VolumeSurfaceNoise::~VolumeSurfaceNoise() {
}

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
    if(inVol.hasMetaData("RealWorldMapping")) {
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

    for(size_t i=0; i<numSurfaceVoxelsToAdd && surface.size() > 0; ++i) {
        auto p = surface.drawRandom(randomEngine);

        if(getVoxel(p, FOREGROUND) == FOREGROUND) {
            setVoxel(p, BACKGROUND);
        } else {
            setVoxel(p, FOREGROUND);
        }

        surface.remove(p);

        tryAddToSurface(tgt::ivec3(p.x  , p.y  , p.z  ));
        tryAddToSurface(tgt::ivec3(p.x+1, p.y  , p.z  ));
        tryAddToSurface(tgt::ivec3(p.x-1, p.y  , p.z  ));
        tryAddToSurface(tgt::ivec3(p.x  , p.y+1, p.z  ));
        tryAddToSurface(tgt::ivec3(p.x  , p.y-1, p.z  ));
        tryAddToSurface(tgt::ivec3(p.x  , p.y  , p.z+1));
        tryAddToSurface(tgt::ivec3(p.x  , p.y  , p.z-1));
    }

    Volume* outvol = new Volume(outputPtr, &inVol);

    outport_.setData(outvol);
}
VoreenSerializableObject* VolumeSurfaceNoise::create() const {
    return new VolumeSurfaceNoise();
}

} // namespace voreen
