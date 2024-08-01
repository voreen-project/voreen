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

#include "lymphatictestvesselgenerator.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include <boost/functional/hash.hpp>

#include <unordered_map>
#include <vector>
#include <random>

namespace voreen {

const std::string LymphaticTestVesselGenerator::loggerCat_("voreen.vesselnetworkanalysisextra.lymphatictestvesselgenerator");

LymphaticTestVesselGenerator::LymphaticTestVesselGenerator()
    : VolumeProcessor()
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output", false, Processor::VALID /* do not update if something connects to this port*/)
    , enabledProp_("enabledProp","Enabled",true)
    , dimensions_("dimensions", "Dimensions", tgt::ivec3(5), tgt::ivec3(5), tgt::ivec3(1000))
    , spacing_("spacing", "Spacing", tgt::vec3(1), tgt::vec3(0.01), tgt::vec3(100))
    , usePredeterminedSeed_("usePredeterminedSeed", "Use Predetermined Seed", false)
    , predeterminedSeed_("predeterminedSeed", "Seed", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , randomDevice_()
{
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(dimensions_);
    addProperty(spacing_);
    addProperty(usePredeterminedSeed_);
    addProperty(predeterminedSeed_);
}

LymphaticTestVesselGenerator::~LymphaticTestVesselGenerator() {
}
void LymphaticTestVesselGenerator::process() {
    if(!enabledProp_.get()) {
        outport_.setData(nullptr, false);
        return;
    }
    const tgt::ivec3 dim = dimensions_.get();

    voreen::VolumeRAM* outputPtr = VolumeFactory().create("uint8", dim);
    voreen::VolumeRAM& output = *outputPtr;

    auto isInside = [&] (tgt::vec3 vec) {
        struct F { float zfreq; float zoffset; float roffset; float scale; };
        std::vector<F> fs;
        fs.push_back({0.5, 0.7, 1.0, 0.3});
        fs.push_back({0.5, 0.2, 0.3, 0.3});
        fs.push_back({0.45, 0.6, 0.5, 0.3});
        fs.push_back({1.0, 0.3, 0.4, 0.2});
        fs.push_back({1.5, 0.2, 0.2, 0.1});
        float displacement = 0;
        float angle = std::atan2(vec.x, vec.y);
        for(auto f : fs) {
            float angle_multiplier = (sin(2.0*(angle+(f.roffset*tgt::PIf)))+1)*0.5;
            displacement += angle_multiplier*f.scale*sin(tgt::PIf * 2.0f * (f.zfreq * vec.z + f.zoffset));
        }
        float r = tgt::length(vec.xy())+displacement*0.5;
        return r < 0.5;
    };

    tgt::mat4 toNormalized = tgt::mat4::createTranslation(tgt::vec3(-1.0))*tgt::mat4::createScale(tgt::vec3(2.0)/(tgt::vec3(dim - tgt::ivec3(1))));
    //tgt::mat4 fromNormalized;
    //bool res = toNormalized.invert(fromNormalized);
    //tgtAssert(res, "Inversion failed");

    for(int z=0; z<dim.z; ++z) {
        for(int y=0; y<dim.y; ++y) {
            for(int x=0; x<dim.x; ++x) {
                tgt::vec3 normalized = toNormalized.transform(tgt::vec3(x,y,z));
                float val = isInside(normalized) ? 1.0 : 0.0;
                output.setVoxelNormalized(val, tgt::svec3(x, y, z));
            }
        }
    }

    Volume* outvol = new Volume(outputPtr, spacing_.get(), tgt::vec3(0,0,0) /*offset*/);

    outport_.setData(outvol);
}
VoreenSerializableObject* LymphaticTestVesselGenerator::create() const {
    return new LymphaticTestVesselGenerator();
}

} // namespace voreen
