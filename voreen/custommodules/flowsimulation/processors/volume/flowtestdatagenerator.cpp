/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can regetRribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is getRributed in the hope that it will be useful, but WITHOUT ANY       *
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

#include <include/voreen/core/datastructures/volume/volumeatomic.h>
#include "flowtestdatagenerator.h"

#include "voreen/core/datastructures/volume/volumeram.h"

namespace voreen {

const std::string FlowTestDataGenerator::loggerCat_("voreen.flowsimulation.FlowTestDataGenerator");

FlowTestDataGenerator::FlowTestDataGenerator()
    : VolumeProcessor()
    , outport_(Port::OUTPORT, "flowtestdatagenerator.outport", "Volume Output")
    , type_("type", "Type")
{
    addPort(outport_);
    addProperty(type_);
    type_.addOption("zero", "zero");
    type_.addOption("constant", "constant");
    type_.addOption("poiseuille", "poiseuille");
    type_.addOption("rotating", "rotating");
    type_.addOption("helical", "helical");
    type_.addOption("coaxial", "coaxial");
    //type_.addOption("chaotic", "chaotic");
    //type_.addOption("wavy", "wavy");
    type_.addOption("stenotic", "stenotic");
}

FlowTestDataGenerator::~FlowTestDataGenerator() {
}

Processor* FlowTestDataGenerator::create() const {
    return new FlowTestDataGenerator();
}

float getR(size_t diameter, const tgt::svec3& pos) {
    float r = diameter * 0.5f;
    float dx = pos.x - r;
    float dy = pos.y - r;
    return std::sqrt(dx*dx + dy*dy) / r;
}

float getZ(size_t length, const tgt::vec3& pos) {
    return static_cast<float>(pos.z) / length;
}

tgt::vec3 zeroFlow(size_t diameter, size_t length, const tgt::svec3& pos, float t) {
    return tgt::vec3::zero;
}

tgt::vec3 constantFlow(size_t diameter, size_t length, const tgt::svec3& pos, float t) {
    if(getR(diameter, pos) <= 1.0)
        return tgt::vec3(0.0f, 0.0f, 1.0f);
    else
        return tgt::vec3::zero;
}

tgt::vec3 poiseuilleFlow(size_t diameter, size_t length, const tgt::svec3& pos, float t) {
    float d = getR(diameter, pos);

    if(d <= 1.0f)
        return tgt::vec3(0.0f, 0.0f, (1.0f - d) * (1.0f - d));
    else
        return tgt::vec3::zero;
}

tgt::vec3 rotatingFlow(size_t diameter, size_t length, const tgt::svec3& pos, float t) {
    float r = diameter * 0.5f;
    float d = getR(diameter, pos);

    float x = (pos.x - r) / (r * std::sqrt(2.0f));
    float y = (pos.y - r) / (r * std::sqrt(2.0f));

    if(d <= 1.0f)
        return tgt::vec3(d*y, -d*x, 0.0f);
    else
        return tgt::vec3::zero;
}

tgt::vec3 helicalFlow(size_t diameter, size_t length, const tgt::svec3& pos, float t) {
    float r = diameter * 0.5f;
    float d = getR(diameter, pos);

    float sqr3 = std::sqrt(3.0f);
    float x = (pos.x - r) / (r * sqr3);
    float y = (pos.y - r) / (r * sqr3);
    float z = 1.0f / sqr3;

    if(d <= 1.0f)
        return tgt::vec3(d*y, -d*x, d*z);
    else
        return tgt::vec3::zero;
}

tgt::vec3 coaxialFlow(size_t diameter, size_t length, const tgt::svec3& pos, float t) {

    const float split = 0.7f;
    const float zone = 0.05f;

    float r = diameter * 0.5f;
    float d = getR(diameter, pos);

    float x = (pos.x - r) / r;
    float y = (pos.y - r) / r;

    if(d <= 1.0f) {
        if(d <= split - zone) {
            return tgt::vec3(0.0f, 0.0f,1.0f);
        }
        else if(d >= split + zone) {
            return tgt::vec3(0.0f, 0.0f, -1.0f);
        }
        else {
            float e = (d - (split-zone) / (2.0f * zone));
            float n = (std::cos(e * 2.0f * tgt::PIf) + 1.0f) / 2.0f;
            return tgt::normalize(tgt::vec3(x, y, n));
        }
    }
    else
        return tgt::vec3::zero;
}

tgt::vec3 stenoticFlow(size_t diameter, size_t length, const tgt::svec3& pos, float t) {

    float r = diameter * 0.5f;
    float x = (pos.x - r) / r;
    float y = (pos.y - r) / r;
    float z = getZ(length, pos);

    float n = -std::sin(z * 2.0f * tgt::PIf);
    float e = (std::cos(z * 2.0f * tgt::PIf) + 2.0f) / 3.0f;

    float d = getR(diameter, pos);
    if(d <= e)
        return tgt::normalize(tgt::vec3(x*n, y*n, 1.0f))*(1.0f - e);
    else
        return tgt::vec3::zero;
}

void FlowTestDataGenerator::process() {

    const size_t diameter = 50;
    const size_t length = 100;

    tgt::svec3 dimensions(diameter, diameter, length);
    VolumeRAM_3xFloat* output = new VolumeRAM_3xFloat(dimensions);

    auto getVelocityProfile = [&] {
        if(type_.get() == "constant") {
            return &constantFlow;
        }
        else if(type_.get() == "poiseuille") {
            return &poiseuilleFlow;
        }
        else if(type_.get() == "rotating") {
            return &rotatingFlow;
        }
        else if(type_.get() == "helical") {
            return &helicalFlow;
        }
        else if(type_.get() == "coaxial") {
            return &coaxialFlow;
        }
        else if(type_.get() == "stenotic") {
            return &stenoticFlow;
        }
        else /*if(type_.get() == "zero")*/ { // Fallback
            return &zeroFlow;
        }
    };

    auto vel = getVelocityProfile();

    tgt::svec3 pos = tgt::svec3::zero;
    for(pos.z=0; pos.z<dimensions.z; pos.z++) {
        for(pos.y=0; pos.y<dimensions.y; pos.y++) {
            for(pos.x=0; pos.x<dimensions.x; pos.x++) {
                output->voxel(pos) = vel(diameter, length, pos, 0.0f);
            }
        }
    }

    VolumeBase* outputVolume = new Volume(output, tgt::vec3::one, tgt::vec3::zero);
    outport_.setData(outputVolume);
}

} // namespace voreen
