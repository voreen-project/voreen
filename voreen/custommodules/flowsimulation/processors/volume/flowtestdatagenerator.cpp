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

#include <include/voreen/core/datastructures/volume/volumeatomic.h>
#include "flowtestdatagenerator.h"

#include "voreen/core/datastructures/volume/volumeram.h"

namespace voreen {

const std::string FlowTestDataGenerator::loggerCat_("voreen.flowsimulation.FlowTestDataGenerator");

FlowTestDataGenerator::FlowTestDataGenerator()
    : VolumeProcessor()
    , outport_(Port::OUTPORT, "flowtestdatagenerator.outport", "Volume Output")
    , type_("type", "Type")
    , resolutionLength_("resolutionLength", "Resolution Length", 50, 10, 500)
    , resolutionDiameter_("resolutionDiameter", "Resolution Diameter", 50, 10, 500)
{
    addPort(outport_);
    addProperty(type_);
    type_.addOption("zero", "zero");
    type_.addOption("constant", "constant");
    type_.addOption("poiseuille", "poiseuille");
    type_.addOption("rotating", "rotating");
    type_.addOption("helical", "helical");
    type_.addOption("stenotic", "stenotic");
    type_.addOption("soteloVelocity", "Velocity by Sotelo et al.");
    type_.addOption("soteloVorticity", "Vorticity by Sotelo et al.");
    addProperty(resolutionLength_);
    addProperty(resolutionDiameter_);
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
        return tgt::vec3(0.0f, 0.0f, 1.0f - d*d);
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

namespace sotelo {

const float e = 2.72f;
const float c = std::sqrt(5 / std::pow(e, 6));
const float k = 0.01f;
const float gamma = 0.5f * tgt::PIf;
const float nu = 0.004;
const float dP = 1.0f; // TODO:

const tgt::vec3 X(1, 0, 0);
const tgt::vec3 Y(0, 1, 0);
const tgt::vec3 Z(0, 0, 1);

tgt::vec3 velocity(size_t diameter, size_t L, const tgt::svec3& pos, float t) {

    float r = getR(diameter, pos);
    float R = diameter * 0.5f;
    float x = pos.x - R;
    float y = pos.y - R;
    size_t z = pos.z;

    tgt::vec3 v = tgt::vec3::zero;

    // Poisueille term.
    v += dP / (4 * nu * L) * (R * R - r * r) * Z;

    // Lamb-Oseen term.
    v += gamma / (2 * tgt::PIf * r * r) * (1.0f - std::exp(-r * r / c * c)) * k * std::sin(2 * tgt::PIf / L * z) * (-X * y + Y * x);

    return v;
}

tgt::vec3 vorticity(size_t diameter, size_t L, const tgt::svec3& pos, float t) {

    float r = getR(diameter, pos);
    float R = diameter * 0.5f;
    float x = pos.x - R;
    float y = pos.y - R;
    size_t z = pos.z;

    tgt::vec3 w = tgt::vec3::zero;

    w.x = x * k * gamma / (L * r * r) * (1 - std::exp(-r * r / c * c)) * std::cos(2 * tgt::PIf / L * z) + dP / (2 * nu * L) * y;
    w.y = y * k * gamma / (L * r * r) * (1 - std::exp(-r * r / c * c)) * std::cos(2 * tgt::PIf / L * z) - dP / (2 * nu * L) * x;
    w.z = -k * gamma / (tgt::PIf*c*c) * std::exp(-r * r / c * c) * std::sin(2 * tgt::PIf / L * z);

    return w;
}

}

void FlowTestDataGenerator::process() {

    const size_t diameter = resolutionDiameter_.get();
    const size_t length = resolutionLength_.get();

    const tgt::vec3 spacing(2.0f / diameter, 2.0f / diameter, 1.0f / length);
    const tgt::vec3 offset(-spacing.x * diameter * 0.5f, -spacing.y * diameter * 0.5f, 0.0f);

    tgt::svec3 dimensions(diameter, diameter, length);
    VolumeRAM_3xFloat* output = new VolumeRAM_3xFloat(dimensions);

    auto getVectorField = [&] {
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
        else if(type_.get() == "stenotic") {
            return &stenoticFlow;
        }
        else if(type_.get() == "soteloVelocity") {
            return &sotelo::velocity;
        }
        else if(type_.get() == "soteloVorticity") {
            return &sotelo::vorticity;
        }
        else /*if(type_.get() == "zero")*/ { // Fallback
            return &zeroFlow;
        }
    };

    auto v = getVectorField();

    tgt::svec3 pos = tgt::svec3::zero;
    for(pos.z=0; pos.z<dimensions.z; pos.z++) {
        for(pos.y=0; pos.y<dimensions.y; pos.y++) {
            for(pos.x=0; pos.x<dimensions.x; pos.x++) {
                output->voxel(pos) = v(diameter, length, pos, 0.0f);
            }
        }
    }

    VolumeBase* outputVolume = new Volume(output, spacing, offset);
    outport_.setData(outputVolume);
}

} // namespace voreen
