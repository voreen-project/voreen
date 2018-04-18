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

#include "seedpointgenerator.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

const std::string SeedpointGenerator::loggerCat_("voreen.SeedpointGenerator");

SeedpointGenerator::SeedpointGenerator()
    : VolumeProcessor(),
    inportVolume_(Port::INPORT, "inport.volume"),
    outportSeeds_(Port::OUTPORT, "outport.seeds"),
    addSeed_("addSeed", "Add Seedpoint"),
    clearLast_("clearLast", "Clear Last"),
    clearAll_("clearAll", "Clear All"),
    numSeeds_("numSeeds", "Number of Seedpoints", 0, 0, 1000, Processor::VALID),
    seedPoint_("seedPoint", "Seedpoint", tgt::vec3(0.f), tgt::vec3(0.f), tgt::vec3(3000.f), Processor::VALID)
{
    addPort(inportVolume_);
    addPort(outportSeeds_);

    addSeed_.onClick(MemberFunctionCallback<SeedpointGenerator>(this, &SeedpointGenerator::addSeed));
    clearLast_.onClick(MemberFunctionCallback<SeedpointGenerator>(this, &SeedpointGenerator::clearLast));
    clearAll_.onClick(MemberFunctionCallback<SeedpointGenerator>(this, &SeedpointGenerator::clearAll));
    numSeeds_.setReadOnlyFlag(true);

    addProperty(seedPoint_);
    addProperty(addSeed_);
    addProperty(clearLast_);
    addProperty(clearAll_);
    addProperty(numSeeds_);
}

Processor* SeedpointGenerator::create() const {
    return new SeedpointGenerator();
}

void SeedpointGenerator::process() {
    if (inportVolume_.hasChanged()) {
        tgt::vec3 dimensions = inportVolume_.getData()->getRepresentation<VolumeRAM>()->getDimensions();
        dimensions.x = dimensions.x - 1.f;
        dimensions.y = dimensions.y - 1.f;
        dimensions.z = dimensions.z - 1.f;
        seedPoint_.setMaxValue(dimensions);
    }

    PointListGeometryVec3* result = new PointListGeometryVec3();
    result->setData(seedPoints_);

    if (!seedPoints_.empty()) {
        outportSeeds_.setData(result);
    }
}

void SeedpointGenerator::serialize(Serializer& s) const  {
    VolumeProcessor::serialize(s);

    s.serialize("seedPoints", seedPoints_);
}

void SeedpointGenerator::deserialize(Deserializer& d) {
    VolumeProcessor::deserialize(d);

    try {
        d.deserialize("seedPoints", seedPoints_);
    }
    catch (XmlSerializationNoSuchDataException&) {
        d.removeLastError();
    }
}

void SeedpointGenerator::addSeed() {
    seedPoints_.push_back(seedPoint_.get());
    updateProperties();

}

void SeedpointGenerator::clearLast() {
    if (seedPoints_.size() > 0) {
        seedPoints_.pop_back();
        updateProperties();
    }
}

void SeedpointGenerator::clearAll() {
    if (seedPoints_.size() > 0) {
        seedPoints_.clear();
        updateProperties();
    }
}

void SeedpointGenerator::updateProperties() {
    numSeeds_.set(static_cast<int>(seedPoints_.size()));
    process();
}

} // voreen namespace
