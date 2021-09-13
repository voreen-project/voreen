/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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
#include "cosmologyfakedatasource.h"

#include <random>

namespace voreen{

CosmologyFakeDataSource::CosmologyFakeDataSource()
    : Processor()
    , outport_(Port::OUTPORT, "particlehandle.output", "Particle Data Output")
    , particleCount_("particleCount", "Particle Count", 10000, 1, 1000000)
{
    addPort(outport_);
    addProperty(particleCount_);
}

CosmologyFakeDataSource::~CosmologyFakeDataSource(){
}

Processor* CosmologyFakeDataSource::create() const{
    return new CosmologyFakeDataSource;
}

void CosmologyFakeDataSource::process(){
    std::vector<CMParticleDataTimeSlice*> timeSlices;
    
    timeSlices.push_back(generateTimeSlice(0.0f));
    CMParticleData *particleData = new CMParticleData(timeSlices);

    outport_.setData(particleData);
}

CMParticleDataTimeSlice* CosmologyFakeDataSource::generateTimeSlice(float time) {
    CMParticleDataTimeSliceVector *timeSlice = new CMParticleDataTimeSliceVector(tgt::Bounds(tgt::vec3(-2000), tgt::vec3(2000)));
    timeSlice->timeStep_ = 1.0f;
    timeSlice->h0_ = 1.0f;

    int numberOfParticles = particleCount_.get();
    
    std::minstd_rand                      generator(5);
    std::uniform_real_distribution<float> posDistribution(-2000, 2000);
    std::uniform_real_distribution<float> velDistribution(-3, 3);

    for(size_t i = 0; i < numberOfParticles; ++i) {
        CMParticle p;
        p.pos.x = posDistribution(generator);
        p.pos.y = posDistribution(generator);
        p.pos.z = posDistribution(generator);
        p.pos   = posDistribution(generator)*tgt::normalize(p.pos);
        p.vel.x = velDistribution(generator);
        p.vel.y = velDistribution(generator);
        p.vel.z = velDistribution(generator);
        p.phi   = 0.0f;
        p.ident = i;

        timeSlice->particles_.push_back(p);
    }
    return timeSlice;
}
}
