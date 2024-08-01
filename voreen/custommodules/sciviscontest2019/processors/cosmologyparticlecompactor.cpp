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

#include "cosmologyparticlecompactor.h"

namespace voreen {


CosmologyParticleCompactor::CosmologyParticleCompactor()
    : Processor()
    , outport_(Port::OUTPORT, "particlehandle.output", "Particle Data Output")
    , inport_( Port::INPORT,  "particlehandle.input",  "Particle Data Input")
{
    addPort(inport_);
    addPort(outport_);

}
namespace{

class LazyTimeSliceVector : public CMParticleDataTimeSliceVector{
public:
    LazyTimeSliceVector(tgt::Bounds universeBounds);
    virtual const CMParticle* startUsingParticles();
    virtual int               getNumberOfParticles();
    virtual float             getTimeStep();
    virtual const int*        getRemappingTable();
    virtual tgt::mat4         getNormalisationTransformation() override;


    void initialize();

    tgt::mat4 normMatrix_;
    bool initialized_;
    std::shared_ptr<std::vector<int> > mapping_;
    CMParticleDataTimeSlice* parent_;

};
LazyTimeSliceVector::LazyTimeSliceVector(tgt::Bounds universeBounds)
    : CMParticleDataTimeSliceVector(universeBounds)
{
}

const CMParticle* LazyTimeSliceVector::startUsingParticles(){
    if (!initialized_){
        initialize();
    }
    return CMParticleDataTimeSliceVector::startUsingParticles();
}

tgt::mat4 LazyTimeSliceVector::getNormalisationTransformation(){
    return normMatrix_;
}

int LazyTimeSliceVector::getNumberOfParticles(){
    return mapping_->size();
}

float LazyTimeSliceVector::getTimeStep(){
    return CMParticleDataTimeSliceVector::getTimeStep();
}

const int* LazyTimeSliceVector::getRemappingTable(){
    if (!initialized_){
        initialize();
    }
    return CMParticleDataTimeSliceVector::getRemappingTable();
}

void LazyTimeSliceVector::initialize(){
    const CMParticle* particles = parent_->startUsingParticles();
    int numberOfParticles_ = mapping_->size();
    int* mapping = mapping_->data();
    const int* inverseMapping = parent_->getRemappingTable();

    particles_.resize(numberOfParticles_);

    for(int i = 0; i != numberOfParticles_; i++){
        CMParticle p = particles[inverseMapping[mapping[i]]];
        p.ident = i;
        particles_[i] = p;
    }
    parent_->finishedUsingParticles();
    initialized_ = true;
}
}

void CosmologyParticleCompactor::process(){
    const CMParticleData* oldParticleData = inport_.getData();
    if (!oldParticleData)
        return;

    std::vector<CMParticleDataTimeSlice*> slices = oldParticleData->particleDataTimeSlices();
    std::vector<CMParticleDataTimeSlice*> newSlices;

    int numberOfParticles = slices.at(0)->getNumberOfParticles();


    std::shared_ptr<std::vector<int> > mapping = std::make_shared<std::vector<int>>();
    const char * enabled = oldParticleData->getEnabledState();
    for(int i = 0; i != numberOfParticles; i++){
        if (enabled[i]){
            mapping->push_back(i);
        }
    }

    for(auto slice : slices){
        LazyTimeSliceVector* sl = new LazyTimeSliceVector(slice->getUniverseBounds());
        sl->timeStep_ = slice->getTimeStep();
        sl->h0_ = slice->geth0();
        sl->initialized_ = false;
        sl->parent_ = slice;
        sl->mapping_ = mapping;
        sl->normMatrix_ = slice->getNormalisationTransformation();
        newSlices.push_back(sl);
    }

    CMParticleData *newParticleData = new CMParticleData(newSlices);

    outport_.setData(newParticleData);

}

}   // end of namespace
