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

#include "cosmologyparticlehalofilter.h"

namespace voreen {


CosmologyParticleHaloFilter::CosmologyParticleHaloFilter()
    : Processor()
    , outport_(Port::OUTPORT, "particlehandle.output", "Particle Data Output")
    , inport_( Port::INPORT,  "particlehandle.input",  "Particle Data Input")
    , haloInport_(Port::INPORT, "haloport.input", "Halo Data Input")
    , autoEval_(   "autoEval",    "Autoevaluate filter", false)
    , evalButton_( "evalButton",  "Evaluate filter")
    , radiusDivProp_("radiusDivProp", "Radius divider" ,10000.0f ,1.0f, 100000.f)
    , selectedHaloIDProp_("selectedHaloIDProp", "ID of focused halo", 257, -1, 1000000)
{
    addPort(inport_);
    addPort(haloInport_);
    addPort(outport_);

    addProperty(autoEval_);
    addProperty(evalButton_);
    addProperty(radiusDivProp_);
    addProperty(selectedHaloIDProp_);

    ON_CHANGE(selectedHaloIDProp_,   CosmologyParticleHaloFilter, onEvalButtonPressed);
    ON_CHANGE(autoEval_,   CosmologyParticleHaloFilter, onAutoEvalChanged);
    ON_CHANGE(evalButton_, CosmologyParticleHaloFilter, onEvalButtonPressed);
    inport_.onNewData(MemberFunctionCallback<CosmologyParticleHaloFilter>(this, &CosmologyParticleHaloFilter::onEvalButtonPressed));
}

void CosmologyParticleHaloFilter::process(){
    if (!autoEval_.get() && !doFilterOnNextEvaluation_)
        return;
    const CMParticleData* oldParticleData = inport_.getData();
    if (!oldParticleData)
        return;

    const CMMergerTree* tree = haloInport_.getData();
    if(!tree) {
        return;
    }

    const CMHalo* halo = tree->haloByID(selectedHaloIDProp_.get());
    float timeStep = halo->scale;
    CMParticleData*           newParticleData = oldParticleData->cloneWithParticles();

    

    CMParticleDataTimeSlice * newTimeSlice    = newParticleData->sliceAtTimeStep(timeStep);

    int                       numOfParticles  = newTimeSlice->getNumberOfParticles();
    char*                     enabledState    = newParticleData->getEnabledStateMutable();


    const CMParticle*         particles       = newTimeSlice->startUsingParticles();

    float radius = halo->radius/radiusDivProp_.get();

    tgt::mat4 norm =newTimeSlice->getNormalisationTransformation();

    for(int i = 0; i != numOfParticles; i++){
        CMParticle p = particles[i];
        tgt::vec4 realPos = norm*tgt::vec4(p.pos, 1.0f);
        bool particleInHalo = distance(tgt::vec3(realPos.x, realPos.y, realPos.z), halo->pos) < radius;        
        enabledState[p.ident] = particleInHalo;
    }
    // Wichtig für Speicherverwaltung der Partikel beim echten Datensatz!
    // Der Array "particles" kann hiernach nicht mehr benutzt werden
    newTimeSlice->finishedUsingParticles();

    outport_.setData(newParticleData);

    doFilterOnNextEvaluation_ = false;
}


void CosmologyParticleHaloFilter::onEvalButtonPressed(){
    doFilterOnNextEvaluation_ = true;
    invalidate();
}

void CosmologyParticleHaloFilter::onAutoEvalChanged(){
    evalButton_.setVisibleFlag(!autoEval_.get());
}
}   // end of namespace
