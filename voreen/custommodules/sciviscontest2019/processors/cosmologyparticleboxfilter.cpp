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

#include "cosmologyparticleboxfilter.h"
#include "../utils/cmmath.h"
namespace voreen {


CosmologyParticleBoxFilter::CosmologyParticleBoxFilter()
    : Processor()
    , outport_(Port::OUTPORT, "particlehandle.output", "Particle Data Output")
    , inport_( Port::INPORT,  "particlehandle.input",  "Particle Data Input")
    , timeStep_(   "timeStep",    "Time Step", 0.0f, 0.0f, 1.0f)
    , autoEval_(   "autoEval",    "Autoevaluate filter", false)
    , evalButton_( "evalButton",  "Evaluate filter")
    , boundingBox_("boundingBox", "Boundingbox")
{
    addPort(inport_);

    addPort(outport_);

    addProperty(timeStep_);
    addProperty(boundingBox_);
    addProperty(autoEval_);
    addProperty(evalButton_);

    ON_CHANGE(timeStep_,   CosmologyParticleBoxFilter, updateRange);
    ON_CHANGE(timeStep_,   CosmologyParticleBoxFilter, onEvalButtonPressed);
    ON_CHANGE(autoEval_,   CosmologyParticleBoxFilter, onAutoEvalChanged);
    ON_CHANGE(evalButton_, CosmologyParticleBoxFilter, onEvalButtonPressed);

    inport_.onNewData(MemberFunctionCallback<CosmologyParticleBoxFilter>(this, &CosmologyParticleBoxFilter::updateRange));
    inport_.onNewData(MemberFunctionCallback<CosmologyParticleBoxFilter>(this, &CosmologyParticleBoxFilter::onEvalButtonPressed));
}

void CosmologyParticleBoxFilter::process(){
    if (!autoEval_.get() && !doFilterOnNextEvaluation_)
        return;
    const CMParticleData* oldParticleData = inport_.getData();
    if (!oldParticleData)
        return;

    // cloneWithParticles erstellt ein neue CMParticleData Struktur, die weiterhin
    // die alten Timeslices enthält(zumindest für alle praktischen Anwendungsfälle)
    // allerdings einen neunen Array für enabledState haben kann
    CMParticleData*           newParticleData = oldParticleData->cloneWithParticles();


    CMParticleDataTimeSlice * newTimeSlice    = newParticleData->sliceAtTimeStep(timeStep_.get());

    int                       numOfParticles  = newTimeSlice->getNumberOfParticles();
    // getEnabledStateMutable holt einen Pointer auf den Array der ausdrückt 
    // welche Partikel aktiviert sind und welche nicht. In diesen Array können
    // wir einfach von [0..numOfParticles) reinschreiben
    char*                     enabledState    = newParticleData->getEnabledStateMutable();


    tgt::Bounds               box             = boundingBox_.get();
    const CMParticle*         particles       = newTimeSlice->startUsingParticles();

    for(int i = 0; i != numOfParticles; i++){
        CMParticle p = particles[i];
        bool particleInBox = box.containsPoint(p.pos);
        // Im enabledState array muss der .ident der Partikels und NICHT der Index im 
        // Array benutzt werden, da dieser über mehrere Zeitschritte garantiert konstant
        // bleibt, was bei dem Index nicht so ist.
        enabledState[p.ident] = particleInBox;
    }
    // Wichtig für Speicherverwaltung der Partikel beim echten Datensatz!
    // Der Array "particles" kann hiernach nicht mehr benutzt werden
    newTimeSlice->finishedUsingParticles();

    outport_.setData(newParticleData);

    doFilterOnNextEvaluation_ = false;
}

void CosmologyParticleBoxFilter::updateRange(){
    const CMParticleData* particleData = inport_.getData();
    if (!particleData) return;
    CMParticleDataTimeSlice* slice = particleData->sliceAtTimeStep(timeStep_.get());
    tgt::Bounds bounds = slice->getBounds();
    bounds = CMtransformBounds(CMinvert(slice->getNormalisationTransformation()), bounds);
    if (bounds.getLLF() != boundingBox_.getMinValue() || bounds.getURB() != boundingBox_.getMaxValue()){
        boundingBox_.setMinValue(bounds.getLLF());
        boundingBox_.setMaxValue(bounds.getURB());
        boundingBox_.set(bounds);
    }
}

void CosmologyParticleBoxFilter::onEvalButtonPressed(){
    doFilterOnNextEvaluation_ = true;
    invalidate();
}

void CosmologyParticleBoxFilter::onAutoEvalChanged(){
    evalButton_.setVisibleFlag(!autoEval_.get());
}
}   // end of namespace