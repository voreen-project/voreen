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

//header file
#include "cmhaloidtimesteplinker.h"

#include <assert.h>

//we are in namespace voreen
namespace voreen {

CMHaloIDTimeStepLinker::CMHaloIDTimeStepLinker()
    : Processor()
    , inport_(Port::INPORT, "haloport.input", "Halo Data Input")
    , selectedHaloIDProp_("selectedHaloIDProp", "ID of focused halo", 257, -1, 1000000)
    , mouseOverHaloIDProp_("mouseOverHaloIDProp", "ID of hovered over halo", CMMergerTree::NO_HALO_ID, -1, 1000000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , timeStep_("timeStep", "Time Step", 0.0f, 0.0f, 1.0f)
    , updateHaloIDProp_("updateHaloIDProp", "Follow Halo", true)
    , targetSelectionLock_(false)
{
    addPort(inport_);

    addProperty(selectedHaloIDProp_);
    addProperty(mouseOverHaloIDProp_);
    addProperty(timeStep_);
    addProperty(updateHaloIDProp_);

    inport_.onNewData(MemberFunctionCallback<CMHaloIDTimeStepLinker>(this, &CMHaloIDTimeStepLinker::selectedHaloChanged));
    ON_CHANGE(selectedHaloIDProp_, CMHaloIDTimeStepLinker, selectedHaloChanged);
    ON_CHANGE(timeStep_, CMHaloIDTimeStepLinker, timeStepChanged);
}

CMHaloIDTimeStepLinker::~CMHaloIDTimeStepLinker(){
}

void CMHaloIDTimeStepLinker::initialize() {
    // call superclass function first
    Processor::initialize();
}

void CMHaloIDTimeStepLinker::deinitialize() {
    Processor::deinitialize();
}
void CMHaloIDTimeStepLinker::calculateHaloPath() {
    const CMHalo* selectedHalo;
    try {
        selectedHalo = getSelectedHalo();
    } catch(...) {
        return;
    }
    if(std::find(selectedPath_.cbegin(), selectedPath_.cend(), selectedHalo)!=selectedPath_.cend()) {
        return;
    }
    selectedPath_.clear();
    selectedPath_.push_back(selectedHalo);
    const CMHalo* parent = selectedHalo->parent();
    while(parent) {
        const CMHalo* current = parent;
        while(current) {
            if(current->mass > parent->mass) {
                parent = current;
            }
            current = current->spouse();
        }
        selectedPath_.push_front(parent);
        parent = parent->parent();
    }
    const CMHalo* child = selectedHalo->descendant();
    while(child) {
        selectedPath_.push_back(child);
        child = child->descendant();
    }
}
const CMHalo* CMHaloIDTimeStepLinker::getSelectedHalo() {
    const CMMergerTree* tree = inport_.getData();
    if(!tree) {
        throw "no tree";
    }
    const CMHalo* selectedHalo = tree->haloByID(selectedHaloIDProp_.get());
    if(!selectedHalo) {
        throw "no halo";
    }
    return selectedHalo;
}

void CMHaloIDTimeStepLinker::selectedHaloChanged() {
    const CMHalo* selectedHalo;
    try {
        selectedHalo = getSelectedHalo();
    } catch(...) {
        return;
    }
    calculateHaloPath();

    targetSelectionLock_ = true;
    timeStep_.set(selectedHalo->scale);
    targetSelectionLock_ = false;
}

void CMHaloIDTimeStepLinker::timeStepChanged() {
    if(targetSelectionLock_) {
        return;
    }
    if(!updateHaloIDProp_.get()) {
        return;
    }
    const CMHalo* selectedHalo;
    try {
        selectedHalo = getSelectedHalo();
    } catch(...) {
        return;
    }
    float targetTimeStep = timeStep_.get();
    auto it = selectedPath_.cbegin();
    const CMHalo* bestCenterHalo = *it++;
    for(; it != selectedPath_.cend(); ++it) {
        if(bestCenterHalo->scale + (*it)->scale > 2*targetTimeStep) {
            break;
        }
        bestCenterHalo = *it;
    }
    //Will call selectedHaloChanged and thus prepare do the rest of the work
    selectedHaloIDProp_.set(bestCenterHalo->ID);
}

void CMHaloIDTimeStepLinker::process() {
}

} // namespace
