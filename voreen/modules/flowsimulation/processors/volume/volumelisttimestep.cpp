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

#include "volumelisttimestep.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeListTimeStep::loggerCat_("voreen.base.VolumeListTimeStep");

VolumeListTimeStep::VolumeListTimeStep()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumelist.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumelist.output", "Volume Output", false)
    , enableProcessing_("enabled", "Enable", false)
    , timeStep_("timeStep", "Time Step (s)", 1.0f, 0.0f, 10000.0f)
    , startTime_("startTime", "Start Time (s)", 0.0f, -100000.0f, 100000.0f)
{
    addPort(inport_);
    inport_.Observable<PortObserver>::addObserver(this);
    addPort(outport_);
    addProperty(enableProcessing_);

    addProperty(timeStep_);
    timeStep_.setNumDecimals(4);
    addProperty(startTime_);
    startTime_.setNumDecimals(4);
}

VolumeListTimeStep::~VolumeListTimeStep() {
    // We need to remove the observer, since outport is destructed before inport (stack hierarchy)
    // and the destruction of inport will clear the (already destructed) outport otherwise.
    inport_.Observable<PortObserver>::removeObserver(this);
}

Processor* VolumeListTimeStep::create() const {
    return new VolumeListTimeStep();
}

void VolumeListTimeStep::initialize() {
    VolumeProcessor::initialize();
}

void VolumeListTimeStep::deinitialize() {
    clearOutput();
    VolumeProcessor::deinitialize();
}

void VolumeListTimeStep::process() {
    clearOutput();

    const VolumeList* inputList = inport_.getData();

    if (!enableProcessing_.get() || !inputList || inputList->empty()) {
        outport_.setData(inputList, false);
        return;
    }

    VolumeList* outputList = new VolumeList();

    for (size_t i = 0; i < inputList->size(); ++i) {
        const VolumeBase* inputVolume = inputList->at(i);

        float timeStep = startTime_.get() + i * timeStep_.get();

        VolumeBase* outputVolume = new VolumeDecoratorReplaceTimestep(inputVolume, timeStep);
        outputList->add(outputVolume);

        decorators_.emplace_back(std::unique_ptr<VolumeBase>(outputVolume));
    }
    outport_.setData(outputList, true);
}

void VolumeListTimeStep::clearOutput() {
    outport_.clear();
    decorators_.clear();
}

void VolumeListTimeStep::dataWillChange(const Port* source) {
    clearOutput();
}

}   // namespace
