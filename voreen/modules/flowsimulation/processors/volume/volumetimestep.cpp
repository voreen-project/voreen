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

#include "volumetimestep.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeTimestep::loggerCat_("voreen.base.VolumeTimestep");

VolumeTimestep::VolumeTimestep()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enabled", "Enable", false)
    , timestep_("timestep", "Timestep", 0.0f, -10000000.0f, 10000000.0f)
    , reset_("reset", "Reset Timestep")
{
    addPort(inport_);
    inport_.Observable<PortObserver>::addObserver(this);
    addPort(outport_);

    addProperty(enableProcessing_);

    timestep_.setNumDecimals(6);
    timestep_.setStepping(0.001f);
    addProperty(timestep_);

    addProperty(reset_);
    ON_CHANGE(reset_, VolumeTimestep, resetTimestep);
}

VolumeTimestep::~VolumeTimestep() {
    // We need to remove the observer, since outport is destructed before inport (stack hierarchy)
    // and the destruction of inport will clear the (already destructed) outport otherwise.
    inport_.Observable<PortObserver>::removeObserver(this);
}

Processor* VolumeTimestep::create() const {
    return new VolumeTimestep();
}

void VolumeTimestep::process() {
    if (!enableProcessing_.get()) {
        outport_.setData(inport_.getData(), false);
    }
    else {
        VolumeBase* outputVolume =
                new VolumeDecoratorReplaceTimestep(inport_.getData(), timestep_.get());
        outport_.setData(outputVolume);
    }
}

void VolumeTimestep::adjustPropertiesToInput() {
    const VolumeBase* inputVolume = inport_.getData();

    // only adjust offset, if processor is disabled (= no effect on output volume)
    if (!enableProcessing_.get()) {
        timestep_.set(inputVolume ? inputVolume->getTimestep() : 0.0f);
    }
}

void VolumeTimestep::resetTimestep() {
    const VolumeBase* inputVolume = inport_.getData();
    if (inputVolume) {
        timestep_.set(inputVolume->getTimestep());
    }
}

void VolumeTimestep::dataWillChange(const Port* source) {
    outport_.clear();
}

}   // namespace
