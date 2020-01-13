/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "timestepduration.h"

#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

const std::string TimestepDuration::loggerCat_("voreen.celltracking.TimestepDuration");

TimestepDuration::TimestepDuration()
    : Processor()
    , inport_(Port::INPORT, "volumecollection", "VolumeList Input", false)
    , outport_(Port::OUTPORT, "textport", "Elapsed Time")
    , startTimestep_("starttimestep", "First timestep (t = 0)", 0, 0, 1, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC)
    , timestepDuration_("stepduration", "Timestep duration (seconds)", 1.f, 0.01f, 1000.f)
    , currentTimestep_("currenttimestep", "Current timestep", 0, 0, 1, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC)
    , elapsedTime_("elapsedtime", "Elapsed time (seconds)", 0.f, -FLT_MAX, FLT_MAX, Processor::VALID)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(startTimestep_);
    addProperty(timestepDuration_);
    addProperty(currentTimestep_);
    elapsedTime_.setReadOnlyFlag(true);
    addProperty(elapsedTime_);
}

Processor* TimestepDuration::create() const {
    return new TimestepDuration();
}

void TimestepDuration::process() {

    outport_.clear();

    if (!inport_.getData() || inport_.getData()->empty()) {
        LWARNING("Empty volume list");
        return;
    }
 
    // compute current time step
    int steps = currentTimestep_.get() - startTimestep_.get();
    float currentTime = timestepDuration_.get() * static_cast<float>(steps);

    elapsedTime_.set(currentTime);
    std::stringstream s;
    s << currentTime << " s";
    outport_.setData(s.str());
}

void TimestepDuration::adjustPropertiesToInput() {

    const VolumeList* collection = inport_.getData();
    int max = ((collection != 0) ? static_cast<int>(collection->size()) : 0);

    if (collection && !collection->empty()) {
        startTimestep_.setMaxValue(max - 1);
        currentTimestep_.setMaxValue(max - 1);
    }

}

} // namespace
