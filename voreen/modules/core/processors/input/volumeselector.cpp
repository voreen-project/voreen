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

#include "volumeselector.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include <limits.h> //for std::numeric_limits

namespace voreen {

const std::string VolumeSelector::loggerCat_("voreen.core.VolumeSelector");

VolumeSelector::VolumeSelector()
    : Processor(),
      volumeID_("volumeID", "Selected volume", -1, -1, std::numeric_limits<int>::max()-1),
      inport_(Port::INPORT, "volumecollection", "VolumeList Input", false),
      outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output", false)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(volumeID_);
}

Processor* VolumeSelector::create() const {
    return new VolumeSelector();
}

void VolumeSelector::process() {
    // processor is only ready, if inport contains a volumelist
    // but the list can be empty
    if(volumeID_.get() == -1) {
        outport_.setData(nullptr);
    } else {
        outport_.setData(inport_.getData()->at(volumeID_.get()), false);
    }
}

void VolumeSelector::adjustPropertiesToInput() {
    const VolumeList* collection = inport_.getData();

    //if inport is empty, do nothing
    if(!collection) return;

    //if we have a list, adapt min and max values
    volumeID_.setMinValue(std::min(0,static_cast<int>(collection->size())-1));
    volumeID_.setMaxValue(static_cast<int>(collection->size())-1);
    // set to first volume if no volume was present earlier
    if (!collection->empty() && volumeID_.get() == -1)
        volumeID_.set(0);
}

} // namespace
