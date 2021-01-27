/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "volumecollector.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/processors/processorwidgetfactory.h"

namespace voreen {

const std::string VolumeCollector::loggerCat_("voreen.VolumeCollector");

VolumeCollector::VolumeCollector()
    : Processor(),
      inport_(Port::INPORT, "inport", "inport", true),
      outport_(Port::OUTPORT, "outport")
{
    addPort(inport_);
    addPort(outport_);
}

Processor* VolumeCollector::create() const {
    return new VolumeCollector();
}

void VolumeCollector::process() {
    update();
}

void VolumeCollector::update() {
    if(inport_.isReady()) {
        VolumeList* vc = new VolumeList();

        std::vector<const VolumeBase*> data = inport_.getAllData();
        for(size_t i=0; i<data.size(); i++)
            vc->add(const_cast<VolumeBase*>(data[i]));

        outport_.setData(vc, true);
    }
    else
        outport_.setData(0);
}

void VolumeCollector::invalidate(int inv) {
    update();

    Processor::invalidate(inv);
}

} // namespace
