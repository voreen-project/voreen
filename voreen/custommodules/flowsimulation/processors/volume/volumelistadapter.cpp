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

#include "volumelistadapter.h"

#include "voreen/core/datastructures/volume/volumedecorator.h"

namespace voreen {

VolumeListAdapter::VolumeListAdapter()
    : Processor()
    , inport_(Port::INPORT, "volume.input", "Volume Input", true)
    , outport_(Port::OUTPORT, "volumelist.output", "Volume Output", false)
{
    addPort(inport_);
    addPort(outport_);
}

VolumeListAdapter::~VolumeListAdapter() {}

Processor* VolumeListAdapter::create() const {
    return new VolumeListAdapter();
}

void VolumeListAdapter::process() {

    // Delete old data.
    decorators_.clear();

    if(!inport_.hasData()) {
        outport_.setData(nullptr);
    }
    else if(inport_.hasChanged()) {
        VolumeList* output = new VolumeList();

        for(const VolumeBase* volume : inport_.getAllData()) {
            VolumeBase* identity = new VolumeDecoratorIdentity(volume);
            output->add(identity);
            decorators_.push_back(std::unique_ptr<VolumeBase>(identity));
        }

        outport_.setData(output, true);
    }
}

}   // namespace
