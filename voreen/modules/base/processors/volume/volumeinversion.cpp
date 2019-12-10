/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "volumeinversion.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorinvert.h"

namespace voreen {

VolumeInversion::VolumeInversion()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , forceUpdate_(true)
{
    addPort(inport_);
    addPort(outport_);

    enableProcessing_.onChange(MemberFunctionCallback<VolumeInversion>(this, &VolumeInversion::forceUpdate));
    addProperty(enableProcessing_);
    inport_.onNewData(MemberFunctionCallback<VolumeInversion>(this, &VolumeInversion::invertVolume));
}

VolumeInversion::~VolumeInversion() {
}

Processor* VolumeInversion::create() const {
    return new VolumeInversion();
}

void VolumeInversion::process() {
    if (!enableProcessing_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
    }
    else if (forceUpdate_) {
        invertVolume();
    }
}

// private methods
//
void VolumeInversion::forceUpdate() {
    forceUpdate_ = true;
}

void VolumeInversion::invertVolume() {
    const VolumeBase* handle = inport_.getData();
    if(!handle)
        return;

    forceUpdate_ = false;

    if (handle->getRepresentation<VolumeRAM>()) {
        //VolumeRAM* v = handle->getRepresentation<VolumeRAM>()->clone();
        Volume* v = VolumeOperatorInvert::APPLY_OP(handle, this);
        outport_.setData(v);
    }
    else {
        outport_.setData(0);
    }
}

}   // namespace
