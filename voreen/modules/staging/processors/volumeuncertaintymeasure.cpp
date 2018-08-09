/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "volumeuncertaintymeasure.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatoruncertaintymeasure.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

const std::string VolumeUncertaintyMeasure::loggerCat_("voreen.staging.VolumeUncertaintyMeasure");

VolumeUncertaintyMeasure::VolumeUncertaintyMeasure()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , forceUpdate_(true)
{
    inport_.addCondition(new PortConditionVolumeTypeFloat());
    inport_.addCondition(new PortConditionVolumeValueRange<float>(0.0f, 1.0f));

    addPort(inport_);
    addPort(outport_);

    enableProcessing_.onChange(MemberFunctionCallback<VolumeUncertaintyMeasure>(this, &VolumeUncertaintyMeasure::forceUpdate));
    addProperty(enableProcessing_);
    inport_.onNewData(MemberFunctionCallback<VolumeUncertaintyMeasure>(this, &VolumeUncertaintyMeasure::measureVolumeUncertainty));
}

VolumeUncertaintyMeasure::~VolumeUncertaintyMeasure() {
}

Processor* VolumeUncertaintyMeasure::create() const {
    return new VolumeUncertaintyMeasure();
}

void VolumeUncertaintyMeasure::process() {
    if (!enableProcessing_.get())
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
    else if (forceUpdate_)
        measureVolumeUncertainty();
}

// private methods
//
void VolumeUncertaintyMeasure::forceUpdate() {
    forceUpdate_ = true;
}

void VolumeUncertaintyMeasure::measureVolumeUncertainty() {
    if (!inport_.isReady())
        return;

    const VolumeBase* handle = inport_.getData();
    if (!handle)
        return;

    forceUpdate_ = false;

    if (handle->getRepresentation<VolumeRAM>()) {
        try {
            Volume* v = VolumeOperatorUncertaintyMeasure::APPLY_OP(handle, this);
            outport_.setData(v);
        }
        catch (const std::bad_alloc&) {
            LERROR("measureVolumeUncertainty(): bad allocation");
            outport_.setData(0);
        }
    }
    else {
        outport_.setData(0);
    }
}

}   // namespace
