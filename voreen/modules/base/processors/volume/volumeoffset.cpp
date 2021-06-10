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

#include "volumeoffset.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeOffset::loggerCat_("voreen.base.VolumeOffset");

VolumeOffset::VolumeOffset()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enabled", "Enable", false)
    , offset_("offset", "Set Offset", tgt::vec3(0.0f), tgt::vec3(-10000000.0), tgt::vec3(10000000.0))
    , reset_("reset", "Reset Offset")
    , offsetDisplay_("offsetDisplay", "Resulting Offset", tgt::vec3(0.0f), tgt::vec3(-10000000.0), tgt::vec3(10000000.0))
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
        ON_CHANGE(enableProcessing_, VolumeOffset, adjustPropertyVisibility);

    offset_.setNumDecimals(6);
    offset_.setStepping(tgt::vec3(0.001f));
    addProperty(offset_);

    addProperty(reset_);
        ON_CHANGE(reset_, VolumeOffset, resetOffset);

    offsetDisplay_.setReadOnlyFlag(true);
    offsetDisplay_.setNumDecimals(6);
    addProperty(offsetDisplay_);
}

Processor* VolumeOffset::create() const {
    return new VolumeOffset();
}

void VolumeOffset::initialize() {
    VolumeProcessor::initialize();
    adjustPropertyVisibility();
}

void VolumeOffset::process() {
    if (!enableProcessing_.get()) {
        outport_.setData(inport_.getData(), false);
    }
    else {
        VolumeBase* outputVolume =
            new VolumeDecoratorReplaceOffset(inport_.getData(), offset_.get());
        outport_.setData(outputVolume);
    }
    if(outport_.getData()) {
        offsetDisplay_.set(outport_.getData()->getOffset());
    }
}

void VolumeOffset::adjustPropertiesToInput() {
    const VolumeBase* inputVolume = inport_.getData();

    // only adjust offset, if processor is disabled (= no effect on output volume)
    if (!enableProcessing_.get()) {
        offset_.set(inputVolume ? inputVolume->getOffset() : tgt::vec3::zero);
    }
}

void VolumeOffset::resetOffset() {
    const VolumeBase* inputVolume = inport_.getData();
    if (inputVolume) {
        offset_.set(inputVolume->getOffset());
    }
}

void VolumeOffset::adjustPropertyVisibility() {
    bool enabled = enableProcessing_.get();
    offset_.setReadOnlyFlag(!enabled);
    reset_.setReadOnlyFlag(!enabled);
}

}   // namespace
