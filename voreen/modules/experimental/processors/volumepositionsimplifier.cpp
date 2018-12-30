/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "volumepositionsimplifier.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

using tgt::vec3;

const std::string VolumePositionSimplifier::loggerCat_("voreen.VolumePositionSimplifier");

VolumePositionSimplifier::VolumePositionSimplifier()
    : VolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input"),
    outport_(Port::OUTPORT, "volumehandle.output", "volumehandle.output", false),
    offsetHandling_("offsetHandling", "Offset"),
    normalizeSpacing_("normalizeSpacing", "Normalize Spacing", false),
    removeTransformation_("removeTransformation", "Remove Transformation", false)
{
    addPort(inport_);
    addPort(outport_);

    offsetHandling_.addOption("none", "None");
    offsetHandling_.addOption("removeOffset", "Remove Offset");
    offsetHandling_.addOption("center", "Center");
    offsetHandling_.select("none");

    addProperty(normalizeSpacing_);
    addProperty(removeTransformation_);
    addProperty(offsetHandling_);
}

VolumePositionSimplifier::~VolumePositionSimplifier() {}

Processor* VolumePositionSimplifier::create() const {
    return new VolumePositionSimplifier();
}

void VolumePositionSimplifier::process() {
    const VolumeBase* inputVolume = inport_.getData();

    if(removeTransformation_.get()) {
        inputVolume = new VolumeDecoratorReplaceTransformation(inputVolume, tgt::mat4::identity);
    }

    if(normalizeSpacing_.get()) {
        vec3 sp = inputVolume->getSpacing();
        vec3 cubeSize = sp * vec3(inputVolume->getDimensions());

        float scale = 2.0f / max(cubeSize);
        sp *= scale;

        inputVolume = new VolumeDecoratorReplaceSpacing(inputVolume, sp);
    }

    if(offsetHandling_.get() == "removeOffset") {
        inputVolume = new VolumeDecoratorReplaceOffset(inputVolume, tgt::vec3::zero);
    }
    else if(offsetHandling_.get() == "center") {
        vec3 sp = inputVolume->getSpacing();

        //set origin to center volume:
        vec3 cubeSize = sp * vec3(inputVolume->getDimensions());
        inputVolume = new VolumeDecoratorReplaceOffset(inputVolume, -cubeSize/2.0f);
    }

    if(inputVolume == inport_.getData())
        outport_.setData(inputVolume, false);
    else
        outport_.setData(inputVolume);
}

}   // namespace
