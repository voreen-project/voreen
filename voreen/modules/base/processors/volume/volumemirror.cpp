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

#include "volumemirror.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormirror.h"

namespace voreen {

const std::string VolumeMirror::loggerCat_("voreen.VolumeMirror");

VolumeMirror::VolumeMirror()
    : CachingVolumeProcessor(),
    mirrorX_("mirrorX", "Mirror X", false),
    mirrorY_("mirrorY", "Mirror Y", false),
    mirrorZ_("mirrorZ", "Mirror Z", false),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
{
    addProperty(mirrorX_);
    addProperty(mirrorY_);
    addProperty(mirrorZ_);

    addPort(inport_);
    addPort(outport_);
}

VolumeMirror::~VolumeMirror() {}

Processor* VolumeMirror::create() const {
    return new VolumeMirror();
}

void VolumeMirror::process() {
    const VolumeBase* inputVolume = inport_.getData();
    Volume* outputVolume = inputVolume->clone();

    if(mirrorX_.get()) {
        Volume* t = VolumeOperatorMirrorX::APPLY_OP(outputVolume);
        delete outputVolume;
        outputVolume = t;
    }
    if(mirrorY_.get()) {
        Volume* t = VolumeOperatorMirrorY::APPLY_OP(outputVolume);
        delete outputVolume;
        outputVolume = t;
    }
    if(mirrorZ_.get()) {
        Volume* t = VolumeOperatorMirrorZ::APPLY_OP(outputVolume);
        delete outputVolume;
        outputVolume = t;
    }

    // assign computed volume to outport
    if (outputVolume)
        outport_.setData(outputVolume);
    else
        outport_.setData(0);
}

}   // namespace
