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

#include "volumenormalization.h"
#include "../operators/volumeoperatornormalize.h"

namespace voreen {

VolumeNormalization::VolumeNormalization()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volume.inport")
    , outport_(Port::OUTPORT, "volume.outport", "volume.outport", true)
    , enableProcessing_("enableProcessing", "Enable")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
}

void VolumeNormalization::process() {
    const VolumeBase* inputVolume = inport_.getData();

    if (!enableProcessing_.get()) {
        outport_.setData(inputVolume, false);
        return;
    }

    Volume* outputVolume = inputVolume->clone();

    Volume* t = VolumeOperatorNormalize::APPLY_OP(outputVolume);
    delete outputVolume;
    outputVolume = t;

    if (outputVolume)
        outport_.setData(outputVolume);
    else
        outport_.setData(0);
}

} // namespace
