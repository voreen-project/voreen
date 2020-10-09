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

#include "helicitydensity.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

HelicityDensity::HelicityDensity()
    : Processor()
    , velocityInport_(Port::INPORT, "velocity", "Velocity Input")
    , vorticityInport_(Port::INPORT, "vorticity", "Vorticity Input")
    , helicityDensityOutport_(Port::OUTPORT, "helicityDensity", "Helicity Density Output")
{
    addPort(velocityInport_);
    velocityInport_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(vorticityInport_);
    vorticityInport_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(helicityDensityOutport_);
}

HelicityDensity::~HelicityDensity() {}

Processor* HelicityDensity::create() const {
    return new HelicityDensity();
}

bool HelicityDensity::isReady() const {

    if(!velocityInport_.isReady()) {
        setNotReadyErrorMessage("Velocity Inport not ready");
        return false;
    }

    if(!vorticityInport_.isReady()) {
        setNotReadyErrorMessage("Vorticity Inport not ready");
        return false;
    }

    if(velocityInport_.getData()->getDimensions() != vorticityInport_.getData()->getDimensions()) {
        setNotReadyErrorMessage("Input dimensions must match");
        return false;
    }

    if(velocityInport_.getData()->getSpacing() != vorticityInport_.getData()->getSpacing()) {
        setNotReadyErrorMessage("Input spacing must match");
        return false;
    }

    if(velocityInport_.getData()->getOffset() != vorticityInport_.getData()->getOffset()) {
        setNotReadyErrorMessage("Input offset must match");
        return false;
    }

    return true;
}

void HelicityDensity::process() {

    VolumeRAMRepresentationLock velocity(velocityInport_.getData());
    RealWorldMapping velocityRwm = velocityInport_.getData()->getRealWorldMapping();
    VolumeRAMRepresentationLock vorticity(vorticityInport_.getData());
    RealWorldMapping vorticityRwm = vorticityInport_.getData()->getRealWorldMapping();

    auto volume = new VolumeRAM_Float(velocity->getDimensions());
    for(size_t i=0; i<volume->getNumVoxels(); i++) {
        // Calculate the dot product.
        float helicityDensity = 0.0f;
        for(size_t channel=0; channel<3; channel++) {
            helicityDensity += velocityRwm.normalizedToRealWorld(velocity->getVoxelNormalized(i, channel)) *
                               vorticityRwm.normalizedToRealWorld(vorticity->getVoxelNormalized(i, channel));
        }

        volume->voxel(i) = helicityDensity;
    }

    Volume* output = new Volume(volume, velocityInport_.getData());
    helicityDensityOutport_.setData(output);
}

}   // namespace
