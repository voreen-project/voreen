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

#include "vectordecompose.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

const std::string VectorDecompose::loggerCat_("voreen.VectorDecompose");

VectorDecompose::VectorDecompose()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , radius_(Port::OUTPORT, "radius", "Radius", false)
    , phi_(Port::OUTPORT, "phi", "Phi (Zenith angle)", false)
    , theta_(Port::OUTPORT, "theta", "Theta (Azimuth angle)", false)
{
    addPort(inport_);
    //inport_.addCondition(new PortConditionVolumeChannelCount(3));
    inport_.addCondition(new PortConditionVolumeType3xFloat());
    addPort(radius_);
    addPort(phi_);
    addPort(theta_);
}

VectorDecompose::~VectorDecompose() {}

Processor* VectorDecompose::create() const {
    return new VectorDecompose();
}

void VectorDecompose::process() {

    const VolumeBase* inputHandle = inport_.getData();
    VolumeRAMRepresentationLock inputVolume(inputHandle);
    const VolumeRAM_3xFloat* vector = dynamic_cast<const VolumeRAM_3xFloat*>(*inputVolume);

    VolumeRAM_Float* radius = new VolumeRAM_Float(inputVolume->getDimensions());
    VolumeRAM_Float* phi = new VolumeRAM_Float(inputVolume->getDimensions());
    VolumeRAM_Float* theta = new VolumeRAM_Float(inputVolume->getDimensions());

    for(size_t i=0; i<inputVolume->getNumVoxels(); i++) {
        tgt::vec3 vec = vector->voxel(i);

        if(vec == tgt::vec3::zero) {
            radius->voxel(i) = 0.0f;
            phi->voxel(i) = 0.0f;//std::numeric_limits<float>::quiet_NaN();
            theta->voxel(i) = 0.0f;//std::numeric_limits<float>::quiet_NaN();
        }
        else {
            radius->voxel(i) = tgt::length(vec);
            phi->voxel(i) = tgt::rad2deg(std::acos(vec.y / radius->voxel(i)));
            theta->voxel(i) = tgt::rad2deg(std::atan2(vec.y, vec.x));
        }
    }

    radius_.setData(new Volume(radius, inputHandle));

    auto vol = new Volume(phi, inputHandle);
    vol->setRealWorldMapping(RealWorldMapping());
    vol->addDerivedData(new VolumeMinMax(0, 180, 0, 180));
    phi_.setData(vol);

    vol = new Volume(theta, inputHandle);
    vol->setRealWorldMapping(RealWorldMapping());
    vol->addDerivedData(new VolumeMinMax(-180, 180, -180, 180));
    theta_.setData(vol);
}

}   // namespace
