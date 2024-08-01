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

#include "materialderivative.h"

namespace voreen {

const std::string MaterialDerivative::loggerCat_ = "MaterialDerivative";

MaterialDerivative::MaterialDerivative() 
    : Processor()
    , v0_( Port::INPORT, "v0", "Volume v for timestep t-1" )
    , v1_( Port::INPORT, "v1", "Volume v for timestep t" )
    , velocity_( Port::INPORT, "flowfield", "3-component flow field" )
    , jacobian_( Port::INPORT, "jacobian_of_v", "Jacobian of the volume v for timestep t" )
    , outport_( Port::OUTPORT, "material_derivative", "Material derivative of v" )
{

    addPort(v0_);
    addPort(v1_);
    addPort(velocity_);
    addPort(jacobian_);

    addPort(outport_);
}

void MaterialDerivative::process() {
    auto v0PortData = v0_.getData();
    auto v1PortData = v1_.getData();
    auto velocityPortData = velocity_.getData();
    auto jacobianPortData = jacobian_.getData();
    
    VolumeRAMRepresentationLock v0Volume(v0PortData);
    VolumeRAMRepresentationLock v1Volume(v1PortData);
    VolumeRAMRepresentationLock velocityVolume(velocityPortData);
    VolumeRAMRepresentationLock jacobianVolume(jacobianPortData);

    const auto v0VolumeData = dynamic_cast<const VolumeRAM_3xFloat*>(*v0Volume);
    const auto v1VolumeData = dynamic_cast<const VolumeRAM_3xFloat*>(*v1Volume);
    const auto velocityVolumeData = dynamic_cast<const VolumeRAM_3xFloat*>(*velocityVolume);
    const auto jacobianVolumeData = dynamic_cast<const VolumeRAM_Mat3Float*>(*jacobianVolume);

    if (!v0VolumeData) {
        throw std::runtime_error("Volume v0 not connected properly");
    }
    if (!v1VolumeData) {
        throw std::runtime_error("Volume v1 not connected properly");
    }
    if (!velocityVolumeData) {
        throw std::runtime_error("Flowfield not connected properly");
    }
    if (!jacobianVolumeData) {
        throw std::runtime_error("Jacobian not connected properly");
    }

    tgt::Vector3<long> dimensions = v0Volume->getDimensions();

    auto materialDerivative = new VolumeRAM_3xFloat(dimensions);

    // For each voxel, compute 
    // Dv/Dt = dv/dt + u * grad v

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (long z = 0; z < dimensions.z; z++) {
        for (long y = 0; y < dimensions.y; y++) {
            for (long x = 0; x < dimensions.x; x++) {
                tgt::svec3 pos(x, y, z);
                
                // Approximate derivative by looking at neighbouring timesteps
                auto dvdt = v1VolumeData->voxel(pos) - v0VolumeData->voxel(pos);

                auto u = velocityVolumeData->voxel(pos);

                auto gradv = jacobianVolumeData->voxel(pos);

                materialDerivative->voxel(pos) = dvdt + u * gradv;
            }
        }
    }

    auto* volume = new Volume(materialDerivative, v0PortData);
    volume->setRealWorldMapping(RealWorldMapping()); // Override to default rwm.
    volume->setModality(Modality("pyrogenic_mapping"));
    outport_.setData(volume);
}

}
