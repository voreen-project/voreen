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

#include "acceleration.h"

#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

Acceleration::Acceleration()
    : Processor()
    , inportJacobianVolume_(Port::INPORT, "inportJacobianVolume", "Jacobian")
    , inportVelocityVolume_(Port::INPORT, "inportVelocityVolume", "Velocity")
    , outport_(Port::OUTPORT, "outport", "3-Float Vector Volume (Acceleration)")
{
    addPort(inportJacobianVolume_);
    inportJacobianVolume_.addCondition(new PortConditionVolumeType("Matrix3(float)", "Volume_Mat3Float"));
    addPort(inportVelocityVolume_);
    addPort(outport_);
}

bool Acceleration::isReady() const {
    bool ready = Processor::isReady();
    if(!ready) {
        return false;
    }

    if(inportJacobianVolume_.getData()->getDimensions() != inportVelocityVolume_.getData()->getDimensions()) {
        setNotReadyErrorMessage("Input dimensions must match");
        return false;
    }

    if(inportJacobianVolume_.getData()->getSpacing() != inportVelocityVolume_.getData()->getSpacing()) {
        setNotReadyErrorMessage("Input spacing must match");
        return false;
    }

    if(inportJacobianVolume_.getData()->getOffset() != inportVelocityVolume_.getData()->getOffset()) {
        setNotReadyErrorMessage("Input offset must match");
        return false;
    }

    return true;
}

template<typename T>
void optimizedProcess(const VolumeRAM_Mat3Float& jacobi, const VolumeAtomic<tgt::Vector3<T>>& velocity, VolumeRAM_3xFloat& outAcceleration) {
#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for (long i = 0; i < static_cast<long>(jacobi.getNumVoxels()); ++i) {
        outAcceleration.voxel(i) = jacobi.voxel(i) * velocity.voxel(i);
    }
}

void Acceleration::Process(const VolumeRAM_Mat3Float& jacobian, const VolumeRAM& velocity, VolumeRAM_3xFloat& outAcceleration, RealWorldMapping rwm) {

    // If volume contains default real world mapping, we can use optimized code.
    if (rwm == RealWorldMapping()) {
        if (velocity.getFormat() == VolumeGenerator3xFloat().getFormat()) {
            optimizedProcess(jacobian, dynamic_cast<const VolumeRAM_3xFloat&>(velocity), outAcceleration);
            return;
        }
        else if (velocity.getFormat() == VolumeGenerator3xDouble().getFormat()) {
            optimizedProcess(jacobian, dynamic_cast<const VolumeRAM_3xDouble&>(velocity), outAcceleration);
            return;
        }
    }

#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for (long i = 0; i < static_cast<long>(jacobian.getNumVoxels()); ++i) {
        tgt::vec3 voxel;

        // Retrieve real world values before multiplication.
        for (size_t channel = 0; channel < tgt::vec3::size; channel++) {
            voxel[channel] = rwm.normalizedToRealWorld(velocity.getVoxelNormalized(i, channel));
        }

        // Perform the multiplication.
        voxel = jacobian.voxel(i) * voxel;

        // Since we output a float volume anyways, we discard the real world mapping.
        for (size_t channel = 0; channel < tgt::vec3::size; channel++) {
            outAcceleration.setVoxelNormalized(voxel[channel], i, channel);
        }
    }
}

void Acceleration::process() {

    VolumeRAMRepresentationLock jacobianVolume(inportJacobianVolume_.getData());
    const auto* jacobianVolumeData = dynamic_cast<const VolumeRAM_Mat3Float*>(*jacobianVolume);

    const VolumeBase* velocity = inportVelocityVolume_.getData();
    VolumeRAMRepresentationLock velocityVolume(velocity);

    std::unique_ptr<VolumeRAM_3xFloat> accelerationVolume(new VolumeRAM_3xFloat(velocity->getDimensions()));
    Process(*jacobianVolumeData, **velocityVolume, *accelerationVolume, velocity->getRealWorldMapping());

    Volume* output = new Volume(accelerationVolume.release(), velocity->getSpacing(), velocity->getOffset());
    output->setRealWorldMapping(RealWorldMapping()); // Set default real world mapping.
    output->setModality(Modality("acceleration"));
    outport_.setData(output);
}

} // namespace voreen
