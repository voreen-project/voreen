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

#include "accelerationprocessor.h"

#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

AccelerationProcessor::AccelerationProcessor()
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

Processor* AccelerationProcessor::create() const {
    return new AccelerationProcessor();
}


bool AccelerationProcessor::isReady() const {
    bool ready = Processor::isReady();
    if(!ready) {
        return false;
    }

    if(inportJacobianVolume_.getData()->getDimensions() != inportVelocityVolume_.getData()->getDimensions()) {
        setNotReadyErrorMessage("Input dimensions must match");
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

void AccelerationProcessor::Process(const VolumeRAM_Mat3Float& jacobi, const VolumeRAM& velocity, VolumeRAM_3xFloat& outAcceleration, RealWorldMapping rwm) {

    // If volume contains no real world mapping, we can use optimized code.
    if (rwm == RealWorldMapping()) {
        if (velocity.getFormat() == VolumeGenerator3xFloat().getFormat()) {
            optimizedProcess(jacobi, dynamic_cast<const VolumeRAM_3xFloat&>(velocity), outAcceleration);
            return;
        }
        else if (velocity.getFormat() == VolumeGenerator3xDouble().getFormat()) {
            optimizedProcess(jacobi, dynamic_cast<const VolumeRAM_3xDouble&>(velocity), outAcceleration);
            return;
        }
    }

#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for (long i = 0; i < static_cast<long>(jacobi.getNumVoxels()); ++i) {
        tgt::vec3 voxel;
        for (size_t channel = 0; channel < 3; channel++) {
            voxel[channel] = rwm.normalizedToRealWorld(velocity.getVoxelNormalized(i, channel));
        }

        // Perform the multiplication.
        voxel = jacobi.voxel(i) * voxel;

        for (size_t channel = 0; channel < 3; channel++) {
            outAcceleration.setVoxelNormalized(rwm.realWorldToNormalized(voxel[channel]), i, channel);
        }
    }
}

void AccelerationProcessor::process() {

    VolumeRAMRepresentationLock jacobianVolume(inportJacobianVolume_.getData());
    const auto* jacobianVolumeData = dynamic_cast<const VolumeRAM_Mat3Float*>(*jacobianVolume);

    const VolumeBase* velocity = inportVelocityVolume_.getData();
    VolumeRAMRepresentationLock velocityVolume(velocity);

    RealWorldMapping rwm = velocity->getRealWorldMapping();
    
    std::unique_ptr<VolumeRAM_3xFloat> accelerationVolume(new VolumeRAM_3xFloat(velocity->getDimensions()));
    Process(*jacobianVolumeData, **velocityVolume, *accelerationVolume, rwm);

    Volume* output = new Volume(accelerationVolume.release(), velocity->getSpacing(), velocity->getOffset());
    output->setRealWorldMapping(velocity->getRealWorldMapping());
    output->setMetaDataValue<StringMetaData>("name", "acceleration");
    outport_.setData(output);
}

} // namespace voreen
