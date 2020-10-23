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
    auto* inputCondition = new PortConditionLogicalOr;
    inputCondition->addLinkedCondition(new PortConditionVolumeType3xFloat);
    inputCondition->addLinkedCondition(new PortConditionVolumeType3xDouble);
    inportVelocityVolume_.addCondition(inputCondition);

    addPort(outport_);
}

Processor *AccelerationProcessor::create() const {
    return new AccelerationProcessor();
}

std::string AccelerationProcessor::getClassName() const {
    return "AccelerationProcessor";
}

std::string AccelerationProcessor::getCategory() const {
    return "Volume Processing";
}

void AccelerationProcessor::setDescriptions() {
    setDescription("Computes acceleration volume by multiplying jacobian and velocity for each voxel");
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

void AccelerationProcessor::process() {

    VolumeRAMRepresentationLock jacobianVolume(inportJacobianVolume_.getData());
    const auto* jacobianVolumeData = dynamic_cast<const VolumeRAM_Mat3Float*>(*jacobianVolume);

    VolumeRAMRepresentationLock velocityVolume(inportVelocityVolume_.getData());

    std::unique_ptr<VolumeRAM> accelerationVolume( new VolumeRAM_3xFloat(jacobianVolume->getDimensions()) );

    if(velocityVolume->getFormat() == VolumeGenerator3xFloat().getFormat()) {
        using Type = VolumeRAM_3xFloat;
        const auto* velocityVolumeData = dynamic_cast<const Type*>(*velocityVolume);
        auto* output = new Type(jacobianVolumeData->getDimensions());
        Process( *jacobianVolumeData, *velocityVolumeData, *output );
        accelerationVolume.reset(output);
    }
    else if(velocityVolume->getFormat() == VolumeGenerator3xDouble().getFormat()) {
        using Type = VolumeRAM_3xDouble;
        const auto* velocityVolumeData = dynamic_cast<const Type*>(*velocityVolume);
        auto* output = new Type(jacobianVolumeData->getDimensions());
        Process( *jacobianVolumeData, *velocityVolumeData, *output );
        accelerationVolume.reset(output);
    }
    else {
        tgtAssert(false, "Unsupported format");
        return;
    }

    Volume* output = new Volume(accelerationVolume.release(), inportVelocityVolume_.getData()->getSpacing(), inportVelocityVolume_.getData()->getOffset());
    output->setMetaDataValue<StringMetaData>("name", "acceleration");
    outport_.setData(output);
}

} // namespace voreen
