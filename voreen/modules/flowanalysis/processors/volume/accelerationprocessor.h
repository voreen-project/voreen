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

#ifndef VRN_ACCELERATIONPROCESSOR_H
#define VRN_ACCELERATIONPROCESSOR_H

//add base class header
#include "voreen/core/processors/processor.h"

//add used port headers
#include "voreen/core/ports/volumeport.h"

//add used datastructure headers
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class VRN_CORE_API AccelerationProcessor : public Processor {
public:
    AccelerationProcessor();

    virtual Processor* create() const { return new AccelerationProcessor(); }
    virtual std::string getClassName() const { return "AccelerationProcessor"; }
    virtual std::string getCategory() const { return "Volume Processing"; }
    
    virtual bool isReady() const;

    static void Process(const VolumeRAM_Mat3Float& jacobi, const VolumeRAM& velocity, VolumeRAM_3xFloat& outAcceleration, RealWorldMapping rwm = RealWorldMapping());

protected:
    virtual void setDescriptions() {
        setDescription("Computes acceleration volume by multiplying jacobian and velocity for each voxel");
    }

    virtual void process();

private:
    VolumePort inportJacobianVolume_;
    VolumePort inportVelocityVolume_;
    VolumePort outport_;
};


} // namespace voreen

#endif // VRN_ACCELERATIONPROCESSOR_H
