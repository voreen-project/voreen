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

#ifndef VRN_ACCELERATIONPROCESSOR_H
#define VRN_ACCELERATIONPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class VRN_CORE_API Acceleration : public Processor {
public:
    Acceleration();

    virtual Processor* create() const { return new Acceleration(); }
    virtual std::string getClassName() const { return "Acceleration"; }
    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }
    virtual bool isReady() const;

    static void Process(const VolumeRAM_Mat3Float& jacobian, const VolumeRAM& velocity, VolumeRAM_3xFloat& outAcceleration, RealWorldMapping rwm = RealWorldMapping());

protected:
    virtual void setDescriptions() {
        setDescription("Computes acceleration volume by voxel-wise multiplication of jacobian and velocity volume");
    }

    virtual void process();

private:
    VolumePort inportJacobianVolume_;
    VolumePort inportVelocityVolume_;
    VolumePort outport_;
};


} // namespace voreen

#endif // VRN_ACCELERATIONPROCESSOR_H
