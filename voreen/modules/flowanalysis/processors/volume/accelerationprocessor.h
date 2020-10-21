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

namespace voreen
{

class VRN_CORE_API AccelerationProcessor : public Processor {
public:
    AccelerationProcessor();

    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual std::string getCategory() const;
    virtual bool isReady() const;

    template<typename T>
    static void Process(const VolumeRAM_Mat3Float& jacobi, const VolumeAtomic<tgt::Vector3<T>>& velocity, VolumeAtomic<tgt::Vector3<T>>& outAcceleration);

protected:
    virtual void setDescriptions();
    virtual void process();

private:
    VolumePort inportJacobianVolume_;
    VolumePort inportVelocityVolume_;
    VolumePort outport_;
};

template<typename T>
void AccelerationProcessor::Process(const VolumeRAM_Mat3Float& jacobi, const VolumeAtomic<tgt::Vector3<T>>& velocity, VolumeAtomic<tgt::Vector3<T>>& outAcceleration) {
#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for( long i = 0; i < static_cast<long>( jacobi.getNumVoxels() ); ++i ) {
        outAcceleration.voxel( i ) = jacobi.voxel( i ) * velocity.voxel( i );
    }
}

} // namespace voreen

#endif // VRN_ACCELERATIONPROCESSOR_H
