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

#ifndef VRN_VOLUMEOPERATORMAGNITUDE_H
#define VRN_VOLUMEOPERATORMAGNITUDE_H

//#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "tgt/vector.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorMagnitude {
public:
    VolumeOperatorMagnitude(){};

    /**
     * Calculates vector magnitudes from a vector field volume.
     */
    template<typename U>
    Volume* apply(const VolumeBase* srcVolume);
};

//---------------------------------------------------------------------------------------------
//      apply function
//---------------------------------------------------------------------------------------------
template<typename U>
Volume* VolumeOperatorMagnitude::apply(const VolumeBase* srcVolume) {
    // get RAM representation of the input volume and relevant meta data
    VolumeRAMRepresentationLock input(srcVolume);
    if(!*input) {
        return nullptr;
    }

    RealWorldMapping rwm = srcVolume->getRealWorldMapping();
    tgt::svec3 dim = input->getDimensions();
    size_t numChannels = input->getNumChannels();

    // create the output volume
    VolumeAtomic<U>* result = new VolumeAtomic<U>(dim);

    // we need to find out the max vector magnitude of the input first
    float maxMagnitude = srcVolume->getDerivedData<VolumeMinMaxMagnitude>()->getMaxMagnitude();

    // now we can normalize magnitudes to the range [0, 1] for storing floats
    tgt::svec3 pos;
    for (pos.z = 0; pos.z < dim.z; pos.z++) {
        for (pos.y = 0; pos.y < dim.y; pos.y++) {
            for (pos.x = 0; pos.x < dim.x; pos.x++) {

                tgt::vec4 vector = tgt::vec4::zero;
                for(size_t channel = 0; channel < numChannels; channel++) {
                    vector[channel] = rwm.normalizedToRealWorld(input->getVoxelNormalized(pos, channel));
                }

                float vectorMagnitude = tgt::length(vector);

                // magnitude is always positive, so we can just normalize using the maximum
                result->setVoxelNormalized(vectorMagnitude / maxMagnitude, pos);
            }
        }
    }

    Volume* magnitudeVolume = new Volume(result, srcVolume);

    // Overwrite real-world mapping to rescale the normalized values
    RealWorldMapping rescaleMapping(maxMagnitude, 0.f, "");
    magnitudeVolume->setRealWorldMapping(rescaleMapping);

    return  magnitudeVolume;
}

} // namespace

#endif // VRN_VOLUMEOPERATORMAGNITUDE_H
