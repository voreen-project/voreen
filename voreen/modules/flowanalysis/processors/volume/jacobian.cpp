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

#include "jacobian.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

const std::string Jacobian::loggerCat_("voreen.flowanalysis.Jacobian");

Jacobian::Jacobian()
    : Processor()
    , inputVolume_(Port::INPORT, "Jacobian.inputVolume", "3D Vector Volume Input")
    , outputJacobian_(Port::OUTPORT, "Jacobian.outputJacobian", "Jacobian matrix")
{
    addPort(inputVolume_);
    inputVolume_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(outputJacobian_);
}

Processor* Jacobian::create() const {
    return new Jacobian();
}

void Jacobian::process() {

    auto inputVolume = inputVolume_.getData();
    VolumeRAMRepresentationLock inputVolumeData(inputVolume);
    tgt::Vector3<long> dimensions = inputVolumeData->getDimensions();
    tgt::vec3 spacing = inputVolume->getSpacing();
    RealWorldMapping rwm = inputVolume->getRealWorldMapping();

    auto jacobianVolume = new VolumeRAM_Mat3Float(inputVolumeData->getDimensions());

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (long z = 0; z < dimensions.z; z++) {
        for (long y = 0; y < dimensions.y; y++) {
            for (long x = 0; x < dimensions.x; x++) {
                tgt::svec3 pos(x, y, z);

                tgt::mat3 jacobian(
                        VolumeOperatorGradient::calcGradientCentralDifferences(*inputVolumeData, spacing, pos, 0),
                        VolumeOperatorGradient::calcGradientCentralDifferences(*inputVolumeData, spacing, pos, 1),
                        VolumeOperatorGradient::calcGradientCentralDifferences(*inputVolumeData, spacing, pos, 2)
                        );

                // Apply real world mapping to each element.
                for(float& i : jacobian.elem) {
                    i = rwm.normalizedToRealWorld(i);
                }

                jacobianVolume->voxel(pos) = jacobian;
            }
        }
    }

    auto* volume = new Volume(jacobianVolume, inputVolume);
    volume->setRealWorldMapping(RealWorldMapping()); // Override to default rwm.
    volume->setModality(Modality("jacobian"));
    outputJacobian_.setData(volume);
}

} // namespace voreen
