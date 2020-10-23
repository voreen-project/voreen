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

#include "curlprocessor.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include <algorithm>


namespace voreen {

CurlProcessor::CurlProcessor()
    : Processor()
    , inport_(Port::INPORT, "inport", "Jacobian")
    , outport_(Port::OUTPORT, "outport", "Curl Vector Field")
{
    inport_.addCondition(new PortConditionVolumeType("Matrix3(float)", "VolumeRAM_Mat3Float"));
    addPort(inport_);
    addPort(outport_);
}

Processor* CurlProcessor::create() const {
    return new CurlProcessor();
}

void CurlProcessor::Process( const VolumeRAM_Mat3Float& jacobi, VolumeRAM_3xFloat& outCurl ) {
    const auto dim = jacobi.getDimensions();

#pragma omp parallel for
    for (long x = 0; x < static_cast<long>(dim.x); ++x) {
        for (auto y = 0; y < dim.y; ++y) {
            for (auto z = 0; z < dim.z; ++z) {
                auto a = jacobi.voxel(x, y, z)[2][1] - jacobi.voxel(x, y, z)[1][2];
                auto b = jacobi.voxel(x, y, z)[0][2] - jacobi.voxel(x, y, z)[2][0];
                auto c = jacobi.voxel(x, y, z)[1][0] - jacobi.voxel(x, y, z)[0][1];
                outCurl.voxel(x, y, z) = tgt::vec3(a, b, c);
            }
        }
    }
}

void CurlProcessor::process() {

    auto input = inport_.getData();
    auto jacobian = dynamic_cast<const VolumeRAM_Mat3Float*>(input->getRepresentation<VolumeRAM>()); // VolumeRAM_Mat3Float
    auto curl = new VolumeRAM_3xFloat(jacobian->getDimensions());
    auto dim = jacobian->getDimensions();

    CurlProcessor::Process(*jacobian, *curl);

    Volume* output = new Volume(curl, input->getSpacing(), input->getOffset());
    output->setMetaDataValue<StringMetaData>("name", "curl");
    outport_.setData(output);
}

} // namespace voreen
