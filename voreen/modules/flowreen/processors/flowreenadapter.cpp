/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "flowreenadapter.h"
#include "modules/flowreen/datastructures/deprecated/volumeflow3d.h"
#include "modules/flowreen/datastructures/deprecated/volumeoperatorflowmagnitude.h"

#include "voreen/core/datastructures/volume/volume.h"

namespace voreen {

FlowreenAdapter::FlowreenAdapter()
    : Processor(),
    processedVolume_(0),
    volInport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    volOutport_(Port::OUTPORT, "volumehandle.output", "Volume Output", true)
{
    addPort(volInport_);
    addPort(volOutport_);
}

FlowreenAdapter::~FlowreenAdapter() {
    /*if ((processedVolume_ != 0) && (processedVolume_ != currentVolume_))
        delete processedVolume_;*/
}

void FlowreenAdapter::process() {
    currentVolume_ = volInport_.getData();
    if (dynamic_cast<const VolumeRAM_3xFloat*>(currentVolume_->getRepresentation<VolumeRAM>()) != 0)
        calculateMagnitudes();
    else
        LERROR("supplied Volume seems to contain no flow data! Cannot proceed.");
}

void FlowreenAdapter::initialize() {
    FlowreenProcessor::init();
    Processor::initialize();
}

// private methods
//

void FlowreenAdapter::calculateMagnitudes() {
    const VolumeRAM_3xFloat* input = dynamic_cast<const VolumeRAM_3xFloat*>(currentVolume_->getRepresentation<VolumeRAM>());
    VolumeRAM_3xFloat* output = 0;
    if(input) {
        tgt::ivec3 dimensions = input->getDimensions();
        const size_t numVoxels = dimensions.x * dimensions.y * dimensions.z;

        tgt::vec3* voxels = 0;
        try {
            voxels = new tgt::vec3[numVoxels];
        } catch (std::bad_alloc) {
            throw;
        }

        const tgt::ivec3 xyz(0, 1, 2);  // target voxel order is XYZ

        tgt::ivec3 pos(0, 0, 0);
        float maxValue = std::numeric_limits<float>::min();
        float minValue = std::numeric_limits<float>::max();
        float maxMagnitude = 0.0f;

        for (size_t i = 0; i < numVoxels; ++i) {
            voxels[i] = input->voxel(i);

            // calculate max and min of current voxel and dataset's min and max.
            float voxelMax = tgt::max(voxels[i]);
            float voxelMin = tgt::min(voxels[i]);
            float magnitude = tgt::length(voxels[i]);

            if (voxelMax > maxValue)
                maxValue = voxelMax;
            if (voxelMin < minValue)
                minValue = voxelMin;
            if (magnitude > maxMagnitude)
                maxMagnitude = magnitude;
        }

        LINFO("min. among all components of all voxels: " << minValue);
        LINFO("max. among all components of all voxels: " << maxValue);
        LINFO("max. magnitude among all  voxels: " << maxMagnitude);

        output = new VolumeFlow3D(voxels, dimensions, minValue, maxValue, maxMagnitude);
    }

    /*if ((processedVolume_ != 0)
        && processedVolume_ != currentVolume_)
    {
            delete processedVolume_;
    }*/

    if (output != 0)
        processedVolume_ = new Volume(output, currentVolume_);
    else
        processedVolume_ = 0;

    volOutport_.setData(processedVolume_);
}

}   // namespace
