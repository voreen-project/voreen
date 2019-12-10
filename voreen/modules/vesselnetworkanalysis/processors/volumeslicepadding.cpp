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

#include "volumeslicepadding.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

namespace voreen {

const std::string VolumeSlicePadding::loggerCat_("voreen.vesselnetworkanalysis.volumefloodfill");

VolumeSlicePadding::VolumeSlicePadding()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output", false, Processor::VALID /* do not update if something connects to this port*/)
    , enabledProp_("enabledProp","Enabled",true)
    , sliceVoxelValue_("sliceVoxelValue", "Voxel value of additional slices", 0, 0, 1)
    , numSlicesBefore_("numSlicesBefore", "Number of slices to prepend", 1, 0, 100)
    , numSlicesAfter_("numSlicesAfter", "Number of slices to append", 1, 0, 100)
{
    addPort(inport_);
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(sliceVoxelValue_);
    addProperty(numSlicesBefore_);
    addProperty(numSlicesAfter_);
}

VolumeSlicePadding::~VolumeSlicePadding() {
}

void VolumeSlicePadding::process() {
    if(!enabledProp_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    const VolumeBase* invol = inport_.getData();
    if(!invol) {
        return;
    }

    const size_t numSlicesBefore = numSlicesBefore_.get();
    const size_t numSlicesAfter = numSlicesAfter_.get();
    const float sliceVoxelValue = sliceVoxelValue_.get();

    const VolumeRAM* input = invol->getRepresentation<VolumeRAM>();
    tgt::svec3 dim = invol->getDimensions();
    tgt::svec3 outputdim = dim;
    outputdim.z += numSlicesBefore + numSlicesAfter;

    voreen::VolumeRAM* output = VolumeFactory().create(input->getFormat(), outputdim);


    for(size_t z = 0; z < numSlicesBefore; ++z) {
        for(size_t y = 0; y < outputdim.y; ++y) {
            for(size_t x = 0; x < outputdim.x; ++x) {
                output->setVoxelNormalized(sliceVoxelValue, x, y, z);
            }
        }
    }
    //VRN_FOR_EACH_VOXEL(p, tgt::svec3(0,0,0), dim) {
    //    float val = input->getVoxelNormalized(p);
    //    output->setVoxelNormalized(val, p.x, p.y, p.z+numSlicesBefore);
    //}
    for(size_t z = 0; z < dim.z; ++z) {
        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                float val = input->getVoxelNormalized(x, y, z);
                output->setVoxelNormalized(val, x, y, z+numSlicesBefore);
            }
        }
    }
    for(size_t z = dim.z + numSlicesBefore; z < outputdim.z; ++z) {
        for(size_t y = 0; y < outputdim.y; ++y) {
            for(size_t x = 0; x < outputdim.x; ++x) {
                output->setVoxelNormalized(sliceVoxelValue, x, y, z);
            }
        }
    }

    Volume* outvol = new Volume(output, invol);
    tgt::vec3 offset = invol->getOffset();
    offset.z -= invol->getSpacing().z*numSlicesBefore;
    outvol->setOffset(offset);

    outport_.setData(outvol);
}
VoreenSerializableObject* VolumeSlicePadding::create() const {
    return new VolumeSlicePadding();
}

} // namespace voreen
