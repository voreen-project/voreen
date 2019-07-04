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

#include "volumefloodfill.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

namespace voreen {

const std::string VolumeFloodFill::loggerCat_("voreen.vesselnetworkanalysisextra.volumefloodfill");

VolumeFloodFill::VolumeFloodFill()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , seedsInport_(Port::INPORT, "volumethinning.seeds.inport", "Seed Mask Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output", false, Processor::VALID /* do not update if something connects to this port*/)
    , enabledProp_("enabledProp","Enabled",true)
{
    addPort(inport_);
    addPort(seedsInport_);
    addPort(outport_);
    addProperty(enabledProp_);
}

VolumeFloodFill::~VolumeFloodFill() {
}

void VolumeFloodFill::process() {
    if(!enabledProp_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    const VolumeBase* invol = inport_.getData();
    const VolumeBase* seedvol = seedsInport_.getData(); //May be null
    if(!invol || !seedvol) {
        return;
    }
    if(invol->getDimensions() != seedvol->getDimensions()) {
        LERROR("Seed mask and input volume dimensions differ");
        tgtAssert(false, "Seed mask and input volume dimensions differ");
        return;
    }
    const VolumeRAM* input = invol->getRepresentation<VolumeRAM>();
    const VolumeRAM* seeds = seedvol->getRepresentation<VolumeRAM>();
    tgt::svec3 dim = invol->getDimensions();

    voreen::VolumeRAM* output = new VolumeAtomic<uint8_t>(dim);
    output->clear();

    std::vector<tgt::svec3> stack;
    VRN_FOR_EACH_VOXEL(p, tgt::svec3(0,0,0), dim) {
        if(seeds->getVoxelNormalized(p) > 0) {
            stack.push_back(p);
        }
    }
    while(!stack.empty()) {
        tgt::svec3 p = stack.back();
        stack.pop_back();

        if(output->getVoxelNormalized(p) == 0 && input->getVoxelNormalized(p) > 0) {
            for(int dz = -1; dz <= 1; ++dz) {
                for(int dy = -1; dy <= 1; ++dy) {
                    for(int dx = -1; dx <= 1; ++dx) {
                        tgt::svec3 q(
                                std::min(p.x + dx, dim.x-1),
                                std::min(p.y + dy, dim.y-1),
                                std::min(p.z + dz, dim.z-1)
                                );
                        if(p != q) {
                            stack.push_back(q);
                        }
                    }
                }
            }
            output->setVoxelNormalized(1.0f, p);
        }
    }

    VolumeBase* outvol = new Volume(output, invol->getSpacing(), invol->getOffset());

    outport_.setData(outvol);
}
VoreenSerializableObject* VolumeFloodFill::create() const {
    return new VolumeFloodFill();
}

} // namespace voreen
