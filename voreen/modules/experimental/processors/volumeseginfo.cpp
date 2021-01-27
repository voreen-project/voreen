/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "volumeseginfo.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <sstream>

namespace voreen {

const std::string VolumeSegInfo::loggerCat_("voreen.experimental.VolumeSegInfo");

VolumeSegInfo::VolumeSegInfo()
    : VolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input"),
    outport_(Port::OUTPORT, "textport.output")

{
    addPort(inport_);
    addPort(outport_);
}

VolumeSegInfo::~VolumeSegInfo() {}

void VolumeSegInfo::process() {
    const VolumeRAM_UInt8* inputVolume = dynamic_cast<const VolumeRAM_UInt8*>(inport_.getData()->getRepresentation<VolumeRAM>());
    std::stringstream outputString;

    tgt::ivec3 dimensions = inputVolume->getDimensions();

    std::map<uint8_t, std::vector<tgt::ivec3> > segments = std::map<uint8_t, std::vector<tgt::ivec3> >();

    for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
        for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
            for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {

                tgt::ivec3 pos = tgt::ivec3(voxel_x, voxel_y, voxel_z);
                uint8_t v = inputVolume->voxel(pos);
                std::map<uint8_t, std::vector<tgt::ivec3> >::iterator it = segments.find(v);
                if(it == segments.end()) {
                    std::vector<tgt::ivec3> positions = std::vector<tgt::ivec3>();
                    positions.push_back(pos);
                    segments.insert(std::make_pair(v, positions));
                } else {
                    it->second.push_back(pos);
                }
            }
        }
    }

    outputString << "No of Segments: " << segments.size() << std::endl;

    for(std::map<uint8_t, std::vector<tgt::ivec3> >::iterator it = segments.begin(); it != segments.end(); it++) {
        tgt::dvec3 avg = tgt::dvec3(0.0);
        std::vector<tgt::ivec3> positions = it->second;
        for(size_t i = 0; i < positions.size(); i++) {
            avg += tgt::dvec3(positions[i]);
        }
        outputString << "Segment " << (int)it->first << ": " << std::endl;
        outputString << "\tSize " << positions.size() << std::endl;
        outputString << "\tPercent " << (100.0 * (double)positions.size()) / (double)(dimensions.x*dimensions.y*dimensions.z) << std::endl;
        outputString << "\tAvg Pos " << avg / (double)positions.size() << std::endl;
        //outputString << "\tnormalized " << 2.0*(avg / (double)positions.size()) / tgt::max(tgt::dvec3(dimensions)) - tgt::dvec3(1.0) << std::endl;
        outputString << "\tnormalized " << 2.0*(avg / (double)positions.size()) / tgt::dvec3(dimensions) - tgt::dvec3(1.0) << std::endl;
    }

    // assign computed volume to outport
    outport_.setData(outputString.str());
}

}   // namespace
