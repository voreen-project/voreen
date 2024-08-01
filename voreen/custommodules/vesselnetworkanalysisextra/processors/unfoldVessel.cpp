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

#include "unfoldVessel.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include <boost/math/constants/constants.hpp>

using tgt::vec3;

namespace voreen {

const std::string unfoldVessel::loggerCat_("voreen.vesselnetworkanalysisextra.unfoldVessel");
unfoldVessel::unfoldVessel()
    : VolumeProcessor()

    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output",Processor::VALID)

{
    addPort(outport_);
    addPort(inport_);


}

unfoldVessel::~unfoldVessel() {
}


void mapToPolar(const VolumeRAM* input, VolumeRAM_UInt8* target){

    tgt::vec3 dimensions = input->getDimensions();
    tgt::vec2 center(dimensions.x/2,dimensions.y/2);
    double factor = (double)input->getDimensions().y / 360.0;

    for (size_t voxel_z=0; voxel_z<(input->getDimensions().z); voxel_z++) {
        for (size_t voxel_y=0; voxel_y<(input->getDimensions().y); voxel_y++) {
            for (size_t voxel_x=0; voxel_x<(input->getDimensions().x); voxel_x++) {

                tgt::vec3 point(voxel_x,voxel_y,voxel_z);
                float value = input->getVoxelNormalized(point);
                if(value > 0){

                    tgt::vec2 xyPoint(point.x,point.y);
                    xyPoint = xyPoint - center;

                    int radius = (int) std::sqrt(xyPoint.x*xyPoint.x + xyPoint.y*xyPoint.y);
                    double radian = std::atan2(xyPoint.x ,xyPoint.y ); //0-360
                    if(radian <= 0)
                        radian+=2*boost::math::constants::pi<double>();
                    double roh = radian *(180.0/boost::math::constants::pi<double>());



                    roh = (int) roh * factor;

                    if(roh >= dimensions.y)
                        roh--;

                    tgt::vec3 newPoint(radius,roh,voxel_z);
                    target->setVoxelNormalized(value,newPoint);
                }
            }
        }
    }

}



void unfoldVessel::process() {

    const VolumeBase* invol = inport_.getData();
    const VolumeRAM* input = invol->getRepresentation<VolumeRAM>();

    tgt::ivec3 dimensions = input->getDimensions();
    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dimensions);

    VolumeRAM* outputVolume = target;
    outputVolume->clear();

    mapToPolar(input,target);


    Volume* vh = new Volume(target, invol->getSpacing(), invol->getOffset() );
    outport_.setData(vh);

}
VoreenSerializableObject* unfoldVessel::create() const {
    return new unfoldVessel();
}

} // namespace voreen
