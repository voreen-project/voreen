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

#include "createVesselAroundPoints.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include <math.h>

using tgt::vec3;

namespace voreen {


const std::string createVesselAroundPoints::loggerCat_("voreen.vesselnetworkanalysisextra.createVesselAroundPoints");

createVesselAroundPoints::createVesselAroundPoints()
    : VolumeProcessor()

    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output",Processor::VALID)
    ,enableProcessing_("enableProcessing","enable processing",true)
    ,backgroundValue_("backgroundValue","backgroundValue",0,0,1.0)

{
    addPort(outport_);
    addPort(inport_);

    addProperty(enableProcessing_);
    addProperty(backgroundValue_);

}

createVesselAroundPoints::~createVesselAroundPoints() {
}

/*
 * finds max x-position of a blub in x direction
 */
void createVesselAroundPoints::findLimit(const VolumeRAM* input, int &maxPoint){

    for (size_t voxel_z=0; voxel_z<(input->getDimensions().z); voxel_z++) {
        for (size_t voxel_y=0; voxel_y<(input->getDimensions().y); voxel_y++) {
            for (size_t voxel_x=0; voxel_x<(input->getDimensions().x); voxel_x++) {

                tgt::ivec3 current(voxel_x,voxel_y,voxel_z);
                float value = input->getVoxelNormalized(current);
                if(value > 0.f){
                    if(voxel_x > static_cast<unsigned int>(maxPoint))
                        maxPoint = voxel_x;

                }

            }
        }
    }

}
/*
 * set all values from input volume to new target volume
 */
void createVesselAroundPoints::setVolume(const VolumeRAM* input, VolumeRAM_UInt8* target){

    VRN_FOR_EACH_VOXEL(p, tgt::svec3(0,0,0), input->getDimensions()) {
        float value = input->getVoxelNormalized(p);
        if(value > 0) {
            target->setVoxelNormalized(value,p);
        }
    }
}

/*
 * creates constant Object (area with constant radius) in Volume(radius,angle,z)
 */
void createVesselAroundPoints::createConstantObject(VolumeRAM_UInt8* input, size_t radius){
    for (size_t voxel_z=0; voxel_z<(input->getDimensions().z); voxel_z++) {
        for(size_t voxel_y=0;voxel_y<(input->getDimensions().y); voxel_y++){
            tgt::ivec3 p(radius,voxel_y,voxel_z);
            input->setVoxelNormalized(1.0f,p);
        }

    }
}

/*
void createHornVesselObject(VolumeRAM_UInt8* input, size_t radius){
    for (size_t voxel_z=0; voxel_z<(input->getDimensions().z); voxel_z++) {
        if(voxel_z%2==0)
            radius--;
        for(size_t voxel_y=0;voxel_y<(input->getDimensions().y); voxel_y++){
            tgt::ivec3 p(radius,voxel_y,voxel_z);
            input->setVoxelNormalized(0.5f,p);
        }

    }
}


void createconstrictionVesselObject(VolumeRAM_UInt8* input, size_t radius){
    for (size_t voxel_z=0; voxel_z<(input->getDimensions().z); voxel_z++) {
        if(voxel_z < input->getDimensions().z/2 && voxel_z%2==0)
            radius--;
        else if(voxel_z%2==0)
            radius++;
        for(size_t voxel_y=0;voxel_y<(input->getDimensions().y); voxel_y++){
            tgt::ivec3 p(radius,voxel_y,voxel_z);
            input->setVoxelNormalized(0.5f,p);
        }

    }
}


void createVesselObject(VolumeRAM_UInt8* input, size_t radius){

        for (size_t voxel_y=0; voxel_y<(input->getDimensions().y); voxel_y++) {
            tgt::vec3 p(radius,voxel_y,0);
            tgt::vec3 p1(radius-1,voxel_y,50);
            input->setVoxelNormalized(0.5f,p);
            input->setVoxelNormalized(0.5f,p1);
        }


}
*/

/*
 * creates constant noise in Volume
 */
void createVesselAroundPoints::setBackground(VolumeRAM_UInt8* input, float value){
    if(value > 0){
        for (size_t voxel_z=0; voxel_z<(input->getDimensions().z); voxel_z++) {
            for (size_t voxel_y=0; voxel_y<(input->getDimensions().y); voxel_y++) {
                for (size_t voxel_x=0; voxel_x<(input->getDimensions().x); voxel_x++) {
                    tgt::vec3 point(voxel_x,voxel_y,voxel_z);
                    if(input->getVoxelNormalized(point) == 0.f)
                        input->setVoxelNormalized(value,point);
                }
            }
        }


    }
}

/*

void createVesselAroundPoints::createCenterLine(VolumeRAM_UInt8* target, std::vector<tgt::vec3>& line){
    //int x_pos = target->getDimensions().x/2;
    //int y_pos = target->getDimensions().y/2;

    for(size_t zi = 0; zi < target->getDimensions().z; zi++){
        tgt::ivec3 middle(0,0,zi);
        line.push_back(middle);
    }

    for(auto const& point: line) {
        target->setVoxelNormalized(1.0f,point);
    }

}
 */

void createVesselAroundPoints::process() {

    const VolumeBase* invol = inport_.getData();
    const VolumeRAM* input = invol->getRepresentation<VolumeRAM>();
    tgt::svec3 dimensions = input->getDimensions();

    if (!enableProcessing_.get()) {
        outport_.setData(invol, false);
        return;
    }


    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dimensions);
    target->clear();

    //set target with input values
    setVolume(input,target);

    //find max x position for a blub -> to be sure that the cylinder is always around a point
    int maxPoint=0;
    findLimit(input,maxPoint);
    maxPoint+=1;
    if(static_cast<size_t>(maxPoint) >= input->getDimensions().x)
        maxPoint = input->getDimensions().x-1;

    //set background noise with given value
    setBackground(target,backgroundValue_.get()/100);

    /*std::vector<tgt::vec3> line;
    line.clear();
    createCenterLine(target,line);*/

    //create cylinder as constant area
    createConstantObject(target,maxPoint);

    //createHornVesselObject(target,maxPoint);
    //createconstrictionVesselObject(target,maxPoint);
    //createVesselObject(target,maxPoint);

    /*VolumeRAM_UInt8* bezierVol= new VolumeRAM_UInt8(dimensions);
    bezierVol->clear();

    std::vector<tgt::vec3> line;
    line.clear();

    tgt::ivec3 minPoints(dimensions.x-1,dimensions.y-1,dimensions.z-1);
    tgt::ivec3 maxPoints(0,0,0);
    findLimits(input,minPoints,maxPoints);


    tgt::ivec3 diameter(maxPoints.x - minPoints.x, maxPoints.y - minPoints.y, maxPoints.z - minPoints.z);
    float maxRadius = (std::max(diameter.x,diameter.y) / 2) + 1 ;


    createBezierLine(bezierVol,dimensions,line);
    cylinderAroundLine(target,dimensions,line,maxRadius);
    convertBezierToPolar(target,bezierVol);

    //vesselAroundPoints(target,maxPoint);
*/


    Volume* vh = new Volume(target, invol->getSpacing(), invol->getOffset() ); //target
    outport_.setData(vh);


}
VoreenSerializableObject* createVesselAroundPoints::create() const {
    return new createVesselAroundPoints();
}

} // namespace voreen
