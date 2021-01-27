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

#include "createTestVolume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include <boost/math/constants/constants.hpp>

using tgt::vec3;

namespace voreen {



const std::string createTestVolume::loggerCat_("voreen.vesselnetworkanalysisextra.createTestVolume");

createTestVolume::createTestVolume()
    : VolumeProcessor()
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output",Processor::VALID)
    , volX_("volX", "Volume XYZ Dimension", 100, 3, 1000)

    , stdevx_("stddevx", "standard deviation x", 1, 1, 25)
    , amount_("amount","amount of blubs",50,0,1000)
    , predeterminedSeed_("predeterminedSeed", "Seed", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , randomDevice()
{
    addPort(outport_);

    addProperty(volX_);
   /* ON_CHANGE_LAMBDA(volX_, [this] {
        meanx_.setMaxValue(int(volX_.get()/4));
    }); */

    addProperty(stdevx_);
    addProperty(amount_);

    addProperty(predeterminedSeed_);



}

createTestVolume::~createTestVolume() {
}


/*
 * help function for densityFunction
 */
 void createTestVolume::qMeanDistance(int &x, int mean, int &result){
     result = (x-mean)*(x-mean);
 }

/*
 * set normal density function
 */
 void createTestVolume::densityFunction(int x, int y, int z, int xMean, int yMean, int zMean, int stddev , float &value){

     int xValue, yValue, zValue;
     qMeanDistance(x,xMean,xValue);
     qMeanDistance(y,yMean,yValue);
     qMeanDistance(z,zMean,zValue);

     float stddev_2 = stddev*stddev;
     float factor = 1;//(1/((stddev_2)*2*M_PI));  // (1/((stddev_2)*2*boost::math::constants::pi<double>()));
     float eValue = - ( (xValue + yValue + zValue) / (2*stddev_2) );
     value = (1) * factor * std::exp(eValue);
 }

 /*
  * set different means for amount of given blubs
  */
 void createTestVolume::setMean(random_engine &randomEngine,std::uniform_int_distribution<uint32_t> dist,std::uniform_int_distribution<uint32_t> distYZ,std::vector<int>& meanX,std::vector<int>& meanY,std::vector<int>& meanZ , int iteration){

     for(int i=0;i < iteration ;i++){
         meanX[i] = dist(randomEngine) ;
         meanY[i] = distYZ(randomEngine) ;
         meanZ[i] = distYZ(randomEngine) ;
     }

 }

/*
 *  Create randon blubs in Volume
 */
 void createTestVolume::setRandomPoints(VolumeRAM* input, random_engine &randomEngine){

     //arrays with random values for x,y,z coordinates
     std::vector<int> xMeanArray(amount_.get(), 0);
     std::vector<int> yMeanArray(amount_.get(), 0);
     std::vector<int> zMeanArray(amount_.get(), 0);
     float value;
     int meanX = 0;
     int meanY = 0;
     int meanZ = 0;


     //uniform int dist for random mean values
     //dist for y,z values from 0 to max size
     //dist for x values from 0 to 1/2 max size -> Else the new created cylinder wont fit into the volume
     std::uniform_int_distribution<uint32_t> distX(0,int(volX_.get()*1/2)); // int(volX_.get()*1/3),int(volX_.get()*2/3));
     std::uniform_int_distribution<uint32_t> distYZ(0,int(volX_.get()-1));
     setMean(randomEngine,distX,distYZ,xMeanArray,yMeanArray,zMeanArray,amount_.get());

     //loop for a given amount of blubs
     for(int i=0;i < amount_.get();i++){
         meanX = xMeanArray[i];
         meanY = yMeanArray[i];
         meanZ = zMeanArray[i];

         for (size_t voxel_z=0; voxel_z<(input->getDimensions().z); voxel_z++) {
             for (size_t voxel_y=0; voxel_y<(input->getDimensions().y); voxel_y++) {
                 for (size_t voxel_x=0; voxel_x<(input->getDimensions().x); voxel_x++) {

                         //set random values with normal density function
                         densityFunction(voxel_x,voxel_y,voxel_z,meanX,meanY,meanZ, stdevx_.get(),value);
                         tgt::vec3 point(voxel_x,voxel_y,voxel_z);
                         if(value > input->getVoxelNormalized(point) ){
                            input->setVoxelNormalized(value,point);
                         }

                 }
             }
        }
     }




}


void createTestVolume::process() {

    //Set dimensions
    tgt::ivec3 dimensions(volX_.get(),volX_.get(),volX_.get());

    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dimensions);

    VolumeRAM* outputVolume = target;
    outputVolume->clear();

    //create random engine and fill with seed
    random_engine randomEngine(randomDevice());
    randomEngine.seed(predeterminedSeed_.get());

    /*
    switch (vesselType_.get()) {
    case 1:
        createConstantObject(target,maxPoint);
        break;
    case 2:
        createBezierObject(target,maxPoint);
    default:
        break;
    } */

    setRandomPoints(target,randomEngine);

    Volume* vh = new Volume(target, vec3(1.0f), vec3(0.0f));
    outport_.setData(vh);


}
VoreenSerializableObject* createTestVolume::create() const {
    return new createTestVolume();
}

} // namespace voreen



