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

#include "lp_measure.h"

#include <sstream>

namespace voreen {

bool LP_Measure::measureIntensity(const tgt::vec3 samplePos, std::vector<const VolumeBase* >& vvh, std::vector<plot_t>& intensity, fetchOption fo) {
    tgt::vec3 voxelPos; // voxel position
    int j = 0; //index for intensity vector j = i+1
    for (size_t i = 0; i < vvh.size(); ++i) {
        j++;
        voxelPos = vvh[i]->getWorldToVoxelMatrix() * samplePos;
        if ( (min(voxelPos) >= 0) && tgt::hand(tgt::lessThan(tgt::svec3(voxelPos), vvh[i]->getDimensions())) ) {
            //switch interpolationmode
            switch(fo) {
            case FO_NORMAL:
                intensity[j] = vvh[i]->getRepresentation<VolumeRAM>()->getVoxelNormalized(tgt::iround(voxelPos));
                break;
            case FO_LINEAR:
                intensity[j] = vvh[i]->getRepresentation<VolumeRAM>()->getVoxelNormalizedLinear(voxelPos);
                break;
            case FO_CUBIC:
                intensity[j] = vvh[i]->getRepresentation<VolumeRAM>()->getVoxelNormalizedCubic(voxelPos);
                break;
            }

            if(vvh[i]->hasMetaData("RealWorldMapping")) {
                const MetaDataBase* mdb = vvh[i]->getMetaData("RealWorldMapping");
                if(mdb) {
                    const RealWorldMappingMetaData* rwm = dynamic_cast<const RealWorldMappingMetaData*>(mdb);
                    if(rwm) {
                        intensity[j] = static_cast<plot_t>(rwm->getValue().normalizedToRealWorld(static_cast<float>(intensity[j])));
                    }
                }
            }
        }
        else
            //FIXME: if one volume is wrong, every volume will be not measured
            return false;
    }
    return true;
}

bool LP_Measure::calculate(const tgt::vec3 start, const tgt::vec3 end, float samplingRate, std::vector<const VolumeBase* >& vvh, fetchOption fo, PlotData* pData) {
    //return, if we have no volume
    if (vvh.size() == 0) return false;
//get the sampling step size
    float voxelSizeWorld = 9999.0f;
    for(size_t i = 0; i < vvh.size() ; i++) {
        tgt::ivec3 volDim = vvh.at(i)->getDimensions();
        tgt::vec3 cubeSizeWorld = vvh.at(i)->getCubeSize() * vvh.at(i)->getPhysicalToWorldMatrix().getScalingPart();
        float tVoxelSizeWorld = tgt::max(cubeSizeWorld / tgt::vec3(volDim));
            if (tVoxelSizeWorld < voxelSizeWorld)
                voxelSizeWorld = tVoxelSizeWorld;
    }
    //modify rate with properties
    float samplingStepSizeWorld = voxelSizeWorld / samplingRate;

//set new plotdata
    pData->reset(1,static_cast<int>(vvh.size()));
    tgt::vec3 link = end - start; //vector of the arrow
    float length = tgt::length(link);                            //length of arrow
    int steps = static_cast<int>(length / samplingStepSizeWorld +1);                //number of steps
    tgt::vec3 inc = link/tgt::vec3(static_cast<float>(steps)); float incLength = tgt::length(inc);
    std::vector<plot_t> intensity(vvh.size()+1); //result vector
    tgt::vec3 position = start-inc; float distanceScale = 0.0f;
    for(int i = 0; i <= steps ;i++) {
        //calculate pos
        position += inc;
        //calculate intensity
        if(measureIntensity(position, vvh, intensity, fo)) {
            intensity[0] = distanceScale;
            pData->insert(intensity);
            distanceScale += incLength;
        }
    }

    //set label + units
    std::string unit;
    pData->setColumnLabel(0,"x");
    for (size_t i = 0; i < vvh.size(); ++i) {
        std::stringstream ss;
        std::string title = vvh[i]->getMetaDataValue<StringMetaData, std::string>("Description", "");
        //get unit
        if(vvh[i]->hasMetaData("RealWorldMapping")) {
            const MetaDataBase* mdb = vvh[i]->getMetaData("RealWorldMapping");
            if(mdb) {
                const RealWorldMappingMetaData* rwm = dynamic_cast<const RealWorldMappingMetaData*>(mdb);
                if(rwm) {
                    unit = rwm->getValue().getUnit();
                }
            }
        }
        if (title.empty())
            ss << "Volume " << i+1;
        else
            ss << title;

        if(unit != ""){
            ss << " [" << unit << "]";
            unit = "";
        }

        pData->setColumnLabel(static_cast<int>(i+1),ss.str());
    }
    return true;
}

} // namespace voreen
