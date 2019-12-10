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

#ifndef VRN_LP_MEASURE_H
#define VRN_LP_MEASURE_H

//volumehandle
#include "voreen/core/datastructures/volume/volume.h"
//plotdata
#include "modules/plotting/datastructures/plotdata.h"
//vectors
#include "tgt/vector.h"
#include <vector>

namespace voreen {

/**
 * Utility class beeing used in lineprofile.h
 * This class encapsulate the functions working on the volume itself.
 */
class LP_Measure {
public:
    //enum of the voxel fetch
    enum fetchOption {
        FO_NORMAL,  //< fetch nearest
        FO_LINEAR,  //< fetch linear
        FO_CUBIC    //< fetch cubic
    };

    /**
     * Function for calculating the plotdata.
     * This function measures the intensity of the voxels from start to end in all volumes in vhh and stores the results in outData.
     *
     * @param start         beginning point of the calculation (worldcoordinates)
     * @param end           end point of the calculation (worldcoordinates)
     * @param samlingRate   sampling step = voxelsize / samplingrate
     * @param vhh           vector of volumes on which the calculation is done
     * @param fo            fetch option, how the voxels are interpolated
     * @param outData       plotData for the output
     * @return              true, if everthing went right
     */
    static bool calculate(const tgt::vec3 start,const tgt::vec3 end, float samplingRate, std::vector<const VolumeBase* >& vvh, fetchOption fo, PlotData* outData);
private:
    /**
     * Measures the voxel intensity at one point.
     * Measures the voxel intensity at one point for all volumes.
     *
     * @param samplePos     position of measurement
     * @param vhh           vector of volumes
     * @param intensity     vector of results
     * @param fo            fetch option
     * @return              true, if it worked all fine
     */
    static bool measureIntensity(const tgt::vec3 samplePos, std::vector<const VolumeBase* >& vvh, std::vector<plot_t>& intensity, fetchOption fo);
};

}   //namespace

#endif // VRN_LP_MEASURE_H
