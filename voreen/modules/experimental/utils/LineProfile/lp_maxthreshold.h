/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_LP_MAXTHRESHOLD_H
#define VRN_LP_MAXTHRESHOLD_H

#include "modules/plotting/datastructures/plotdata.h"

namespace voreen {

/**
 * Utility class beeing used in lineprofile.h
 * This class encapsulate the functions for threshold calculation.
 */
class LP_MaxThreshold {
public:
    /**
     * Set the maxValue value and position value of the LineProfile for a given column
     * @param maxValue  the maxValue member of the LineProfile
     * @param position  the position of the maximum
     * @param pData     the plotdata
     * @param column    the maximum of this column
     */
    static void setMaxValue(float& maxValue, float& position, const PlotData* pData, int column);
    /**
     * Gets all Threshold intervals.
     * This function gets all intervals (x-axes) which y-value is above a given value in a given column.
     * @param threshold     the threshold value y must be greater
     * @param vec           result vector if intervals
     * @param pData         plotdata we are looking in
     * @param column        column we are looking in
     */
    static void setThresholdArea(float threshold, std::vector<std::pair<float,float> >& vec, const PlotData* pData, int column);

};

}   //namespace

#endif // VRN_LP_MAXTHRESHOLD_H
