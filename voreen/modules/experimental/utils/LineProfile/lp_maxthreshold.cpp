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

#include "lp_maxthreshold.h"
//for navigation through the plotdata
#include "modules/plotting/datastructures/plotrow.h"


namespace voreen {

//******************************************************************************************
//                    Max Line
//******************************************************************************************
void LP_MaxThreshold::setMaxValue(float& maxValue, float& position, const PlotData* pData, int column){
    if(!pData || (pData->getColumnCount() <= column) || (pData->getRowsCount() == 0))
        return;
    maxValue = static_cast<float>(pData->getRow(0).getValueAt(column));
    position = static_cast<float>(pData->getRow(0).getValueAt(0));
    for (int i = 1; i < pData->getRowsCount(); i++){
        if(pData->getRow(i).getValueAt(column) > maxValue){
            maxValue = (float)pData->getRow(i).getValueAt(column);
            position = (float)pData->getRow(i).getValueAt(0);
        }
    }
}


void LP_MaxThreshold::setThresholdArea(float threshold, std::vector<std::pair<float,float> >& vec, const PlotData* pData, int column) {
    if(!pData || (column >= pData->getColumnCount()))
        return;
    //set variables
    bool inInter = false; std::pair<float,float> interval;
    //clear threshold vector
    vec.clear();
    //for each row
    for(int i = 0; i < pData->getRowsCount(); i++){
        if(inInter){
            if(pData->getRow(i).getValueAt(column) >= threshold) {
                interval.second = static_cast<float>(pData->getRow(i).getValueAt(0));
            } else {
                inInter = false;
                vec.push_back(interval);
            }
        } else {//not in interval
            if(pData->getRow(i).getValueAt(column) >= threshold) {
                interval.first = static_cast<float>(pData->getRow(i).getValueAt(0));
                inInter = true;
                interval.second = static_cast<float>(pData->getRow(i).getValueAt(0));
            }
        }
    }
    if(inInter)
        vec.push_back(interval);
}

} // namespace voreen
