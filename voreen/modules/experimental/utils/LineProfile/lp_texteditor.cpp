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

#include "lp_texteditor.h"

namespace voreen {
//******************************************************************************************
//            Text Output
//******************************************************************************************
void LP_TextEditor::setTextMaxValue(std::stringstream& strstr, float value){
    strstr << "Maximum Value: " << value << std::endl;
}

void LP_TextEditor::setTextThresholdArea(std::stringstream& strstr,const std::vector<std::pair<float,float> >& vec){
    //setting of the interval borders
    int i = 1;
    for(std::vector<std::pair<float,float> >::const_iterator it = vec.begin(); it != vec.end(); it++)
        strstr <<"Length of " << i++ <<  ". Interval: " << (it->second-it->first) << std::endl;
}

void LP_TextEditor::setTextThresholdValue(std::stringstream& strstr, float value){
    strstr << "Threshold Value:" << value << std::endl;
}

} // namespace voreen
