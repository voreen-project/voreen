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

#ifndef VRN_LP_TEXTEDITOR_H
#define VRN_LP_TEXTEDITOR_H

//stream
#include <sstream>
//vector
#include <vector>

namespace voreen {

/**
 * Utility class beeing used in lineprofile.h
 * This class encapsulate the text output functions used in lineprofile::setTextOutport.
 */
class LP_TextEditor {
public:
    /**
     * Function to to give out the max value.
     *
     * @param strstr    stream for textoutput
     * @param value     max value for the output
     */
    static void setTextMaxValue(std::stringstream& strstr, float value);

    /**
     * Function to to give out the length of the threshold itervals.
     *
     * @param strstr    stream for textoutput
     * @param vec       vector of all threshold intervals
     */
    static void setTextThresholdArea(std::stringstream& strstr, const std::vector<std::pair<float,float> >& vec);

    /**
     * Function to to give out the threshold value.
     *
     * @param strstr    stream for textoutput
     * @param value     threshold value for the output
     */
    static void setTextThresholdValue(std::stringstream& strstr, float value);
};

}   //namespace

#endif // VRN_LP_TEXTEDITOR_H
