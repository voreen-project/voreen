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

#ifndef VRN_LP_PLOTMODIFIER_H
#define VRN_LP_PLOTMODIFIER_H

//includes for plotdata modifying
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"
//vector
#include <vector>

namespace voreen {

/**
 * Utility class beeing used in lineprofile.h
 * This class encapsulate the functions modifying plotdata.
 */
class LP_PlotModifier {
public:
    /**
     * Finds the maximum in a given column
     * @param pData     the plotdata we are searching in
     * @param column    the column we are searching in
     * @return          the maximum value
     */
    static float findMax(const PlotData* pData, int column);
    /**
     * Finds the minimum in a given column
     * @param pData     the plotdata we are searching in
     * @param column    the column we are searching in
     * @return          the minimum value
     */
    static float findMin(const PlotData* pData, int column);
    /**
     * Finds the value of x coordinate.
     * This function returns the x value (the value in column 0) at percent% from [begin,end] of the column.
     * For instance (1,4,5,9) and 25% will return 1+(9-1)/4 = 3.
     * @param pData     the plotdata we are searching in
     * @param percent   the percentage value (can be 0.25 or 25)
     * @return          value at the the point specified by percent
     */
    static float findXValue(const PlotData* pData, float percent);
    /**
     * Finds the value of y coordinate.
     * This function returns the y value at percent% from [begin,end] of the column "column".
     * For instance (1,4,5,9) and 25% will return 1+(9-1)/4 = 3.
     * @param pData     the plotdata we are searching in
     * @param percent   the percentage value (can be 0.25 or 25)
     * @param column    the column we are interested in
     * @return          value at the the point specified by percent
     */
    static float findYValue(const PlotData* pData, float percent, int column);
    /**
     * Adds a column at the end of a given plotdata.
     * This function adds a column with values of vec at the end of a given plotdata with the given label
     * @param pData     plotdata we like to modify
     * @param vec       vector of values to add (must have same size as pData->getRowsCount())
     * @param label     the label of teh new column
     */
    static void addColumn(PlotData* pData, const std::vector<float>& vec,const std::string& label);
    /**
     * Adds a column at the end of a given plotdata.
     * This function adds a column with a constant value at the end of a given plotdata with the given label
     * @param pData     plotdata we like to modify
     * @param value     value of the new column
     * @param label     the label of the new column
     */
    static void addFloatColumn(PlotData* pData, float value, const std::string& label);
    /**
     * Stores a column into a vector.
     * This function stores a given column into a vector
     * @param pData     plotdata we like to modify
     * @param vec       vector the column is stored in
     * @param column    the column index we want to store
     */
    static void getColumn(const PlotData* pData, std::vector<float>& vec, int column);
    /**
     * Removes a column from a given plotdata.
     * @param pData     plotdata we like to modify
     * @param column    the column we want to remove
     */
    static void removeColumn(PlotData* pData, int column);
    /**
     * Merges two plotdatas into one.
     * This function merges two plotdatas into one. The new one is pOutput or if not specified pInput1
     * @param pInput1   first plotdata (and result plotdata if pOutput is 0)
     * @param pInput2   second plotdata
     * @param pOutput   the merged plotdata (can be null)
     */
    static void mergePlotData(PlotData* pInput1, const PlotData* pInput2, PlotData* pOutput = 0);
     /**
     * Adds a vector of rows to a given plotdata.
     * The row count of the plotData and the vector size have to be the same.
     * @param pData     plotdata, which will be modified
     * @param vec       vector of rows (will be modified too)
     */
    static void addRowVector(PlotData* pData, std::vector< std::vector<PlotCellValue> > & vec);
};

}   //namespace

#endif // VRN_LP_PLOTMODIFIER_H
