/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTROW_H
#define VRN_PLOTROW_H

#include "plotcell.h"
#include "plotbase.h"
#include "interval.h"
#include "aggregationfunction.h"

#include <vector>

namespace voreen {

class PlotData;

/**
* Single row of a PlotBase table storing specific values.
*
* Held PlotCellValues store either plot_t or string values.
*/
class VRN_CORE_API PlotRowValue {

// PlotData is our friend for two reasons:
//  - We want to grant PlotData access to our private constructor.
//  - Unfortunately PlotData needs some non-const pointers to cells in here. As this is
//    the one and only special case, where we want to give non-const access to our
//    internals, we declare PlotData as our friend.
friend class PlotData;

public:
    /// Destructor
    ~PlotRowValue();


    /// returns all stored PlotCells
    const std::vector<PlotCellValue>& getCells() const;

    /**
     * Returns the PlotCell at column \a column.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    const PlotCellValue& getCellAt(int column) const;

    /**
     * Returns the value in PlotCell at column \a column.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    plot_t getValueAt(int column) const;

    /**
     * Returns a string-pointer to the tag of PlotCell at column \a column.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    std::string getTagAt(int column) const;

    /**
     * Sets hightlighted-flag of cell in column \a column to \a value.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    void setHighlighted(int column, bool value);

    /**
    * \property bool operator<(PlotRowValue& rhs) const
    *
    * \brief  compares two PlotRowValues lexicographically by their key values
    *
    * \param rhs  second PlotRowValue to compare
    *
    * \return     true if lhs < rhs (lexicographic order)
    **/
    bool operator<(const PlotRowValue& rhs) const;


private:
    /**
    * Constructor for a PlotRow containing the given PlotCellsValues.
    *
    * This constructor is private because only PlotData should call this.
    *
    * \param  parent      parent PlotData - must be valid
    * \param  cells       array of const PlotCell pointers assuming dimensions are
    *                     parent->keyColumnCount + parent->dataColumnCount
    *
    * \note    PlotRowValue does not take ownership of \a parent and as long as this object lives
    *          it is absolutely necessary that \a parent exists.
    *          The easiest way to ensure that is not to call this constructor directly and let the
    *          parent PlotData handle its rows.
    */
    PlotRowValue(const PlotData* parent, const std::vector<PlotCellValue>& cells);

    /**
    * Constructor for a PlotRow containing the given PlotCellsValues.
    *
    * This constructor is private because only PlotData should call this.
    *
    * \param  parent      parent PlotData - must be valid
    * \param  cells       array of const PlotCell pointers assuming dimensions are
    *                     parent->keyColumnCount + parent->dataColumnCount
    *
    * \note    PlotRowValue does not take ownership of \a parent and as long as this object lives
    *          it is absolutely necessary that \a parent exists.
    *          The easiest way to ensure that is not to call this constructor directly and let the
    *          parent PlotData handle its rows.
    */
    PlotRowValue(const PlotData* parent, std::vector<PlotCellValue>& cells);


    /// held PlotCells, dimensions will be (parent.getKeyColumnCount() + parent.getDataColumnCount())
    std::vector<PlotCellValue> cells_;

    const PlotData* parent_;       ///< parent PlotData

};

/**
 * Single row of a PlotBase table storing AggregationFunctions.
 *
 * Data values are given given implicit and compute by the aggregation function on calling getValueAt().
 */
class VRN_CORE_API PlotRowImplicit {
// grant PlotData access to private constructor
friend class PlotData;

public:
    /// Destructor
    ~PlotRowImplicit();


    /// returns all stored PlotCells
    const std::vector<PlotCellImplicit>& getCells() const;

    /**
     * Returns the PlotCell at column \a column.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    const PlotCellImplicit& getCellAt(int column) const;

    /**
     * returns the value in PlotCell at column \a column.
     *
     * \note    Be sure \a column is within bounds (>= 0 && < getColumnCount()) or this function might crash!
     **/
    plot_t getValueAt(int column) const;

private:
    /**
     * Constructor for a PlotRow containing PlotCellImplicits with the given AggregationFunctions.
     *
     *
     * This constructor is private because only PlotData should call this.
     *
     * \param   parent      parent PlotData - must be valid
     * \param   functions   functions AggregationFunction
     *
     * \note    This Constructor make deep copies of the given AggregationFunctions.
     *          PlotRowValue does not take ownership of \a parent and as long as this object lives
     *          it is absolutely necessary that \a parent exists.
     *          The easiest way to ensure that is not to call this constructor directly and let the
     *          parent PlotData handle its rows.
     */
    PlotRowImplicit(const PlotData* parent, const std::vector<AggregationFunction*>& functions);

    /// held PlotCells, dimensions will be (parent.getKeyColumnCount() + parent.getDataColumnCount())
    std::vector<PlotCellImplicit> cells_;

    const PlotData* parent_;       ///< parent PlotData

};

} // namespace voreen

#endif // VRN_PLOTROW_H
