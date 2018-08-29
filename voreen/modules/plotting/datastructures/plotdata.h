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

#ifndef VRN_PLOTDATA_H
#define VRN_PLOTDATA_H

#include "plotbase.h"
#include "interval.h"

#include <vector>
#include <set>

namespace voreen {

class AggregationFunction;
class PlotPredicate;
class PlotCellValue;
class PlotCellImplicit;
class PlotRowValue;
class PlotRowImplicit;

/**
 * \brief   A PlotData stores rows which can be sorted by key columns. Insert and find calls should be fast.
 *
 * \note    As we are handling pointers in private member highlightedCells_, it is crucial to
 *          pay attention on implementing new member functions. Especially if you plan to
 *          modify rows_ make sure all pointers in this set remain valid or are removed!
 *
 **/
class VRN_CORE_API PlotData : public PlotBase {
    friend class PlotFunction;
public:
    /// possible types of the Plotdatamerge
    enum MergeType {
        NOSELECTION,
        IGNORECOLUMNLABELS,
        NEWCOLUMNS,
        NEWROWS
    };

    /**
     * \brief initializes keyColumnCount, dataColumnCount by given values and columnLabels with empty strings.
     *
     * \param keyColumnCount   count of key columns for this PlotBase
     * \param dataColumnCount  count of data columns for this PlotBase
     **/
    PlotData(int keyColumnCount, int dataColumnCount);

    /// copy constructor
    PlotData(const PlotData& rhs);

    /// default destructor
    virtual ~PlotData();

    /// assignment operator
    PlotData& operator=(const PlotData& rhs);

    /**
     * \brief  Selects the subset of stored rows matching the PlotPredicates in \a predicates
     *         and stores them in \a target
     *
     * \param predicates   Vector of pointers to PlotPredicates combined in a pair with an integer
     *                     indicating the column to apply the predicate to.
     **/
    void select(const std::vector<std::pair< int, PlotPredicate* > >& predicates, PlotData& target) const;

    /**
     * \brief  Selects all stored rows restricted to the columns given in \a columns vector.
     *
     * \param columns          Vector of integers indicating the columns to select.
     *                         Dimension is keyColumnCount + dataColumnCount.
     * \param keyColumnCount   Count of key columns for this PlotBase.
     * \param dataColumnCount  Count of data columns for this PlotBase.
     **/
    void select(const std::vector<int>& columns, int keyColumnCount, int dataColumnCount, PlotData& target) const;

    /**
     * \brief  Selects all stored rows restricted to the columns given in \a columns vector
     *         and matching the PlotPredicates in \a predicates.
     *
     * \param columns          Vector of integers indicating the columns to select.
     *                         Dimension is keyColumnCount + dataColumnCount.
     * \param keyColumnCount   Count of key columns for this PlotBase.
     * \param dataColumnCount  Count of data columns for this PlotBase.
     * \param predicates       Vector of pointers to PlotPredicates combined in a pair with an integer
     *                         indicating the column to apply the predicate to.
     **/
    void select(const std::vector<int>& columns, int keyColumnCount, int dataColumnCount, const std::vector<std::pair< int, PlotPredicate* > >& predicates, PlotData& target) const;

    /**
     * \brief  Selects all stored rows given in \a rows vector restricted to the columns given in \a columns vector.
     *
     * \param columns          Vector of integers indicating the columns to select.
     *                         Dimension is keyColumnCount + dataColumnCount.
     * \param rows             Vector of integers indicating the rows to select.
     * \param keyColumnCount   Count of key columns for this PlotBase.
     * \param dataColumnCount  Count of data columns for this PlotBase.
     **/
    void select(const std::vector<int>& columns, const std::vector<int>& rows, int keyColumnCount, int dataColumnCount, PlotData& target) const;

    /**
     * \brief Applies the given AggregationFunction \a function to all values in column \a column.
     *
     * \param column       destination column of \a function
     * \param function     pointer to AggregationFunction to apply
     *
     * \return     plot_t result of the given AggregationFunction
     **/
    plot_t aggregate(int column, const AggregationFunction* function) const;

    /**
     * \brief  Inserts a new PlotRowValue into PlotData with the values given in \a values. Assumes
     *         that the values are given in ascending column order (first key columns, then data columns)
     *
     * \param values    vector of values of the new row given given in ascending column order as
     *                  defined in this PlotData (first key columns, then data columns). If
     *                  values.size() < getColumnCount() missing columns will be filled with
     *                  empty PlotCellValues.
     *
     * \return      true upon success, false on failure (e.g. too many values/wrong indices)
     **/
    bool insert(const std::vector<plot_t>& values);

    /**
     * \brief  Inserts a new PlotRowValue into PlotData with the values given in the left part of the
     *         \a values pair. Each value is mapped to the column given by the right part of \a values.
     *
     * \param values   vector of values, each paired with the destination column index
     *
     * \return      true upon success, false on failure (e.g. too many values/wrong indices)
     **/
    bool insert(std::vector<std::pair<int, plot_t> >& values);

    /**
     * \brief  Inserts a new PlotRowValue into PlotData with the tags of the PlotCells given in \a tags. Assumes
     *         that the values are given in ascending column order (first key columns, then data columns)
     *
     * \param tags      vector of tags of the new row given given in ascending column order as
     *                  defined in this PlotData (first key columns, then data columns). If
     *                  values.size() < getColumnCount() missing columns will be filled with
     *                  empty PlotCellValues.
     *
     * \return      true upon success, false on failure (e.g. too many values/wrong indices)
     **/
    bool insert(const std::vector<std::string>& tags);

    /**
     * \brief  Inserts a new PlotRowValue into PlotData with the tags given in the left part of the
     *         \a values pair. Each tag is mapped to the column given by the right part of \a values.
     *
     * \param tags     vector of tags, each paired with the destination column index
     *
     * \return      true upon success, false on failure (e.g. too many values/wrong indices)
     **/
    bool insert(const std::vector<std::pair<int, std::string> >& tags);

    /**
     * \brief   Inserts a new PlotRowValue into PlotData containing the PlotCells given in \a cells.
     *
     * Assumes that the cells are given in ascending column order (first key columns, then data columns).
     *
     * \param cells     Vector of PlotCellValues of the new row given given in ascending column order as
     *                  defined in this PlotData (first key columns, then data columns). If
     *                  values.size() < getColumnCount() missing columns will be filled with
     *                  empty PlotCellValues.
     *
     * \return      true upon success, false on failure (e.g. too many values/wrong indices)
     *
     * \note    As \a cells is const this insert function is reasonably slower that the non-const version.
     *          Consider calling that instead.
     **/
    bool insert(const std::vector<PlotCellValue>& cells);

    /**
     * \brief  Inserts a new PlotRowValue into PlotData containing the PlotCells given in \a cells. Assumes
     *         that the cells are given in ascending column order (first key columns, then data columns)
     *
     * \param cells     Vector of PlotCellValues of the new row given given in ascending column order as
     *                  defined in this PlotData (first key columns, then data columns). If
     *                  values.size() < getColumnCount() missing columns will be filled with
     *                  empty PlotCellValues. Caution: As this parameter is a non-const reference the contained
     *                  cells might be changed if they do not fit to the column layout of this PlotData!
     *
     * \return      true upon success, false on failure (e.g. too many values/wrong indices)
     *
     * \note    This insert function is reasonably faster that the the version accepting only a const reference
     *          of cells. Thereforce \a cells might be changed if they do not fit to the column layout of this
     *          PlotData. If you cannot live with that, call the other version.
     **/
    bool insert(std::vector<PlotCellValue>& cells);

    /**
     * \brief   Inserts a new PlotRowValue into PlotData with the cells given in the left part of the
     *          \a values pair. Each cell is mapped to the column given by the right part of \a values.
     *
     * \param cells    vector of cells, each paired with the destination column index
     *
     * \return      true upon success, false on failure (e.g. too many values/wrong indices)
     **/
    bool insert(const std::vector<std::pair<int, PlotCellValue> >& cells);

    /**
     * \brief   Returns an iterator to the first row which does not compare less than \a values.
     *
     *          If the rows aren't yet sorted this will be done first. Comparison will be done
     *          lexicographically for the first <keyColumnCount_> elements. If \a values has less
     *          than <keyColumnCount_> elements it will be expanded with NULL-cells.
     *
     * \param   values  vector of PlotCellValues for comparison.
     *
     * \note    Caution: All iterators to rows may be invalidated due to sorting process!
     **/
    std::vector<PlotRowValue>::const_iterator lower_bound(const std::vector<PlotCellValue>& values) const;

    /**
     * \brief      removes all rows matching PlotPredicates in \a predicates.
     *
     * \param predicates       Array of PlotPredicates associated with an integer value indicating
     *                         the column to apply the predicate to.
     * \param predicateCount   count of elements in \a predicates array
     *
     * \return     the number of removed rows
     **/
    int remove(const std::vector<std::pair<int, PlotPredicate*> >& predicates);

    /**
     * \brief   Copies the rows with rearranged columns into target PlotData.
     *          The key columns will be the columns with the indices as defined in \a newKeyColumns, the data
     *          columns the columns with the indices as defined in \a newDataColumns. \a target will be
     *          be reset before inserting the rows: column count will be newKeyColumns.size() + newDataColumns.size().
     *
     * \param   newKeyColumns   indices of the new key columns, must be less than this->getColumnCount()
     * \param   newDataColumns  indices of the new data columns, must be less than this->getColumnCount()
     * \param   target          target PlotData where the rearranged PlotRows will be stored
     *
     * \return      true upon success, false on failure (e.g. index out of bounds)
     **/
    bool rearrangeColumns(std::vector<int>& newKeyColumns, std::vector<int>& newDataColumns, PlotData& target) const;

    /**
     * \brief   executes SQL-like GROUP BY operation on this PlotBase.
     *          The PlotRows will be grouped by the column with index \a groupColumn, that means all PlotRows containing
     *          the same value/tag in column with index \a groupColumn will be combined in one new PlotRow and inserted
     *          into target. All other columns will be rejected if not referred as second part in \a functions: The values
     *          in the columns with indices referred as second part in \a functions will be aggregated with the according
     *          AggregationFunction and inserted into the corresponding grouped PlotRow.
     *          \a target will be reset before inserting the rows: column count will be (1+functions.size()).
     *
     * \param   groupColumn     column index to apply grouping to
     * \param   functions       set of pairs of an AggregationFunction and a column index which specifies where to apply
     *                          the AggregationFunction to
     * \param   target          target PlotData where the grouped rows will be stored
     *
     * \return  true upon success, false on failure (e.g. index out of bounds)
     **/
    bool groupBy(int groupColumn, const std::vector<std::pair<int, AggregationFunction*> >& functions, PlotData& target) const;

    /**
     * \brief   merges this PlotData with \a otherPlotData into \a target
     *          Merges all rows of *this and all rows of \a otherPlotData into \a target according to the column mapping
     *          defined in \a keyColumns and \a dataColumns: for each pair in these lists the target PlotData will have
     *          a column: the first part defines the source column of this, the second part the source column of
     *          \a otherPlotData.
     *          \a target will be reset before inserting the rows: column count will be keyColumns.size() + dataColumns.size().
     *
     * \param   otherPlotData   the second source PlotData for the merging process
     * \param   keyColumns      defines the source columns for the new key columns: first indices the column of *this,
     *                          second indices the column of \a otherPlotData
     * \param   dataColumns     defines the source columns for the new data columns: first indices the column of *this,
     *                          second indices the column of \a otherPlotData
     * \param   target          target PlotData where the merged rows will be stored
     *
     * \return  true upon success, false on failure (e.g. index out of bounds)
     **/
    bool mergeWith(const PlotData& otherPlotData,
                    std::vector<std::pair<int, int> >& keyColumns,
                    std::vector<std::pair<int, int> >& dataColumns,
                    PlotData& target) const;
    /**
     * \brief   merges this PlotData with \a otherPlotData into \a target
     *          \a target will be reset before inserting the rows
     *
     * \param   otherPlotData   the second source PlotData for the merging process
     * \param   target          target PlotData where the merged rows will be stored
     *
     * \return  true upon success, false on failure (e.g. index out of bounds)
     **/
    bool mergeWith(const PlotData& otherPlotData, PlotData& target, MergeType mergeType = NOSELECTION) const;

    /**
     * \brief   inserts an implicit row into this PlotData
     *          Implicit rows do not contain data, their values are automatically calculated by the
     *          AggregationFunctions defined in the first parts in \a functions which are applied to the columns
     *          defined in the second parts in \a functions.
     *
     * \param   functions   list of pairs defining the AggregationFunctions and the columns on which to apply
     *                      the AggregationFunctions to
     *
     * \return  true upon success.
     **/
    bool insertImplicit(const std::vector<std::pair<int, AggregationFunction*> >& functions);

    /**
     * \brief   inserts an implicit row into this PlotData
     *          Implicit rows do not contain data, their values are automatically calculated by the
     *          AggregationFunctions defined in the first parts in \a functions which are applied to the columns
     *          defined in the second parts in \a functions.
     *
     * \param   implicitrow   list of PlotCellImplicit
     *
     * \return  true upon success.
     **/
    bool insertImplicit(const std::vector<PlotCellImplicit>& implicitrow);

    /**
     * \brief   returns the interval of the values in column \a column of this PlotData
     *          Calculates the min/max values of all PlotCells in the column with index \a column and
     *          returns them as Interval<plot_t>. The interval of a string column is the number of rows.
     *
     * \param   column  index of the column of which to calculate the value interval
     *
     * \return  interval with minimum and maximum values
     **/
    Interval<plot_t> getInterval(int column) const;

    /**
     * \brief   returns the interval from 0 to the maximal sum of absolute values in one row using the columns \a column
     *          of this PlotData. String columns will be ignored.
     *
     * \param   column  indices of the columns of which to calculate the value interval
     *
     * \return  interval from 0 to the maximal sum
     **/
    Interval<plot_t> getSumInterval(const std::vector<int>& column) const;

    /// Returns a const iterator to the first PlotRowValue of this PlotData.
    std::vector<PlotRowValue>::const_iterator getRowsBegin() const;
    /// Returns a const iterator to the first PlotRowImplicit of this PlotData.
    std::vector<PlotRowImplicit>::const_iterator getImplicitRowsBegin() const;

    /// Returns a const iterator to the end of the PlotRowValues of this PlotData.
    /// Remember:   Do not try to dereference this iterator as it points to the element
    ///             behind the last element!
    std::vector<PlotRowValue>::const_iterator getRowsEnd() const;
    /// Returns a const iterator to the end of the PlotRowImplicits of this PlotData.
    /// Remember:   Do not try to dereference this iterator as it points to the element
    ///             behind the last element!
    std::vector<PlotRowImplicit>::const_iterator getImplicitRowsEnd() const;

    /// Returns the PlotRowValue in row \a row.
    const PlotRowValue& getRow(int row) const;
    /// Returns the PlotRowImplicit in row \a row.
    const PlotRowImplicit& getImplicitRow(int row) const;

    /// Returns the number of PlotRowValues in this PlotData.
    int getRowsCount() const;
    /// Returns the number of PlotRowImplicits in this PlotData.
    int getImplicitRowsCount() const;

    /// Returns whether this PlotData has PlotRowValues or not.
    bool rowsEmpty() const;
    /// Returns whether this PlotData has PlotRowImplicits or not.
    bool implicitRowsEmpty() const;
    /// Checks if this PlotCell is highlighted
    bool isHighlighted(const tgt::ivec2& cellPosition) const;

    /**
     * \brief   Changes the highlighted-flag state of all PlotCellValues named in \a cellPosition.
     *
     * \note    If \a cellPosition is not valid, no cell will be toggled.
     *
     * \param   cellPosition  Table position (row, column) of the cell to highlight - as in
     *                        PlotPickingManager the value -1 encodes whole row/column.
     * \param   additive      If true, no other cells will be affected, if false the highlighted-
     *                        state of all other cells will be set to false
     **/
    void toggleHighlight(const tgt::ivec2& cellPosition, bool additive);

    /**
     * \brief   The highlighted-flag of all other PlotCellValues will be set to false.
     **/
    void clearHighlights();

    /**
     * \brief   Sets the highlighted-flag of all PlotCellValues named in \a cellPosition to \a value.
     *
     * \note    If \a cellPosition is not valid, no cell will be highlighted.
     *
     * \param   cellPosition Table position (row, column) of the cell to highlight - as in
     *                       PlotPickingManager the value -1 encodes whole row/column.
     * \param   value        highlight state
     * \param   additive     If true, no other cells will be affected, if false the highlighted-
     *                       state of all other cells will be set to false
     **/
    void setHighlight(const tgt::ivec2& cellPosition, bool value, bool additive);

    /**
     * \brief Clears all data and resets key- and data column count.
     *
     * \param   keyColumnCount  the new count of key columns
     * \param   dataColumnCount the new count of data columns
     **/
    void reset(int keyColumnCount, int dataColumnCount);

    /**
     * Ensures that rows_ is sorted lexicographically by key columns.
     *
     * \note    All iterators to PlotRowValues will be invalidated!
     **/
    void sortRows() const;

    /// Returns whether the data is sorted.
    bool sorted() const;

    /**
     * \brief   Returns whether the given column is an index column in \a pData.
     *
     * Checks for column for being 0 (first column), column type NUMBER and column header "Index".
     *
     * \param   pData   PlotData containing the column
     * \param   column  column index
     **/
    static bool isIndexColumn(const PlotData& pData, int column);

private:
    /// updateIntervals
    void updateIntervals(const PlotRowValue& row);
    /// delete implicit rows and their values
    void deleteImplicitRows();

    /**
     * \brief   Removes all pointers to cells of PlotRow \a row from highlightedCells_.
     *
     * \note    Make sure to call this function before removing \a row from rows_. If not
     *          sooner or later everything here will explode!
     **/
    void removePointersToCellsOfRow(PlotRowValue& row);

    /**
     * \brief   all value rows of this PlotData
     *
     * \note    As we might hold pointers to cells in entries of this vector it is crucial
     *          to use only member functions to modify this vector.
     **/
    std::vector<PlotRowValue> rows_;

    std::vector<PlotRowImplicit> implicitRows_; ///< all implicit rows of this PlotData
    std::vector< Interval<plot_t> > intervals_; ///< cached intervals of each column_

    /**
     * \brief   List of pointers to all highlighted cells.
     *
     * This set is needed to be able to clear the highlight-flags of all PlotCellValues
     * without iterating over all existing cells.
     *
     * \note    As we are handling pointers here, it is crucial to maintain this set via
     *          the corresponding accessors to ensure nothing explodes! Also all changes
     *          in rows_, especially removing entries must done using member functions
     *          or with extreme caution!
     *          On implementing new member functions make sure all pointers in this list
     *          remain valid or are removed.
     **/
    std::set<PlotCellValue*> highlightedCells_;

    /// flag whether rows_ is sorted lexicographically by key columns or not
    mutable bool sorted_;

};

} // namespace voreen

#endif // VRN_PLOTDATA_H
