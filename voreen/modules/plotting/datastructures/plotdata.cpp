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

#include "plotdata.h"
#include "plotpredicate.h"
#include "plotrow.h"
#include "plotcell.h"

#include "tgt/assert.h"
#include "tgt/logmanager.h"

#include <map>
#include <vector>
#include <algorithm>
#include <limits>
#include <sstream>

namespace voreen {

PlotData::PlotData(int keyColumnCount, int dataColumnCount)
    : PlotBase(keyColumnCount, dataColumnCount)
{
    intervals_.resize(keyColumnCount_ + dataColumnCount_);
}

PlotData::PlotData(const PlotData& rhs)
    : PlotBase(rhs)
    , rows_(rhs.rows_)
    , intervals_(rhs.intervals_)
    , sorted_(rhs.sorted_)
{
    for (std::vector<PlotRowValue>::iterator it = rows_.begin(); it < rows_.end(); ++it) {
        it->parent_ = this;
        for (std::vector<PlotCellValue>::iterator cit = it->cells_.begin(); cit != it->cells_.end(); ++cit) {
            if (cit->isHighlighted())
                highlightedCells_.insert(&(*cit));
        }
    }

    std::vector<PlotRowImplicit>::const_iterator imit;
    for (imit = rhs.implicitRows_.begin(); imit < rhs.implicitRows_.end(); ++imit) {
        insertImplicit(imit->getCells());
    }
}

PlotData::~PlotData() {
    deleteImplicitRows();
}

PlotData& PlotData::operator=(const PlotData& rhs) {
    // prevent self assignment
    if (this == &rhs)
        return *this;

    PlotBase::operator=(rhs);

    // we do not use the copy-and-swap pattern here for performance reasons (avoid walking through all cells twice)
    // if an exception raises, this object will be empty but valid
    try {
        rows_ = rhs.rows_;
        intervals_ = rhs.intervals_;
        sorted_ = rhs.sorted_;
        highlightedCells_.clear();

        for (std::vector<PlotRowValue>::iterator it = rows_.begin(); it < rows_.end(); ++it) {
            it->parent_ = this;
            for (std::vector<PlotCellValue>::iterator cit = it->cells_.begin(); cit != it->cells_.end(); ++cit) {
                if (cit->isHighlighted())
                    highlightedCells_.insert(&(*cit));
            }
        }

        std::vector<PlotRowImplicit>::const_iterator imit;
        for (imit = rhs.implicitRows_.begin(); imit < rhs.implicitRows_.end(); ++imit) {
            insertImplicit(imit->getCells());
        }
    }
    catch (std::bad_alloc&) {
        // something bad happened, we can't be sure about the current state
        // clear everything to have at least a valid state.
        rows_.clear();
        implicitRows_.clear();
        intervals_.clear();
        highlightedCells_.clear();
        LERRORC("PlotData::operator=()", "bad_alloc occured, object won't contain any data!");
        return *this;
    }
    catch (...) {
        // something bad happened, we can't be sure about the current state
        // clear everything to have at least a valid state.
        rows_.clear();
        implicitRows_.clear();
        intervals_.clear();
        highlightedCells_.clear();
        LERRORC("PlotData::operator=()", "unknown exception occured, object won't contain any data!");
        return *this;
    }

    return *this;
}

void PlotData::select(const std::vector< std::pair< int, PlotPredicate*> >& predicates, PlotData& target) const {
    target.reset(keyColumnCount_, dataColumnCount_);

    std::vector<PlotRowValue>::const_iterator rit;
    bool matches;
    for (rit = rows_.begin(); rit < rows_.end(); ++rit) {
        matches = true;
        for (size_t i = 0; i < predicates.size(); ++i) {
            if (!predicates.at(i).second->check(rit->getCellAt(predicates.at(i).first))) {
                matches = false;
                break;
            }
        }
        if (matches) {
            target.insert(rit->getCells());
        }
    }
    for (int i = 0; i < getColumnCount(); ++i) {
        target.setColumnLabel(i,getColumnLabel(i));
    }
    std::vector<PlotRowImplicit>::const_iterator imit;
    for (imit = implicitRows_.begin(); imit < implicitRows_.end(); ++imit) {
        target.insertImplicit(imit->getCells());
    }
}

void PlotData::select(const std::vector<int>& columns, int keyColumnCount, int dataColumnCount, PlotData& target) const {
    target.reset(keyColumnCount, dataColumnCount);
    int columnCount = 0;
    if (columns.size() != 0) {
        columnCount = keyColumnCount + dataColumnCount;

        std::vector<PlotRowValue>::const_iterator it;
        int i;
        for (it = rows_.begin(); it < rows_.end(); ++it) {
            std::vector<PlotCellValue> cellsToInsert;

            for (i=0; i< columnCount; ++i) {
                cellsToInsert.push_back(it->getCellAt(columns[i]));
            }

            target.insert(cellsToInsert);
        }
        for (i = 0; i < columnCount; ++i) {
            target.setColumnLabel(i,getColumnLabel(columns[i]));
        }
        std::vector<PlotRowImplicit>::const_iterator imit;
        for (imit = implicitRows_.begin(); imit < implicitRows_.end(); ++imit) {
            std::vector<PlotCellImplicit> implicitcellsToInsert;
            for (i=0; i< columnCount; ++i) {
                implicitcellsToInsert.push_back(imit->getCellAt(columns[i]));
            }
            target.insertImplicit(implicitcellsToInsert);
        }
    }
}

void PlotData::select(const std::vector< int >& columns, int keyColumnCount, int dataColumnCount,
                      const std::vector< std::pair< int, voreen::PlotPredicate* > >& predicates, voreen::PlotData& target) const {
    target.reset(keyColumnCount, dataColumnCount);
    int columnCount = 0;
    if (columns.size() != 0) {
        columnCount = keyColumnCount + dataColumnCount;

        std::vector<PlotRowValue>::const_iterator it;
        std::vector<std::pair<int, PlotPredicate*> >::const_iterator pit;
        bool matches;
        int i;
        for (it = rows_.begin(); it < rows_.end(); ++it) {
            std::vector<PlotCellValue> cellsToInsert;
            matches = true;
            for (pit = predicates.begin(); pit < predicates.end(); ++pit) {
                if (! pit->second->check(it->getCellAt(pit->first))) {
                    matches = false;
                    break;
                }
            }
            if (matches) {
                for (i=0; i< columnCount; ++i) {
                    cellsToInsert.push_back(it->getCellAt(columns[i]));
                }
                target.insert(cellsToInsert);
            }
        }
        for (i = 0; i < columnCount; ++i) {
            target.setColumnLabel(i,getColumnLabel(columns[i]));
        }
        std::vector<PlotRowImplicit>::const_iterator imit;
        for (imit = implicitRows_.begin(); imit < implicitRows_.end(); ++imit) {
            std::vector<PlotCellImplicit> implicitcellsToInsert;
            for (i=0; i< columnCount; ++i) {
                implicitcellsToInsert.push_back(imit->getCellAt(columns[i]));
            }
            target.insertImplicit(implicitcellsToInsert);
        }
    }
}

void PlotData::select(const std::vector<int>& columns, const std::vector<int>& rows, int keyColumnCount,
                      int dataColumnCount, PlotData& target) const {
    target.reset(keyColumnCount, dataColumnCount);
    int columnCount = 0;
    if (columns.size() != 0) {
        columnCount = keyColumnCount + dataColumnCount;

        std::vector<PlotRowValue>::const_iterator it;
        int i;
        int zaehler = 0;
        bool match;
        for (it = rows_.begin(); it < rows_.end(); ++it) {
            match = false;
            for (size_t j = 0; j < rows.size(); ++j) {
                if (zaehler == rows.at(j)) {
                    match = true;
                    break;
                }
            }
            ++zaehler;
            if (!match)
                continue;
            std::vector<PlotCellValue> cellsToInsert;

            for (i=0; i< columnCount; ++i) {
                cellsToInsert.push_back(it->getCellAt(columns[i]));
            }

            target.insert(cellsToInsert);
        }
        for (i = 0; i < columnCount; ++i) {
            target.setColumnLabel(i,getColumnLabel(columns[i]));
        }
        std::vector<PlotRowImplicit>::const_iterator imit;
        for (imit = implicitRows_.begin(); imit < implicitRows_.end(); ++imit) {
            std::vector<PlotCellImplicit> implicitcellsToInsert;
            for (i=0; i< columnCount; ++i) {
                implicitcellsToInsert.push_back(imit->getCellAt(columns[i]));
            }
            target.insertImplicit(implicitcellsToInsert);
        }
    }
}

plot_t PlotData::aggregate(int column, const AggregationFunction* function) const {
    std::vector<plot_t> values;
    std::vector<PlotRowValue>::const_iterator it;
    for (it = rows_.begin(); it < rows_.end(); ++it) {
        values.push_back(it->getValueAt(column));
    }

    plot_t toReturn = function->evaluate(values);
    return toReturn;
}

bool PlotData::insert(const std::vector<plot_t>& values) {
    if (static_cast<int>(values.size()) <= getColumnCount()) {
        std::vector<plot_t>::const_iterator it;
        std::vector<PlotCellValue> cellsToInsert;
        int i = 0;
        for (it = values.begin(); it < values.end(); ++it) {
            if (getColumnType(i) == NUMBER) {
                if (*it == *it)
                    cellsToInsert.push_back(PlotCellValue(*it));
                else
                    cellsToInsert.push_back(PlotCellValue());
            }
            else if (getColumnType(i) ==  EMPTY) {
                if (*it == *it) {
                    cellsToInsert.push_back(PlotCellValue(*it));
                    setColumnType(i,NUMBER);
                }
                else
                    cellsToInsert.push_back(PlotCellValue());
            }
            else {
                if (*it == *it) {
                    std::stringstream ss;
                    ss << *it;
                    cellsToInsert.push_back(PlotCellValue(ss.str()));
                }
                else
                    cellsToInsert.push_back(PlotCellValue());
            }
            ++i;
        }

        sorted_ = false;
        rows_.push_back(PlotRowValue(this, cellsToInsert));
        updateIntervals(rows_.back());
        return true;
    }
    return false;
}

bool PlotData::insert(std::vector<std::pair<int, plot_t> >& values) {
    if (static_cast<int>(values.size()) <= getColumnCount()) {
        std::vector<std::pair<int,plot_t> >::const_iterator it;
        std::vector<PlotCellValue> cellsToInsert;

        for (int i=0; i< getColumnCount(); ++i)
            cellsToInsert.push_back(PlotCellValue());

        for (it = values.begin(); it < values.end(); ++it) {
            if (getColumnType(it->first) == NUMBER) {
                if (it->second == it->second)
                    cellsToInsert[it->first] = PlotCellValue(it->second);
                else
                    cellsToInsert[it->first] = PlotCellValue();
            }
            else if (getColumnType(it->first) ==  EMPTY) {
                if (it->second == it->second) {
                    cellsToInsert[it->first] = PlotCellValue(it->second);
                    setColumnType(it->first,NUMBER);
                }
                else
                    cellsToInsert[it->first] = PlotCellValue();
            }
            else {
                if (it->second == it->second) {
                    std::stringstream ss;
                    ss << it->second;
                    cellsToInsert[it->first] = PlotCellValue(ss.str());
                }
                else
                    cellsToInsert[it->first] = PlotCellValue();

            }
        }
        sorted_ = false;
        rows_.push_back(PlotRowValue(this, cellsToInsert));
        updateIntervals(rows_.back());
        return true;
    }
    return false;
}

bool PlotData::insert(const std::vector<std::string>& tags) {
    if (static_cast<int>(tags.size()) <= getColumnCount()) {
        std::vector<std::string>::const_iterator it;
        std::vector<PlotCellValue> cellsToInsert;
        int i =0;
        for (it = tags.begin(); it < tags.end(); ++it) {
            if (getColumnType(i) == STRING) {
                cellsToInsert.push_back(PlotCellValue(*it));
            }
            else if (getColumnType(i) ==  EMPTY) {
                cellsToInsert.push_back(PlotCellValue(*it));
                setColumnType(i,STRING);
            }
            else {
                std::stringstream ss((*it));
                double num;
                if (ss >> num) {
                    cellsToInsert.push_back(PlotCellValue(num));
                }
                else {
                    cellsToInsert.push_back(PlotCellValue());
                }
            }
            ++i;
        }

        sorted_ = false;
        rows_.push_back(PlotRowValue(this, cellsToInsert));
        updateIntervals(rows_.back());
        return true;
    }
    return false;
}

bool PlotData::insert(const std::vector<std::pair<int, std::string> >& tags) {
    if (static_cast<int>(tags.size()) <= getColumnCount()) {
        std::vector<std::pair<int, std::string> >::const_iterator it;
        std::vector<PlotCellValue> cellsToInsert;

        for (int i=0; i< getColumnCount(); ++i)
            cellsToInsert.push_back(PlotCellValue());

        for (it = tags.begin(); it < tags.end(); ++it) {
            if (getColumnType(it->first) == STRING) {
                cellsToInsert[it->first] = PlotCellValue(it->second);
            }
            else if (getColumnType(it->first) ==  EMPTY) {
                cellsToInsert[it->first] = PlotCellValue(it->second);
                setColumnType(it->first,STRING);
            }
            else {
                std::stringstream ss(it->second);
                double num;
                if (ss >> num) {
                    cellsToInsert[it->first] = PlotCellValue(num);
                }
                else {
                    cellsToInsert[it->first] = PlotCellValue();
                }
            }
        }

        sorted_ = false;
        rows_.push_back(PlotRowValue(this, cellsToInsert));
        updateIntervals(rows_.back());
        return true;
    }
    return false;
}

bool PlotData::insert(const std::vector< PlotCellValue >& cells) {
    if (static_cast<int>(cells.size()) <= getColumnCount()) {
        std::vector< PlotCellValue > newcells(0);
        for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
            if ((getColumnType(i) == NUMBER && cells[i].isValue()) ||
                (getColumnType(i) == STRING && cells[i].isTag())) {
                newcells.push_back(cells[i]);
            }
            else if (getColumnType(i) == EMPTY) {
                newcells.push_back(cells[i]);
                if (cells[i].isTag()) {
                    setColumnType(i,STRING);
                }
                else if (cells[i].isValue() ) {
                    setColumnType(i,NUMBER);
                }
            }
            else if (getColumnType(i) == NUMBER && cells[i].isTag()) {
                std::stringstream ss(cells[i].getTag());
                double num;
                if (ss >> num) {
                    newcells.push_back(PlotCellValue(num));
                }
                else {
                    newcells.push_back(PlotCellValue());
                }
            }
            else if (getColumnType(i) == STRING && cells[i].isValue()) {
                std::stringstream ss;
                ss << cells[i].getValue();
                newcells.push_back(PlotCellValue(ss.str()));
            }
            else {
                newcells.push_back(PlotCellValue());
            }
        }
        sorted_ = false;
        rows_.push_back(PlotRowValue(this, newcells));
        updateIntervals(rows_.back());
        return true;
    }
    return false;
}

bool PlotData::insert(std::vector< PlotCellValue >& cells) {
    if (static_cast<int>(cells.size()) <= getColumnCount()) {
        for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
            // check if every cell fits column type (in-place version):
            if (getColumnType(i) == EMPTY) {
                // column type not yet defined - do that now:
                if (cells[i].isTag()) {
                    setColumnType(i, STRING);
                }
                else if (cells[i].isValue() ) {
                    setColumnType(i, NUMBER);
                }
            }
            else if (getColumnType(i) == NUMBER && cells[i].isTag()) {
                // column type is number but cell is tag - try to convert
                std::stringstream ss(cells[i].getTag());
                double num;
                if (ss >> num) {
                    cells[i].setValue(num);
                }
                else {
                    cells[i].clear();
                }
            }
            else if (getColumnType(i) == STRING && cells[i].isValue()) {
                // column type is string but cell is number - convert number to string
                std::stringstream ss;
                ss << cells[i].getValue();
                cells[i].setTag(ss.str());
            }
            // elsewise everything should be fine
        }
        sorted_ = false;
        rows_.push_back(PlotRowValue(this, cells));
        updateIntervals(rows_.back());
        return true;
    }
    return false;
}

bool PlotData::insert(const std::vector< std::pair<int, PlotCellValue > >& cells) {
    if (static_cast<int>(cells.size()) <= getColumnCount()) {
        std::vector<std::pair<int, PlotCellValue> >::const_iterator it;
        std::vector<PlotCellValue> cellsToInsert;

        for (int i=0; i< getColumnCount(); ++i)
            cellsToInsert.push_back(PlotCellValue());

        for (it = cells.begin(); it < cells.end(); ++it) {
            if ((getColumnType(it->first) == NUMBER && it->second.isValue()) ||
                (getColumnType(it->first) == STRING && it->second.isTag())) {
                    cellsToInsert[it->first] = it->second;
            }
            else if (getColumnType(it->first) == EMPTY) {
                cellsToInsert[it->first] = it->second;
                if (it->second.isTag()) {
                    setColumnType(it->first,STRING);
                }
                else if (it->second.isValue()) {
                    setColumnType(it->first,NUMBER);
                }
            }
            else if (getColumnType(it->first) == NUMBER && it->second.isTag()) {
                std::stringstream ss(it->second.getTag());
                double num;
                if (ss >> num) {
                    cellsToInsert[it->first] = PlotCellValue(num);
                }
                else {
                    cellsToInsert[it->first] = PlotCellValue();
                }
            }
            else if (getColumnType(it->first) == STRING && it->second.isValue()) {
                std::stringstream ss;
                ss << it->second.getValue();
                cellsToInsert[it->first] = PlotCellValue(ss.str());
            }
        }
    sorted_ = false;
    rows_.push_back(PlotRowValue(this, cellsToInsert));
    updateIntervals(rows_.back());
    return true;
    }
    return false;
}

std::vector<PlotRowValue>::const_iterator PlotData::lower_bound(const std::vector<PlotCellValue>& values) const {
    sortRows();
    PlotRowValue tester(this, values);
    return std::lower_bound(rows_.begin(), rows_.end(), tester);
}

int PlotData::remove(const std::vector<std::pair<int, PlotPredicate*> >& predicates) {
    std::vector<PlotRowValue>::iterator rit;
    std::vector<std::pair<int, PlotPredicate*> >::const_iterator pit;

    int count = 0;
    bool matched;
    rit = rows_.begin();
    while (rit < rows_.end()) {
        matched = false;
        for (pit = predicates.begin(); pit < predicates.end(); ++pit)
            if (pit->second->check(rit->getCellAt(pit->first))) {
                matched = true;
                break;
            }
        if (matched) {
            removePointersToCellsOfRow(*rit);
            rit = rows_.erase(rit);
            ++count;
        }
        else {
            ++rit;
        }
    }

    return count;
}

bool PlotData::rearrangeColumns(std::vector<int>& newKeyColumns, std::vector<int>& newDataColumns, PlotData& target) const {
    target.reset(static_cast<int>(newKeyColumns.size()), static_cast<int>(newDataColumns.size()));

    std::vector<PlotRowValue>::const_iterator rit;     // iterator for PlotRowValues
    std::vector<PlotRowImplicit>::const_iterator irit; // iterator for PlotRowValues
    std::vector<int>::const_iterator kit;              // iterator for key column indexs
    std::vector<int>::const_iterator dit;              // iterator for data column indexs

    // check that no index is out of bounds
    for (kit = newKeyColumns.begin(); kit < newKeyColumns.end(); ++kit)
        if (*kit < 0 || *kit >= getColumnCount())
            return false;
    for (dit = newKeyColumns.begin(); dit < newKeyColumns.end(); ++dit)
        if (*dit < 0 || *dit >= getColumnCount())
            return false;

    for (rit = rows_.begin(); rit < rows_.end(); ++rit) {
        std::vector<PlotCellValue> cellsToInsert;
        for (kit = newKeyColumns.begin(); kit < newKeyColumns.end(); ++kit) {
            cellsToInsert.push_back(rit->getCellAt(*kit));
        }
        for (dit = newKeyColumns.begin(); dit < newKeyColumns.end(); ++dit) {
            cellsToInsert.push_back(rit->getCellAt(*dit));
        }
        target.insert(cellsToInsert);
    }
    for (irit =  implicitRows_.begin(); irit < implicitRows_.end(); ++irit) {
        std::vector<PlotCellImplicit> cellsToInsert;
        for (kit = newKeyColumns.begin(); kit < newKeyColumns.end(); ++kit) {
            cellsToInsert.push_back(irit->getCellAt(*kit));
        }
        for (dit = newKeyColumns.begin(); dit < newKeyColumns.end(); ++dit) {
            cellsToInsert.push_back(irit->getCellAt(*dit));
        }
        target.insertImplicit(cellsToInsert);
    }
    return true;
}

bool PlotData::groupBy(int groupColumn, const std::vector<std::pair<int, AggregationFunction*> >& functions, PlotData& target) const {
    // check if all indexes are within bounds
    if (groupColumn < 0 || groupColumn >= getColumnCount())
        return false;
    for (std::vector<std::pair<int, AggregationFunction*> >::const_iterator it = functions.begin(); it < functions.end(); ++it)
        if (it->first < 0 || it->first >= getColumnCount())
            return false;

    target.reset(1, static_cast<int>(functions.size()));

    // we will sort and group the data into hashmaps before applying the AggregationFunctions to them
    std::map<plot_t, std::vector<plot_t>*> groupedByValues;
    std::map<std::string, std::vector<plot_t>*> groupedByTags;
    std::map<bool, std::vector<plot_t>*> groupedByNulls;

    std::map<plot_t, std::vector<plot_t>*>::iterator gbvIt;
    std::map<std::string, std::vector<plot_t>*>::iterator gbtIt;
    std::map<bool, std::vector<plot_t>*>::iterator gbnIt;
    // alright, we will start with grouping the data in our hashmaps
    int funcCount = static_cast<int>(functions.size());
    std::vector<PlotRowValue>::const_iterator rowIt;
    for (rowIt = rows_.begin(); rowIt != rows_.end(); ++rowIt) {
        const PlotCellValue& groupedCell = rowIt->getCellAt(groupColumn);

        // check whether to group by value or by tag - here we group by value
        if (groupedCell.isValue()) {
            // check if the data in this cell is already in the hashmap -> a group for it already exists
            gbvIt = groupedByValues.find(groupedCell.getValue());

            // it is not in the hashmap -> create new group
            if (gbvIt == groupedByValues.end()) {
                // first we create an array of vector of plot_t where we will store the data which shall be aggregated:
                std::vector<plot_t>* temp = new std::vector<plot_t>[funcCount]();
                for (int i=0; i< funcCount; ++i) {
                    temp[i].push_back(rowIt->getValueAt(functions[i].first));
                }
                groupedByValues.insert(std::pair<plot_t, std::vector<plot_t>*>(groupedCell.getValue(), temp));
            }
            // there is already an entry in the hashmap, so just push the values into the vector which shall be aggregated
            else {
                std::vector<plot_t>* temp = gbvIt->second;
                for (int i=0; i< funcCount; ++i) {
                    temp[i].push_back(rowIt->getValueAt(functions[i].first));
                }
            }
        }
        else if(groupedCell.isNull()) {
            gbnIt = groupedByNulls.find(groupedCell.isNull());
            if (gbnIt == groupedByNulls.end()) {
                // first we create an array of vector of plot_t where we will store the data which shall be aggregated:
                std::vector<plot_t>* temp = new std::vector<plot_t>[funcCount];
                for (int i=0; i< funcCount; ++i) {
                    temp[i].push_back(rowIt->getValueAt(functions[i].first));
                }
                groupedByNulls.insert(std::pair<bool, std::vector<plot_t>*>(groupedCell.isNull(), temp));
            }
            // there is already an entry in the hashmap, so just push the values into the vector which shall be aggregated
            else {
                std::vector<plot_t>* temp = gbnIt->second;
                for (int i=0; i< funcCount; ++i) {
                    temp[i].push_back(rowIt->getValueAt(functions[i].first));
                }
            }
        }
        // ok, we group by tag
        else {
            // check if the data in this cell is already in the hashmap -> a group for it already exists
            gbtIt = groupedByTags.find(groupedCell.getTag());

            // it is not in the hashmap -> create new group
            if (gbtIt == groupedByTags.end()) {
                // first we create an array of vector of plot_t where we will store the data which shall be aggregated:
                std::vector<plot_t>* temp = new std::vector<plot_t>[funcCount];
                for (int i=0; i< funcCount; ++i) {
                    temp[i].push_back(rowIt->getValueAt(functions[i].first));
                }
                groupedByTags.insert(std::pair<std::string, std::vector<plot_t>*>(groupedCell.getTag(), temp));
            }
            // there is already an entry in the hashmap, so just push the values into the vector which shall be aggregated
            else {
                std::vector<plot_t>* temp = gbtIt->second;
                for (int i=0; i< funcCount; ++i) {
                    temp[i].push_back(rowIt->getValueAt(functions[i].first));
                }
            }
        }
    }

    // grouping is done, now, we will insert the grouped data into the target PlotData
    // first insert all groupedByValues rows
    for (gbvIt = groupedByValues.begin(); gbvIt != groupedByValues.end(); ++gbvIt) {
        std::vector<PlotCellValue> cells;
        cells.push_back(PlotCellValue(gbvIt->first));

        for (int i=0; i<funcCount; ++i) {
            plot_t val = functions[i].second->evaluate(gbvIt->second[i]);
            if (val != val){
                cells.push_back(PlotCellValue());
            }
            else {
                cells.push_back(PlotCellValue(functions[i].second->evaluate(gbvIt->second[i])));
            }
        }

        target.insert(cells);
    }

    // now insert all groupedByTags rows
    for (gbtIt = groupedByTags.begin(); gbtIt != groupedByTags.end(); ++gbtIt) {
        std::vector<PlotCellValue> cells;
        cells.push_back(PlotCellValue(gbtIt->first));

        for (int i=0; i<funcCount; ++i) {
            plot_t val = functions[i].second->evaluate(gbtIt->second[i]);
            if (val != val){
                cells.push_back(PlotCellValue());
            }
            else {
                cells.push_back(PlotCellValue(functions[i].second->evaluate(gbtIt->second[i])));
            }
        }

        target.insert(cells);
    }
        // now insert all groupedByNulls rows
    for (gbnIt = groupedByNulls.begin(); gbnIt != groupedByNulls.end(); ++gbnIt) {
        std::vector<PlotCellValue> cells;
        cells.push_back(PlotCellValue());

        for (int i=0; i<funcCount; ++i) {
            plot_t val = functions[i].second->evaluate(gbnIt->second[i]);
            if (val != val){
                cells.push_back(PlotCellValue());
            }
            else {
                cells.push_back(PlotCellValue(functions[i].second->evaluate(gbnIt->second[i])));
            }
        }

        target.insert(cells);
    }

    // finally free all data which is not necessary anymore:
    while (groupedByValues.size() > 0) {
        gbvIt = groupedByValues.begin();
        delete [] gbvIt->second;
        groupedByValues.erase(gbvIt);
    }
    while (groupedByTags.size() > 0) {
        gbtIt = groupedByTags.begin();
        delete [] gbtIt->second;
        groupedByTags.erase(gbtIt);
    }
    while (groupedByNulls.size() > 0) {
        gbnIt = groupedByNulls.begin();
        delete [] gbnIt->second;
        groupedByNulls.erase(gbnIt);
    }

    for (int i = 0; i < static_cast<int>(functions.size()); ++i) {
        target.setColumnLabel(i+1,getColumnLabel(functions.at(i).first));
    }
    target.setColumnLabel(0,getColumnLabel(groupColumn));
    return true;
}

bool PlotData::mergeWith(
        const PlotData& otherPlotData,
        std::vector<std::pair<int, int> >& keyColumns,
        std::vector<std::pair<int, int> >& dataColumns,
        PlotData& target) const {

    target.reset(static_cast<int>(keyColumns.size()), static_cast<int>(dataColumns.size()));

    std::vector<PlotRowValue>::const_iterator rit;
    std::vector<PlotRowImplicit>::const_iterator irit;
    std::vector<std::pair<int, int> >::const_iterator kit;
    std::vector<std::pair<int, int> >::const_iterator dit;

    // first insert all rows from *this into target
    for (rit = rows_.begin(); rit < rows_.end(); ++rit) {
        std::vector<PlotCellValue> cellsToInsert;
        for (kit = keyColumns.begin(); kit < keyColumns.end(); ++kit)
            cellsToInsert.push_back(rit->getCellAt(kit->first));

        for (dit = dataColumns.begin(); dit < dataColumns.end(); ++dit)
            cellsToInsert.push_back(rit->getCellAt(dit->first));

        target.insert(cellsToInsert);
    }
    for (irit = implicitRows_.begin(); irit < implicitRows_.end(); ++irit) {
        std::vector<PlotCellImplicit> cellsToInsert;
        for (kit = keyColumns.begin(); kit < keyColumns.end(); ++kit)
            cellsToInsert.push_back(irit->getCellAt(kit->first));

        for (dit = dataColumns.begin(); dit < dataColumns.end(); ++dit)
            cellsToInsert.push_back(irit->getCellAt(dit->first));

        target.insertImplicit(cellsToInsert);
    }

    // now insert alls rows from otherPlotData into target
    for (rit = otherPlotData.rows_.begin(); rit < otherPlotData.rows_.end(); ++rit) {
        std::vector<PlotCellValue> cellsToInsert;
        for (kit = keyColumns.begin(); kit < keyColumns.end(); ++kit)
            cellsToInsert.push_back(rit->getCellAt(kit->second));

        for (dit = dataColumns.begin(); dit < dataColumns.end(); ++dit)
            cellsToInsert.push_back(rit->getCellAt(dit->second));

        target.insert(cellsToInsert);
    }
    for (size_t i = 0; i < keyColumns.size(); ++i) {
        target.setColumnLabel(static_cast<int>(i), getColumnLabel(keyColumns.at(i).first));
    }
    for (size_t i = 0; i < dataColumns.size(); ++i) {
        target.setColumnLabel(static_cast<int>(keyColumns.size()+i),getColumnLabel(dataColumns.at(i).first));
    }

    return true;
}

bool PlotData::mergeWith(const PlotData& otherPlotData, PlotData& target, MergeType mergeType) const {
    if (otherPlotData.getColumnCount() == getColumnCount() && (mergeType == NOSELECTION || mergeType == NEWROWS)){
        bool columnMatch = true;
        for (int i=0; i< getColumnCount();  ++i){
            if (otherPlotData.getColumnLabel(i) != getColumnLabel(i)){
                    columnMatch = false;
                    break;
            }
        }
        if(columnMatch){
            std::vector<std::pair<int, int> > keyColumns;
            std::vector<std::pair<int, int> > dataColumns;
            for(int i=0; i< getKeyColumnCount();  ++i){
                keyColumns.push_back(std::pair<int,int>(i,i));
            }
            for(int i=0; i< getDataColumnCount();  ++i){
                dataColumns.push_back(std::pair<int,int>(i+ getKeyColumnCount(),i+ getKeyColumnCount()));
            }
            mergeWith(otherPlotData, keyColumns, dataColumns, target);
            return true;
        }
        else if (mergeType == NEWROWS)
            return false;
    }
    else if (mergeType == NEWROWS)
        return false;

    if (keyColumnCount_ == otherPlotData.getKeyColumnCount() && (mergeType == NOSELECTION || mergeType == IGNORECOLUMNLABELS)) {
        bool columnMatch = true;
        for (int i = 0; i < keyColumnCount_; ++i) {
            if (getColumnType(i) != otherPlotData.getColumnType(i) || (getColumnLabel(i) != otherPlotData.getColumnLabel(i) && mergeType == NOSELECTION)) {
                columnMatch = false;
                break;
            }
        }
        if (columnMatch) {
            sortRows();
            otherPlotData.sortRows();
            target.reset(getKeyColumnCount(),getDataColumnCount() + otherPlotData.getDataColumnCount());
            for (int i = 0; i < getColumnCount(); ++i) {
                target.setColumnLabel(i,getColumnLabel(i));
            }
            for (int i = 0; i < otherPlotData.getDataColumnCount(); ++i) {
                target.setColumnLabel(getColumnCount() + i, otherPlotData.getColumnLabel(getKeyColumnCount() + i));
            }
            int selfcounter = 0;
            int othercounter = 0;
            std::vector<PlotCellValue> row;
            std::vector<PlotCellValue> otherrow;
            std::vector<PlotCellValue> cells;
            row.resize(getKeyColumnCount());
            otherrow.resize(getKeyColumnCount());
            while (selfcounter < getRowsCount() || othercounter < otherPlotData.getRowsCount()) {
                cells.clear();
                cells.resize(getColumnCount() + otherPlotData.getDataColumnCount());
                for (int i = 0; i < getKeyColumnCount(); ++i) {
                    if (selfcounter < getRowsCount())
                        row[i] = rows_[selfcounter].cells_[i];
                    if (othercounter < otherPlotData.getRowsCount())
                        otherrow[i] = otherPlotData.getRow(othercounter).cells_[i];
                }
                if (row == otherrow) {
                    for (int i = 0; i < getColumnCount(); ++i) {
                        cells[i] = rows_[selfcounter].cells_[i];
                    }
                    for (int i = 0; i < otherPlotData.getDataColumnCount(); ++i) {
                        cells[i+getColumnCount()] = otherPlotData.getRow(othercounter).cells_[i+getKeyColumnCount()];
                    }
                    target.insert(cells);
                    ++selfcounter;
                    ++othercounter;
                }
                else if ((row > otherrow && othercounter < otherPlotData.getRowsCount()) || (selfcounter >= getRowsCount())) {
                    for (int i = 0; i < getKeyColumnCount(); ++i) {
                        cells[i] = otherPlotData.getRow(othercounter).cells_[i];
                    }
                    for (int i = 0; i < otherPlotData.getDataColumnCount(); ++i) {
                        cells[i+getColumnCount()] = otherPlotData.getRow(othercounter).cells_[i+getKeyColumnCount()];
                    }
                    target.insert(cells);
                    ++othercounter;
                }
                else if ((row < otherrow && selfcounter < getRowsCount()) || (othercounter >= otherPlotData.getRowsCount())) {
                    for (int i = 0; i < getColumnCount(); ++i) {
                        cells[i] = rows_[selfcounter].cells_[i];
                    }
                    target.insert(cells);
                    ++selfcounter;
                }
            }
            return true;
        }
    }
    else if (mergeType == IGNORECOLUMNLABELS)
        return false;

    if (mergeType == NOSELECTION || mergeType == NEWCOLUMNS) {
        std::vector<PlotCellValue> cells;
        std::vector<PlotCellValue> row;
        std::vector<PlotCellValue> otherRow;
        int count = std::max<int>(getRowsCount(),otherPlotData.getRowsCount());
        target.reset(getKeyColumnCount(), getDataColumnCount() + otherPlotData.getColumnCount());
        for (int i=0; i< count; ++i){ // Anzahl der rows
            cells.clear();
            if (i < getRowsCount()) {
                row = getRow(i).getCells();
                for( size_t j=0; j< row.size(); ++j){
                    cells.push_back(row.at(j));
                }
            }
            else {
                for(int j=0; j< getColumnCount(); ++j){
                    cells.push_back(PlotCellValue());
                }
            }
            if (i < otherPlotData.getRowsCount()) {
                otherRow = otherPlotData.getRow(i).getCells();
                for( size_t j=0; j< otherRow.size(); ++j){
                    cells.push_back(otherRow.at(j));
                }
            }
            target.insert(cells);
        }
        for (int i = 0; i < getColumnCount(); ++i) {
            target.setColumnLabel(i,getColumnLabel(i));
        }
        for (int i = 0; i < otherPlotData.getColumnCount(); ++i) {
            target.setColumnLabel(i + getColumnCount(),otherPlotData.getColumnLabel(i));
        }
        std::vector<PlotCellImplicit> cellsToInsert;
        count = std::max<int>(getImplicitRowsCount(),otherPlotData.getImplicitRowsCount());
        for (int j = 0; j < count; ++j) {
            cellsToInsert.clear();
            if (getImplicitRowsCount() > j) {
                for (int i = 0; i < getColumnCount(); ++i) {
                    cellsToInsert.push_back(implicitRows_[j].getCellAt(i));
                }
            }
            else {
                for (int i = 0; i < getColumnCount(); ++i) {
                    cellsToInsert.push_back(PlotCellImplicit(0,&target,i));
                }
            }
            if (otherPlotData.getImplicitRowsCount() > j) {
                for (int i = 0; i < otherPlotData.getColumnCount(); ++i) {
                    cellsToInsert.push_back(otherPlotData.getImplicitRow(j).getCellAt(i));
                }
            }
            else {
                for (int i = 0; i < otherPlotData.getColumnCount(); ++i) {
                    cellsToInsert.push_back(PlotCellImplicit(0,&target,getColumnCount()+i));
                }
            }
            target.insertImplicit(cellsToInsert);
        }
        return true;
    }
    else if (mergeType == NEWCOLUMNS)
        return false;
    return false;
}

bool PlotData::insertImplicit(const std::vector<std::pair<int, AggregationFunction*> >& functions) {
    std::vector<AggregationFunction*> sortedFuncs;
    for (int i=0; i<getColumnCount(); ++i) {
        sortedFuncs.push_back(NULL);
    }

    std::vector<std::pair<int, AggregationFunction*> >::const_iterator it;
    for (it = functions.begin(); it < functions.end(); ++it) {
        sortedFuncs[it->first] = it->second;
    }

    PlotRowImplicit pri = PlotRowImplicit(this, sortedFuncs);
    implicitRows_.push_back(pri);
    return true;
}

bool PlotData::insertImplicit(const std::vector<PlotCellImplicit> &implicitrow) {
    if (getColumnCount() > static_cast<int>(implicitrow.size()))
        return false;
    std::vector<AggregationFunction*> sortedFuncs;
    for (int i=0; i<getColumnCount(); ++i) {
        sortedFuncs.push_back(implicitrow.at(i).getAggregationFunction());
    }

    PlotRowImplicit pri = PlotRowImplicit(this, sortedFuncs);
    implicitRows_.push_back(pri);
    return true;
}

Interval<plot_t> PlotData::getInterval(int column) const {
    if (getColumnCount() <= column || rows_.empty() || columns_[column].type_ == EMPTY) {
        Interval<plot_t> toReturn(0, 0, true, true);
        return toReturn;
    }
    return intervals_[column];
}

Interval<plot_t> PlotData::getSumInterval(const std::vector< int >& column) const {
    if (rows_.empty()) {
        return Interval<plot_t>(0, 0, true, true);
    }

    std::vector<int> numberColumn;
    //only keep number columns
    std::vector<int>::const_iterator colIt;
    for (colIt = column.begin(); colIt < column.end(); ++colIt) {
        if (columns_[*colIt].type_ == NUMBER)
            numberColumn.push_back(*colIt);
    }

    std::vector<PlotRowValue>::const_iterator rowIt = rows_.begin();
    plot_t min = 0;
    plot_t max = 0;

    for ( ; rowIt < rows_.end(); ++rowIt) {
        plot_t temp = 0;
        for (colIt = numberColumn.begin(); colIt < numberColumn.end(); ++colIt) {
            if (!rowIt->getCellAt(*colIt).isNull())
                temp += fabs(rowIt->getValueAt(*colIt));
        }
        if (temp > max)
            max = temp;
    }

    return Interval<plot_t>(min, max);
}

std::vector<PlotRowValue>::const_iterator  PlotData::getRowsBegin() const {
//    sortRows();
    return rows_.begin();
}

std::vector<PlotRowImplicit>::const_iterator PlotData::getImplicitRowsBegin() const {
    return implicitRows_.begin();
}

std::vector<PlotRowValue>::const_iterator  PlotData::getRowsEnd() const {
//    sortRows();
    return rows_.end();
}

std::vector<PlotRowImplicit>::const_iterator PlotData::getImplicitRowsEnd() const {
    return implicitRows_.end();
}

const PlotRowValue& PlotData::getRow(int row) const {
    tgtAssert((row >= 0 && row < static_cast<int>(rows_.size())), "PlotData::getRow(): row out of bounds");
//    sortRows();
    return rows_[row];
}

const PlotRowImplicit& PlotData::getImplicitRow(int row) const {
    tgtAssert((row >= 0 && row < static_cast<int>(rows_.size())), "PlotData::getImplicitRow(): row out of bounds");
    return implicitRows_[row];
}

int PlotData::getRowsCount() const {
    return static_cast<int>(rows_.size());
}

int PlotData::getImplicitRowsCount() const {
    return static_cast<int>(implicitRows_.size());
}

bool PlotData::rowsEmpty() const {
    return rows_.empty();
}

bool PlotData::implicitRowsEmpty() const {
    return rows_.empty();
}

bool PlotData::isHighlighted(const tgt::ivec2& cellPosition) const {
    // check if cellPosition is valid
    if ((cellPosition.x >= -1 && cellPosition.x < getRowsCount()) && (cellPosition.y >= -1 && cellPosition.y < getColumnCount())) {
        // (-1,-1) is not a valid cellPosition
        if (cellPosition.x == -1 && cellPosition.y == -1)
            return false;
        // check the highlight_state in one cell
        else if (cellPosition.x != -1 && cellPosition.y != -1)
            return rows_[cellPosition.x].getCellAt(cellPosition.y).isHighlighted();
        // check the highlight_state in all cells from that column
        else if (cellPosition.x == -1) {
            bool highlighted = true;
            // the column is highlighted if all cells are highlighted
            for (std::vector<PlotRowValue>::const_iterator it = rows_.begin(); it < rows_.end(); ++it) {
                if (!it->getCellAt(cellPosition.y).isHighlighted()) {
                    highlighted = false;
                    break;
                }
            }
            return highlighted;
        }
        // check the highlight_state in all cells from that row
        else {  // if (cellPosition.y == -1) {
            bool highlighted = true;
            // the row is highlighted if all cells are highlighted
            for (int i = 0; i < getColumnCount(); ++i) {
                if (!rows_[cellPosition.x].getCellAt(i).isHighlighted()) {
                    highlighted = false;
                    break;
                }
            }
            return highlighted;
        }
    }
    else
        return false;
}

void PlotData::toggleHighlight(const tgt::ivec2& cellPosition, bool additive) {
    setHighlight(cellPosition, !isHighlighted(cellPosition), additive);
}

void PlotData::clearHighlights() {
    for (std::set<PlotCellValue*>::iterator it = highlightedCells_.begin(); it != highlightedCells_.end(); ++it) {
        (*it)->setHighlighted(false);
    }
    highlightedCells_.clear();
}

void PlotData::setHighlight(const tgt::ivec2& cellPosition, bool value, bool additive) {
    if (! additive)
        clearHighlights();
    // check for sanity of position
    if ((cellPosition.x >= -1 && cellPosition.x < getRowsCount())
        && (cellPosition.y >= -1 && cellPosition.y < getColumnCount())
        && !(cellPosition.x == -1 && cellPosition.y == -1)) {
        // check if we have to change the state
        if (value != isHighlighted(cellPosition)) {
            // single cell
            if (cellPosition.x != -1 && cellPosition.y != -1) {
                rows_[cellPosition.x].setHighlighted(cellPosition.y, value);
                if (value)
                    highlightedCells_.insert(&(rows_[cellPosition.x].cells_[cellPosition.y]));
                else
                    highlightedCells_.erase(&(rows_[cellPosition.x].cells_[cellPosition.y]));
            }
            // entire column
            else if (cellPosition.x == -1) {
                for (std::vector<PlotRowValue>::iterator it = rows_.begin(); it < rows_.end(); ++it) {
                    // check if we have to change the state of single cell
                    if (it->getCellAt(cellPosition.y).isHighlighted() != value) {
                        it->setHighlighted(cellPosition.y, value);
                        if (value)
                            highlightedCells_.insert(&(it->cells_[cellPosition.y]));
                        else
                            highlightedCells_.erase(&(it->cells_[cellPosition.y]));
                    }
                }
            }
            // entire row
            else if (cellPosition.y == -1) {
                for (int i = 0; i < getColumnCount(); ++i) {
                    // check if we have to change the state of single cell
                    if (rows_[cellPosition.x].getCellAt(i).isHighlighted() != value) {
                        rows_[cellPosition.x].setHighlighted(i, value);
                        if (value)
                            highlightedCells_.insert(&(rows_[cellPosition.x].cells_[i]));
                        else
                            highlightedCells_.erase(&(rows_[cellPosition.x].cells_[i]));
                    }
                }
            }
        }
    }
}

void PlotData::deleteImplicitRows() {
    std::vector<PlotRowImplicit>::const_iterator imit;
    std::vector<PlotCellImplicit> implicitCells;
    AggregationFunction* aggf;
    for (imit = implicitRows_.begin(); imit < implicitRows_.end(); ++imit) {
        implicitCells = imit->getCells();
        for (size_t i = 0; i < implicitCells.size(); ++i) {
            aggf = implicitCells.at(i).getAggregationFunction();
            delete aggf;
        }
    }
}

void PlotData::removePointersToCellsOfRow(PlotRowValue& row) {
    for (std::vector<PlotCellValue>::iterator cit = row.cells_.begin(); cit < row.cells_.end(); ++cit) {
        PlotCellValue* pCell = &(*cit);
        highlightedCells_.erase(pCell);
    }
}

void PlotData::reset(int keyColumnCount, int dataColumnCount) {
    highlightedCells_.clear();
    rows_.clear();
    deleteImplicitRows();
    implicitRows_.clear();
    intervals_.clear();
    sorted_ = false;
    PlotBase::reset(keyColumnCount, dataColumnCount);
    intervals_.resize(getColumnCount());
}

void PlotData::updateIntervals(const PlotRowValue& row) {
    std::vector< Interval<plot_t> >::iterator iit = intervals_.begin();
    std::vector< PlotCellValue >::const_iterator rit = row.getCells().begin();
    int column = 0;
    for (; rit < row.getCells().end() && iit < intervals_.end(); ++iit, ++rit, ++column) {
        if (columns_[column].type_ == PlotBase::STRING)
            *iit = Interval<plot_t>(0, static_cast<plot_t>(rows_.size()) - 1, false, false);
        else if (rit->isValue()) {
            const plot_t& value = rit->getValue();
            iit->nibble(value);
        }
    }
}

void PlotData::sortRows() const {
    if (! sorted_ && rows_.size() > 0) {
        // in many cases the data is sorted by construction, so we check that
        bool sorted = true;
        for (std::vector< PlotRowValue >::const_iterator it = rows_.begin(); it+1 < rows_.end(); ++it) {
            if ( *(it+1) < *it ) {
                sorted = false;
                break;
            }
        }
        if (sorted)
            sorted_ = true;
        else {
            // sortRows() shall be callable by const member functions so it has to be const itself
            // Unfortunately std::sort does not offer bitwise-constness but in the semantic way of constness
            // (the items of this plot data are still the same) sortRows() is const. So we do a const cast
            // here to be allowed to call std::sort
            PlotData* foo = const_cast<PlotData*>(this);
            std::sort(foo->rows_.begin(), foo->rows_.end());
            sorted_ = true;
        }
    }
}

bool PlotData::sorted() const {
    return sorted_;
}

bool PlotData::isIndexColumn(const PlotData& pData, int column) {
    return (column == 0 && pData.getColumnCount() > 0 && pData.getColumnLabel(0) == "Index" && pData.getColumnType(0) == NUMBER);
}

} // namespace
