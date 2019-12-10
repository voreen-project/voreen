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

#include "plotrow.h"
#include "plotdata.h"
#include "tgt/assert.h"


namespace voreen {

// PlotRowValue methods --------------------------------------------------

PlotRowValue::PlotRowValue(const PlotData* parent, const std::vector<PlotCellValue>& cells)
    : cells_(cells)
    , parent_(parent)
{
    tgtAssert(parent, "PlotRowValue initialized with NULL pointed parent");

    // make sure, cells_ has parent_->getColumnCount() entries
    cells_.resize(parent_->getColumnCount());
}

PlotRowValue::PlotRowValue(const PlotData* parent, std::vector<PlotCellValue>& cells)
    : cells_(cells)
    , parent_(parent)
{
    tgtAssert(parent, "PlotRowValue initialized with NULL pointed parent");

    // make sure, cells_ has parent_->getColumnCount() entries
    cells_.resize(parent_->getColumnCount());
}

PlotRowValue::~PlotRowValue() {
}

const std::vector<PlotCellValue>& PlotRowValue::getCells() const {
    return cells_;
}

const PlotCellValue& PlotRowValue::getCellAt(int column) const {
    tgtAssert((column >= 0 && column < parent_->getColumnCount()), "PlotRowValue::getCellAt: column out of bounds.");
    return cells_[column];
}

plot_t PlotRowValue::getValueAt(int column) const {
    tgtAssert((column >= 0 && column < parent_->getColumnCount()), "PlotRowValue::getValueAt: column out of bounds.");
    return getCellAt(column).getValue();
}

std::string PlotRowValue::getTagAt(int column) const {
    tgtAssert((column >= 0 && column < parent_->getColumnCount()), "PlotRowValue::getTagAt: column out of bounds.");
    return getCellAt(column).getTag();
}

void PlotRowValue::setHighlighted(int column, bool value) {
    tgtAssert((column >= 0 && column < parent_->getColumnCount()), "PlotRowValue::setHighlighted: column out of bounds.");
    cells_[column].setHighlighted(value);
}

bool PlotRowValue::operator<(const PlotRowValue& rhs) const {
    int currentIndex = 0;

    // having some fun with lexicographic comparison:
    while (currentIndex < parent_->getKeyColumnCount()) {
        // most cells will be filled with values, so check that first
        if (cells_[currentIndex].isValue() && rhs.getCellAt(currentIndex).isValue()) {
            if (cells_[currentIndex].getValue() < rhs.getValueAt(currentIndex)) {
                return true;
            }
            if (cells_[currentIndex].getValue() > rhs.getValueAt(currentIndex)) {
                return false;
            }
        }
        // now check for tags
        else if (cells_[currentIndex].isTag() && rhs.getCellAt(currentIndex).isTag()) {
            if (cells_[currentIndex].getTag() < rhs.getTagAt(currentIndex)) {
                return true;
            }
            if (cells_[currentIndex].getTag() > rhs.getTagAt(currentIndex)) {
                return false;
            }
        }
        // check if exactly one cell is null
        // in that case the null cell is considered smaller, if both are null they are considered as equal
        // and we continue with the next column
        else if (cells_[currentIndex].isNull() && !rhs.getCellAt(currentIndex).isNull()) {
            return true;
        }
        else if (!cells_[currentIndex].isNull() && rhs.getCellAt(currentIndex).isNull()) {
            return false;
        }
        // If we enter here none of the cells are null, one is a value cell, the other is a tag cell.
        // The tag-containing cell is considered as smaller.
        else
            return cells_[currentIndex].isTag();

        // if we reach here the current cells in the current column are equal, try the next column
        ++currentIndex;
    }
    return false;
}

// PlotRowImplicit methods --------------------------------------------------

PlotRowImplicit::PlotRowImplicit(const PlotData* parent, const std::vector<AggregationFunction*>& functions)
    : parent_(parent)
{
    tgtAssert(parent, "PlotRowImplicit initialized with NULL pointed parent");

    for (int i = 0; (i < parent->getColumnCount() && i < static_cast<int>(functions.size())); ++i) {
        PlotCellImplicit pci(functions[i], parent, i);
        cells_.push_back(pci);
    }
}

PlotRowImplicit::~PlotRowImplicit() {
}


const std::vector<PlotCellImplicit>& PlotRowImplicit::getCells() const {
    return cells_;
}

const PlotCellImplicit& PlotRowImplicit::getCellAt(int column) const {
    tgtAssert((column >= 0 && column < parent_->getColumnCount()), "PlotCellImplicit& PlotRowImplicit::getCellAt(): Column out of bounds.");
    return cells_[column];
}

plot_t PlotRowImplicit::getValueAt(int column) const {
    tgtAssert((column >= 0 && column < parent_->getColumnCount()), "PlotRowImplicit::getValueAt(): Column out of bounds.");
    return getCellAt(column).getValue();
}

}
