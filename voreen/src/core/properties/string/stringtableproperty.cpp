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

#include "voreen/core/properties/string/stringtableproperty.h"

#include "voreen/core/utils/stringutils.h"

namespace voreen {

StringTableProperty::StringTableProperty(const std::string& id, const std::string& guiText,
                       int columnCount, int invalidationLevel, Property::LevelOfDetail lod)
    : Property(id, guiText, invalidationLevel, lod)
    , columnCount_(columnCount)
    , selectedRow_(-1)
    , neededTableUpdates_(UPDATE_ALL)
{
    tgtAssert(columnCount_ > 0, "at least one column is required");
    for(int i = 0; i < columnCount_; i++)
        columnLabels_.push_back("Column " + itos(i));
}

StringTableProperty::StringTableProperty()
{
}

void StringTableProperty::addWidget(PropertyWidget* widget) {
    Property::addWidget(widget);
    neededTableUpdates_ = StringTableProperty::UPDATE_ALL;
    updateWidgets();
}

void StringTableProperty::updateWidgets() {
    Property::updateWidgets();
    neededTableUpdates_ = StringTableProperty::UPDATE_NONE;
}

void StringTableProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    s.serialize("columnCount",  columnCount_);
    s.serialize("selectedRow", selectedRow_);
    s.serialize("columnLabels", columnLabels_);
    s.serialize("values",       values_);
}

void StringTableProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    int newColumnCount = columnCount_;
    s.deserialize("columnCount",  newColumnCount);
    if(newColumnCount == columnCount_) {
        int selectedRow = -1;
        s.optionalDeserialize("selectedRow", selectedRow,-1);
        setSelectedRowIndex(selectedRow);
        s.deserialize("columnLabels", columnLabels_);
        s.deserialize("values",       values_);
    }
    else {
        LWARNING("Number of columns changed, discarding labels and content");
        selectedRow_ = -1;
    }

    neededTableUpdates_ = UPDATE_ALL;
    updateWidgets();
}

//-------------------------------------------------------------------------------------------
//          Column Functions
//-------------------------------------------------------------------------------------------
int StringTableProperty::getNumColumns() const {
    return columnCount_;
}

void StringTableProperty::setColumnLabel(size_t number, const std::string& label) {
    tgtAssert(number < columnLabels_.size(), "Not enough columns");
    columnLabels_[number] = label;
    neededTableUpdates_ |= UPDATE_LABELS;
    updateWidgets();
}

const std::vector<std::string>& StringTableProperty::getColumnLabels() const {
    return columnLabels_;
}

//-------------------------------------------------------------------------------------------
//          Row    Functions
//-------------------------------------------------------------------------------------------
int StringTableProperty::getNumRows() const {
    return static_cast<int>(values_.size());
}

int StringTableProperty::getSelectedRowIndex() const {
    return selectedRow_;
}

void StringTableProperty::setSelectedRowIndex(int rowIndex) {

    // Ignore invalid selection silently.
    if(rowIndex >= getNumRows()) {
        rowIndex = -1;
    }

    if(rowIndex != selectedRow_) {
        selectedRow_ = rowIndex;
        neededTableUpdates_ |= UPDATE_SELECTION;
        invalidate();
    }
}

void StringTableProperty::reset() {
    values_.clear();
    selectedRow_ = -1;
    neededTableUpdates_ = UPDATE_ALL;
    invalidate();
}

void StringTableProperty::addRow(const std::vector<std::string>& row) {
    tgtAssert(row.size() == static_cast<size_t>(columnCount_), "row elements do not match column count");
    values_.push_back(row);
    neededTableUpdates_ |= UPDATE_ROWS;
    updateWidgets();
}

void StringTableProperty::removeRow(int rowIndex) {
    tgtAssert(rowIndex < getNumRows(), "row index out of bounds");
    values_.erase(values_.begin() + rowIndex);
    neededTableUpdates_ |= UPDATE_ROWS;

    if(selectedRow_ >= 0
       && (selectedRow_ == static_cast<int>(getNumRows())
           || rowIndex > selectedRow_)) {
        selectedRow_--;
        neededTableUpdates_ |= UPDATE_SELECTION;
    }

    updateWidgets();
}

}   // namespace
