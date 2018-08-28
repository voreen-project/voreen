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

#include "voreen/core/properties/string/stringtableproperty.h"

#include "voreen/core/utils/stringutils.h"

namespace voreen {

StringTableProperty::StringTableProperty(const std::string& id, const std::string& guiText,
                       const unsigned int columnCount, int invalidationLevel, Property::LevelOfDetail lod)
    : Property(id, guiText, invalidationLevel, lod)
    , columnCount_(columnCount), selectedRow_(-1)
    , neededTableUpdates_(UPDATE_ALL)
{
    tgtAssert(columnCount_ > 0, "zero columns are not allowed");
    for(size_t i = 0; i < columnCount_; i++)
        columnLabels_.push_back("Column " + itos(i));
}

StringTableProperty::StringTableProperty()
    : columnCount_(0), selectedRow_(-1)
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

    s.deserialize("columnCount",  columnCount_);
    s.optionalDeserialize("selectedRow", selectedRow_,-1);
    s.deserialize("columnLabels", columnLabels_);
    s.deserialize("values",       values_);
    neededTableUpdates_ = UPDATE_ALL;
    updateWidgets();
}

//-------------------------------------------------------------------------------------------
//          Column Functions
//-------------------------------------------------------------------------------------------
unsigned int StringTableProperty::getNumColumns() const {
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
unsigned int StringTableProperty::getNumRows() const {
    return static_cast<unsigned int>(values_.size());
}

int StringTableProperty::getSelectedRowIndex() const {
    return selectedRow_;
}

void StringTableProperty::setSelectedRowIndex(int rowIndex) {
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
    tgtAssert(row.size() == columnCount_, "row elements do not match column count");
    values_.push_back(row);
    neededTableUpdates_ |= UPDATE_ROWS;
    updateWidgets();
}

}   // namespace
