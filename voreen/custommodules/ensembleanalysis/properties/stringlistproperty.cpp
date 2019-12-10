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

#include "stringlistproperty.h"

#include "voreen/core/utils/stringutils.h"

namespace voreen {

StringListProperty::StringListProperty(const std::string& id, const std::string& guiText,
                       int invalidationLevel, Property::LevelOfDetail lod)
    : TemplateProperty<std::vector<int>>(id, guiText, std::vector<int>(), invalidationLevel, lod)
    , neededTableUpdates_(UPDATE_ALL)
{
}

StringListProperty::StringListProperty()
{
}

void StringListProperty::addWidget(PropertyWidget* widget) {
    Property::addWidget(widget);
    neededTableUpdates_ = StringListProperty::UPDATE_ALL;
    updateWidgets();
}

void StringListProperty::updateWidgets() {
    Property::updateWidgets();
    neededTableUpdates_ = StringListProperty::UPDATE_NONE;
}

void StringListProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    s.serialize("selectedRow", selectedRows_);
    s.serialize("values", values_);
}

void StringListProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    s.optionalDeserialize("selectedRow", selectedRows_, std::vector<int>());
    s.deserialize("values", values_);
    neededTableUpdates_ = UPDATE_ALL;
    updateWidgets();
}

//-------------------------------------------------------------------------------------------
//          Row    Functions
//-------------------------------------------------------------------------------------------
int StringListProperty::getNumRows() const {
    return static_cast<int>(values_.size());
}

const std::vector<int>& StringListProperty::getSelectedRowIndices() const {
    return selectedRows_;
}

void StringListProperty::setSelectedRowIndices(const std::vector<int>& rowIndices) {
    selectedRows_.clear();
    // Ignore invalid selection silently.
    for(int rowIdx : rowIndices) {
        if(rowIdx >= 0 && rowIdx < getNumRows()) {
            selectedRows_.push_back(rowIdx);
        }
    }
    neededTableUpdates_ |= UPDATE_SELECTION;
    invalidate();
}

void StringListProperty::reset() {
    values_.clear();
    selectedRows_.clear();
    neededTableUpdates_ = UPDATE_ALL;
    invalidate();
}

void StringListProperty::addRow(const std::string& row, const tgt::vec3& color) {
    values_.push_back(std::make_pair(row, color));
    neededTableUpdates_ |= UPDATE_ROWS;
    updateWidgets();
}

}   // namespace
