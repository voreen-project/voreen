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

#include "plotselection.h"
#include "../properties/plotselectionproperty.h"

#include "voreen/core/io/serialization/serialization.h"

#include <sstream>

namespace voreen {

PlotSelection::PlotSelection()
    : isTablePositionFlag_(false)
    , tablePosition_(-1, -1)
{}

PlotSelection::PlotSelection(tgt::ivec2 tablePosition)
    : isTablePositionFlag_(true)
    , tablePosition_(tablePosition)
{}

PlotSelection::PlotSelection(const PlotSelection& rhs)
    : Serializable()
    , isTablePositionFlag_(rhs.isTablePositionFlag_)
    , tablePosition_(rhs.tablePosition_)
{
    for (std::vector< std::pair<int, PlotPredicate*> >::const_iterator it = rhs.selection_.begin(); it < rhs.selection_.end(); ++it) {
        PlotPredicate* p = 0;
        if (it->second != 0)
            p = it->second->clone();
        selection_.push_back(std::make_pair(it->first, p));
    }
}

PlotSelection::~PlotSelection() {
    for (std::vector< std::pair<int, PlotPredicate*> >::iterator it = selection_.begin(); it < selection_.end(); ++it) {
        delete it->second;
    }
}

PlotSelection& PlotSelection::operator=(PlotSelection rhs) {
    std::swap(isTablePositionFlag_, rhs.isTablePositionFlag_);
    std::swap(tablePosition_, rhs.tablePosition_);
    std::swap(selection_, rhs.selection_); // no worries: all formerly pointed objects will be deleted on ~rhs()
    return *this;
}

void PlotSelection::serialize(Serializer& s) const {
    s.serialize("isTablePositionFlag", isTablePositionFlag_);
    s.serialize("tablePosition", tablePosition_);
    s.serialize("selection", selection_);
}

void PlotSelection::deserialize(Deserializer& s) {
    s.deserialize("isTablePositionFlag", isTablePositionFlag_);
    s.deserialize("tablePosition", tablePosition_);
    s.deserialize("selection", selection_);
}

void PlotSelection::addPredicate(int column, const PlotPredicate* pred) {
    PlotPredicate* p = 0;
    if (pred != 0) {
        p = pred->clone();
    }
    selection_.push_back(std::make_pair(column, p));
}

void PlotSelection::removeRow(size_t index) {
    tgtAssert(index < selection_.size(), "PlotSelection::removeRow(): Index out of bounds");
    if (selection_[index].second != 0)
        delete selection_[index].second;
    selection_.erase(selection_.begin()+index);
}

bool PlotSelection::isTablePosition() const {
    return isTablePositionFlag_;
}

const tgt::ivec2& PlotSelection::getTablePosition() const {
    tgtAssert(isTablePositionFlag_, "No table position");
    return tablePosition_;
}

const std::vector<std::pair<int, PlotPredicate*> >& PlotSelection::getSelection() const {
    tgtAssert(!isTablePositionFlag_, "Is table position");
    return selection_;
}

void PlotSelection::setColumn(size_t index, int column) {
    tgtAssert(index < selection_.size(), "PlotSelection::setColumn(): Index out of bounds");
    selection_[index].first = column;
}

void PlotSelection::setPredicate(size_t index, const PlotPredicate* pred) {
    tgtAssert(index < selection_.size(), "PlotSelection::setPredicate(): Index out of bounds");
    if (selection_[index].second != 0)
        delete selection_[index].second;
    if (pred != 0) {
        PlotPredicate* p = pred->clone();
        selection_[index].second = p;
    }
    else {
        selection_[index].second = 0;
    }
}

std::string PlotSelection::toString(const PlotData& pData) const {
    std::stringstream ss;
    if (isTablePositionFlag_) {
        ss << "Table Position: (" << tablePosition_.x << ", " << tablePosition_.y << ")";
    }
    else {
        for (std::vector< std::pair<int, PlotPredicate*> >::const_iterator it = selection_.begin(); it < selection_.end(); ++it) {
            switch (it->first) {
                case PlotSelectionProperty::X_AXIS_COLUMN:
                    ss << "x Axis: ";
                    break;
                case PlotSelectionProperty::Y_AXIS_COLUMN:
                    ss << "y Axis: ";
                    break;
                case PlotSelectionProperty::Z_AXIS_COLUMN:
                    ss << "z Axis: ";
                    break;
                default:
                    if (pData.getColumnCount() > it->first)
                        ss << pData.getColumnLabel(it->first) << ": ";
                    else
                        ss << it->first << ": ";
                    break;
            }
            if (it->second != 0)
                ss << it->second->toString();
            else
                ss << "No Predicate";
            // add line break if not last predicate
            if (it+1 != selection_.end())
                ss << std::endl;
        }
    }
    return ss.str();
}


bool PlotSelection::operator==(const PlotSelection& rhs) const {
    if (isTablePositionFlag_ != rhs.isTablePositionFlag_)
        return false;
    if (isTablePositionFlag_)
        return (tablePosition_ == rhs.tablePosition_);
    else {
        std::vector< std::pair<int, PlotPredicate*> >::const_iterator it, rhsIt;
        for (it = selection_.begin(), rhsIt = rhs.selection_.begin(); it != selection_.end() && rhsIt != rhs.selection_.end(); ++it, ++rhsIt) {
            if (it->first != rhsIt->first)
                return false;
            if (it->second == 0 && rhsIt->second == 0)
                continue;
            if (it->second == 0 || rhsIt->second == 0)
                return false;
            if (typeid(*(it->second)) != typeid(*(rhsIt->second)))
                return false;
            if (it->second->getThresholdValues() != rhsIt->second->getThresholdValues())
                return false;
        }
    }
    return true;
}

//----------------------------------------------------------------------------

void PlotSelectionEntry::serialize(Serializer& s) const {
    s.serialize("selection", selection_);
    s.serialize("highlight", highlight_);
    s.serialize("renderLabel", renderLabel_);
    s.serialize("zoomTo", zoomTo_);
}

void PlotSelectionEntry::deserialize(Deserializer& d) {
    d.deserialize("selection", selection_);
    d.deserialize("highlight", highlight_);
    d.deserialize("renderLabel", renderLabel_);
    d.deserialize("zoomTo", zoomTo_);
}

PlotSelectionEntry::PlotSelectionEntry() {
}

PlotSelectionEntry::PlotSelectionEntry(PlotSelection selection, bool highlight, bool renderLabel, bool zoomTo)
    : selection_(selection)
    , highlight_(highlight)
    , renderLabel_(renderLabel)
    , zoomTo_(zoomTo)
{}

} // namespace voreen
