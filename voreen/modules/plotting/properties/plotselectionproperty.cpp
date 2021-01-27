/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "plotselectionproperty.h"
#include "../datastructures/plotdata.h"
#include "../datastructures/plotrow.h"

#include "tgt/assert.h"

#include <list>
#include <limits>
#include <algorithm>
#include <sstream>
#include <complex>

namespace voreen {

// PlotSelectionProperty::LabelSelectionIterator-------------------------------

PlotSelectionProperty::LabelSelectionIterator::LabelSelectionIterator(
            const std::vector<PlotSelectionEntry>& selectionEntries)
    : selectionEntries_(selectionEntries)
    , position_(selectionEntries.begin())
{
    while (position_ != selectionEntries_.end() && !position_->renderLabel_)
        ++position_;
}

PlotSelectionProperty::LabelSelectionIterator::LabelSelectionIterator(
            const std::vector<PlotSelectionEntry>& selectionEntries, std::vector<PlotSelectionEntry>::const_iterator position)
    : selectionEntries_(selectionEntries)
    , position_(position)
{}

bool PlotSelectionProperty::LabelSelectionIterator::operator==(LabelSelectionIterator const& rhs) const {
    return (position_ == rhs.position_);
}

bool PlotSelectionProperty::LabelSelectionIterator::operator!=(LabelSelectionIterator const& rhs) const {
    return !(*this == rhs);
}

PlotSelectionProperty::LabelSelectionIterator& PlotSelectionProperty::LabelSelectionIterator::operator++() {
    // check if we've already reached the end
    if (position_ == selectionEntries_.end())
        return *this;
    // otherwise advance iterator until reached entry with renderLabel flag set or end
    do {
        ++position_;
    } while (position_ != selectionEntries_.end() && !position_->renderLabel_);
    return *this;
}

PlotSelectionProperty::LabelSelectionIterator PlotSelectionProperty::LabelSelectionIterator::operator++(int) {
    LabelSelectionIterator tmp (*this);
    ++(*this);
    return tmp;
}

PlotSelectionProperty::LabelSelectionIterator& PlotSelectionProperty::LabelSelectionIterator::operator--() {
    // check if we've already reached the end
    if (position_ != selectionEntries_.begin())
        return *this;
    // otherwise advance iterator backwards until reached entry with renderLabel flag set or end
    do {
        --position_;
    } while (position_ != selectionEntries_.begin() && !position_->renderLabel_);
    return *this;
}

PlotSelectionProperty::LabelSelectionIterator PlotSelectionProperty::LabelSelectionIterator::operator--(int) {
    LabelSelectionIterator tmp (*this);
    --(*this);
    return tmp;
}

const PlotSelection& PlotSelectionProperty::LabelSelectionIterator::operator* () const {
    return position_->selection_;
}

const PlotSelection* PlotSelectionProperty::LabelSelectionIterator::operator-> () const {
    return &(position_->selection_);
}

// PlotSelectionProperty ------------------------------------------------------

const std::string PlotSelectionProperty::loggerCat_("voreen.Properties.PlotSelectionProperty");

PlotSelectionProperty::PlotSelectionProperty(const std::string& id, const std::string& guiText,
                                             PlotData* pData, const PlotEntitiesProperty& entitiesProp,
                                             bool threeDimensional,
                                             Processor::InvalidationLevel invalidationLevel)
    : TemplateProperty< std::vector< PlotSelectionEntry > >(id, guiText, std::vector< PlotSelectionEntry >(), invalidationLevel)
    , plotData_(pData)
    , entitiesProp_(entitiesProp)
    , threeDimensional_(threeDimensional)
    , enableHighlightFlags_(true)
    , enableLabelFlags_(true)
    , enableZoomToFlags_(true)
{
    tgtAssert(pData, "PlotSelectionProperty(): plotData is NULL");
}

PlotSelectionProperty::PlotSelectionProperty()
    : TemplateProperty< std::vector< PlotSelectionEntry > >("", "", std::vector< PlotSelectionEntry >(), Processor::INVALID_RESULT)
    , dummyEntitiesProp_()
    , plotData_(0)
    , entitiesProp_(dummyEntitiesProp_)
    , threeDimensional_(false)
    , enableHighlightFlags_(false)
    , enableLabelFlags_(false)
    , enableZoomToFlags_(false)
{}

Property* PlotSelectionProperty::create() const {
    return new PlotSelectionProperty();
}

void PlotSelectionProperty::add(const PlotSelectionEntry& value, bool additive) {
    // if not additive selection, clear all according flags
    if (value.zoomTo_ && !additive)
        clearZoomFlags();
    if (value.highlight_ && !additive)
        clearHighlightFlags();
    if (value.renderLabel_ && !additive)
        clearLabelFlags();

    // add selection to property
    value_.push_back(value);

    // update internal states and invalidate owners
    if (value.zoomTo_)
        updateZoomState();
    if (value.highlight_)
        applyHighlights();
    invalidate();
}

void PlotSelectionProperty::remove(size_t index) {
    tgtAssert(index < value_.size(), "PlotSelectionProperty::remove(): index out of bounds!");
    value_.erase(value_.begin() + index);
    updateZoomState();
    applyHighlights();
    invalidate();
}

void PlotSelectionProperty::removeUnused() {
    for (size_t i = value_.size(); i > 0; --i) {
        if (!value_[i-1].highlight_ && !value_[i-1].renderLabel_ && !value_[i-1].zoomTo_)
            remove(i-1);
    }
    invalidate();
}

void PlotSelectionProperty::clear() {
    value_.clear();
    updateZoomState();
    applyHighlights();
    invalidate();
}

void PlotSelectionProperty::clearFlags(bool highlight, bool label, bool zoomTo) {
    for (std::vector<PlotSelectionEntry>::iterator it = value_.begin(); it < value_.end(); ++it) {
        if (zoomTo)
            it->zoomTo_ = false;
        if (highlight)
            it->highlight_ = false;
        if (label)
            it->renderLabel_ = false;
    }
    applyHighlights();
    updateZoomState();
    invalidate();
}

void PlotSelectionProperty::clearLastFlag() {
    size_t i = value_.size();
    while (i > 0 && !value_[i-1].zoomTo_ && !value_[i-1].highlight_ && !value_[i-1].renderLabel_)
        --i;
    if (i > 0) {
        value_[i-1].zoomTo_      = false;
        value_[i-1].highlight_   = false;
        value_[i-1].renderLabel_ = false;
        applyHighlights();
        updateZoomState();
        invalidate();
    }
}

void PlotSelectionProperty::clearZoomFlags() {
    for (std::vector<PlotSelectionEntry>::iterator it = value_.begin(); it < value_.end(); ++it)
        it->zoomTo_ = false;
    updateZoomState();
    invalidate();
}

void PlotSelectionProperty::clearLastZoomFlag() {
    size_t i = value_.size();
    while (i > 0 && !value_[i-1].zoomTo_)
        --i;
    if (i > 0) {
        value_[i-1].zoomTo_ = false;
        if (i > 1)
            value_[i-2].zoomTo_ = true;
        updateZoomState();
        invalidate();
    }
}

void PlotSelectionProperty::clearHighlightFlags() {
    for (std::vector<PlotSelectionEntry>::iterator it = value_.begin(); it < value_.end(); ++it)
        it->highlight_ = false;
    applyHighlights();
    invalidate();
}

void PlotSelectionProperty::clearLastHighlightFlag() {
    size_t i = value_.size();
    while (i > 0 && !value_[i-1].highlight_)
        --i;
    if (i > 0) {
        value_[i-1].highlight_ = false;
        if (i > 1)
            value_[i-2].highlight_ = true;
        applyHighlights();
        invalidate();
    }
}

void PlotSelectionProperty::clearLabelFlags() {
    for (std::vector<PlotSelectionEntry>::iterator it = value_.begin(); it < value_.end(); ++it)
        it->renderLabel_ = false;
    invalidate();
}

void PlotSelectionProperty::clearLastLabelFlag() {
    size_t i = value_.size();
    while (i > 0 && !value_[i-1].renderLabel_)
        --i;
    if (i > 0) {
        value_[i-1].renderLabel_ = false;
        if (i > 1)
            value_[i-2].renderLabel_ = true;
        invalidate();
    }
}

const std::vector<PlotSelectionEntry>& PlotSelectionProperty::get() const {
    return value_;
}

const PlotData* PlotSelectionProperty::getPlotData() const {
    return plotData_;
}

void PlotSelectionProperty::setPlotData(PlotData* plotData) {
    tgtAssert(plotData, "PlotSelectionProperty::setPlotData(): plotData is NULL");
    plotData_ = plotData;
    applyHighlights();
    updateZoomState();
}

const PlotEntitiesProperty& PlotSelectionProperty::getEntitiesProperty() const {
    return entitiesProp_;
}

bool PlotSelectionProperty::getThreeDimensional() const {
    return threeDimensional_;
}

const PlotSelectionEntry& PlotSelectionProperty::getSelectionAt(size_t index) {
    tgtAssert(index < value_.size(), "PlotSelectionProperty::getSelectionAt(): index out of bounds!");
    return value_[index];
}

void PlotSelectionProperty::setDimension(bool threeDimensional) {
    if (threeDimensional_ != threeDimensional) {
        threeDimensional_ = threeDimensional;
        clear();
    }
}

void PlotSelectionProperty::setBaseZoomState(PlotZoomState state) {
    baseZoomState_ = state;
    updateZoomState();
}

void PlotSelectionProperty::set(const std::vector<PlotSelectionEntry>& value) {
    clear();
    value_ = value;
    applyHighlights();
    updateZoomState();
    invalidate();
}

void PlotSelectionProperty::set(PlotSelectionEntry entry, int index) {
    if (static_cast<size_t>(index) >= value_.size())
        value_.push_back(entry);
    else {
        value_.erase(value_.begin()+index);
        value_.insert(value_.begin()+index, entry);
    }
    if (entry.highlight_)
        applyHighlights();
    if (entry.zoomTo_)
        updateZoomState();
    invalidate();
}

int PlotSelectionProperty::size() const {
    return static_cast<int>(value_.size());
}

void PlotSelectionProperty::setHighlight(size_t index, bool value) {
    tgtAssert(index < value_.size(), "PlotSelectionProperty::setHighlight(): index out of bounds!");
    value_[index].highlight_ = value;
    applyHighlights();
    invalidate();
}

void PlotSelectionProperty::setRenderLabel(size_t index, bool value) {
    tgtAssert(index < value_.size(), "PlotSelectionProperty::setRenderLabel(): index out of bounds!");
    value_[index].renderLabel_ = value;
    invalidate();
}

void PlotSelectionProperty::setZoomTo(size_t index, bool value) {
    tgtAssert(index < value_.size(), "PlotSelectionProperty::setZoomTo(): index out of bounds!");
    value_[index].zoomTo_ = value;
    updateZoomState();
    invalidate();
}

bool PlotSelectionProperty::getEnableHighlightFlags() const {
    return enableHighlightFlags_;
}

void PlotSelectionProperty::setEnableHighlightFlags(bool value) {
    enableHighlightFlags_ = value;
}

bool PlotSelectionProperty::getEnableLabelFlags() const {
    return enableLabelFlags_;
}

void PlotSelectionProperty::setEnableLabelFlags(bool value) {
    enableLabelFlags_ = value;
}

bool PlotSelectionProperty::getEnableZoomToFlags() const {
    return enableZoomToFlags_;
}

void PlotSelectionProperty::setEnableZoomToFlags(bool value) {
    enableZoomToFlags_ = value;
}

const PlotZoomState& PlotSelectionProperty::getZoom() const {
    return currentZoomState_;
}

void PlotSelectionProperty::applyHighlights() {
    if (!entitiesProp_.dataValid() || !enableHighlightFlags_)
        return;

    // clear highlights first
    plotData_->clearHighlights();

    // If any of the selections has PlotPredicates applied to one of our special column indices
    // (e.g. X_AXIS_COLUMN) we need to decode that special column to a real column index using
    // our entitiesProp_ member reference. This leads to constructing a new vector of
    // PlotSelectionEntries but also quite simplifies the loop the predicates run
    // through later...

    // This vector of vectors of predicates is to understand as follows:
    //  - The predicates in the inner vector are applied as conjunction (each has to match)
    //  - The predicates (more precisely the results) in the outer vector are applied as disjunction
    //    (just one of them has to match)
    std::vector< std::vector< std::pair<int, const PlotPredicate* > > > predicateSelections;

    for (std::vector<PlotSelectionEntry>::const_iterator pse = value_.begin(); pse < value_.end(); ++pse) {
        if (pse->highlight_) {
            if (! pse->selection_.isTablePosition()) {
                // check if we have any Y_AXIS_COLUMN/Z_AXIS_COLUMN predicates
                std::vector<const PlotPredicate*> dataColumnPredicates;
                std::vector< std::pair<int, const PlotPredicate* > > conjunctionSelections;
                std::vector< std::pair<int, PlotPredicate* > >::const_iterator pair;
                for (pair = pse->selection_.getSelection().begin();
                            pair < pse->selection_.getSelection().end(); ++pair) {
                    // we got a valid column index
                    if (pair->first >= 0 && pair->first < plotData_->getColumnCount()) {
                        conjunctionSelections.push_back(*pair);
                    }
                    // we've got the special 'x axis colum' index - add predicate to each according column
                    else if (pair->first == X_AXIS_COLUMN) {
                        conjunctionSelections.push_back(std::make_pair(
                                    entitiesProp_.getXColumnIndex(), pair->second));
                    }
                    // we've got the special 'y axis colum' index and it is not a data column - add
                    // predicate to each according column
                    else if (threeDimensional_ && pair->first == Y_AXIS_COLUMN) {
                        conjunctionSelections.push_back(std::make_pair(
                                    entitiesProp_.getYColumnIndex(), pair->second));
                    }
                    // we've got a special data column index - add the predicate to dataColumnPredicates
                    // list as for each of them an outer predicate will be created
                    else if (   (!threeDimensional_ && pair->first == Y_AXIS_COLUMN)
                        || (threeDimensional_ && pair->first == Z_AXIS_COLUMN))
                        dataColumnPredicates.push_back(pair->second);
                }

                std::vector<int> columns = entitiesProp_.getDataColumnIndices();
                // if we have no data columns or no predicates for them one outer predicate is enough
                if (columns.empty() || dataColumnPredicates.empty()) {
                    predicateSelections.push_back(conjunctionSelections);
                }
                // otherwise we need to create an outer predicate for each data column predicate
                else {
                    for (std::vector<int>::const_iterator col = columns.begin(); col < columns.end(); ++col) {
                        predicateSelections.push_back(conjunctionSelections);
                        std::vector<const PlotPredicate*>::iterator pp;
                        for (pp = dataColumnPredicates.begin(); pp != dataColumnPredicates.end(); ++pp) {
                            predicateSelections.back().push_back(std::make_pair(*col, *pp));
                        }
                    }
                }
            }
            // on our way we already can highlight the cells specified by table positions:
            else {
                plotData_->setHighlight(pse->selection_.getTablePosition(), true, true);
            }
        }
    }

    std::vector< std::vector< std::pair<int, const PlotPredicate* > > >::iterator outerIt;
    std::vector< std::pair<int, const PlotPredicate* > >::iterator innerIt;

    // walk through PlotData and check for each row if one of the outer predicates match
    for (int row = 0; row < plotData_->getRowsCount(); ++row) {
        const PlotRowValue& plotRow = plotData_->getRow(row);

        // now walk through our previously generated predicates
        for (outerIt = predicateSelections.begin(); outerIt != predicateSelections.end(); ++outerIt) {
            // aggregated matched state
            bool matched = true;

            // check each inner predicate
            for (innerIt = outerIt->begin(); matched && innerIt != outerIt->end(); ++innerIt) {
                if (innerIt->second != 0)
                    matched &= innerIt->second->check(plotRow.getCellAt(innerIt->first));
            }

            // if all predicates matched => highlight all cells referenced in inner predicates
            if (matched) {
                for (innerIt = outerIt->begin(); innerIt != outerIt->end(); ++innerIt)
                    plotData_->setHighlight(tgt::ivec2(row, innerIt->first), true, true);
            }
        }
    }
}

void PlotSelectionProperty::updateZoomState() {
    if (!enableZoomToFlags_)
        return;

    // reset intervals first
    currentZoomState_.xZoom_ = Interval<plot_t>();
    currentZoomState_.yZoom_ = Interval<plot_t>();
    currentZoomState_.zZoom_ = Interval<plot_t>();

    // get the column indices of all data columns
    std::vector<int> dataColumns = entitiesProp_.getDataColumnIndices();

    // flag whether to use baseZoomState_
    bool useBaseZoom = true;

    // now obtain new values out of all PlotSelectionEntries with zomm flag set
    for (std::vector<PlotSelectionEntry>::const_iterator it = value_.begin(); it < value_.end(); ++it) {
        if (it->zoomTo_) {
            // we found a zoom selection so we do not use baseZoomState_
            useBaseZoom = false;

            // handle predicate selections
            if (! it->selection_.isTablePosition()) {
                // for each predicate
                for (std::vector< std::pair<int, PlotPredicate* > >::const_iterator pit = it->selection_.getSelection().begin(); pit < it->selection_.getSelection().end(); ++pit) {
                    // no predicate no zoom
                    if (pit->second == 0)
                        continue;

                    // check if associated column is special axis column constant or in axis columns vector
                    // if so union zoom interval with predicate interval
                    if (pit->first == X_AXIS_COLUMN || pit->first == entitiesProp_.getXColumnIndex() || PlotData::isIndexColumn(*plotData_, pit->first))
                        currentZoomState_.xZoom_.unionWith(pit->second->getIntervalRepresentation());
                    else if (  (pit->first == Y_AXIS_COLUMN)
                            || (threeDimensional_ && pit->first == entitiesProp_.getYColumnIndex())
                            || (!threeDimensional_ && std::find(dataColumns.begin(), dataColumns.end(), pit->first) != dataColumns.end()))
                        currentZoomState_.yZoom_.unionWith(pit->second->getIntervalRepresentation());
                    else if (pit->first == Z_AXIS_COLUMN || (threeDimensional_ && std::find(dataColumns.begin(), dataColumns.end(), pit->first) != dataColumns.end()))
                        currentZoomState_.zZoom_.unionWith(pit->second->getIntervalRepresentation());
                }
            }
            // handle table position selections
            else if (entitiesProp_.dataValid()) {
                tgt::ivec2 pos = it->selection_.getTablePosition();
                if (pos.x >= -1 && pos.x < plotData_->getRowsCount() && pos.y >= 0 && pos.y < plotData_->getColumnCount()) {
                    // if pos.x == -1 the entire column is selected
                    int start = 0;
                    int end = plotData_->getRowsCount();
                    if (pos.x != -1) {
                        start = pos.x;
                        end = pos.x+1;
                    }
                    for (int i = start; i < end; ++i) {
                        const PlotRowValue& row = plotData_->getRow(i);
                        // check if table position encodes a data column
                        if (std::find(dataColumns.begin(), dataColumns.end(), pos.y) != dataColumns.end()) {
                            currentZoomState_.xZoom_.nibble(row.getValueAt(entitiesProp_.getXColumnIndex()));
                            if (! threeDimensional_) {
                                currentZoomState_.yZoom_.nibble(row.getValueAt(pos.y));
                            }
                            else {
                                currentZoomState_.yZoom_.nibble(row.getValueAt(entitiesProp_.getYColumnIndex()));
                                currentZoomState_.zZoom_.nibble(row.getValueAt(pos.y));
                            }
                        }
                        else if (pos.y == entitiesProp_.getXColumnIndex())
                            currentZoomState_.xZoom_.nibble(row.getValueAt(pos.y));
                        else if (threeDimensional_ && pos.y == entitiesProp_.getYColumnIndex())
                            currentZoomState_.yZoom_.nibble(row.getValueAt(pos.y));
                    }
                }
            }
        }
    }
    // clamp/enlarge currentZoomState to baseZoomState
    currentZoomState_.xZoom_.clampInfinitesTo(baseZoomState_.xZoom_);
    currentZoomState_.yZoom_.clampInfinitesTo(baseZoomState_.yZoom_);
    currentZoomState_.zZoom_.clampInfinitesTo(baseZoomState_.zZoom_);

    if (useBaseZoom)
        currentZoomState_ = baseZoomState_;
}

PlotSelectionProperty::LabelSelectionIterator PlotSelectionProperty::getLabelsBegin() const {
    return LabelSelectionIterator(value_);
}

PlotSelectionProperty::LabelSelectionIterator PlotSelectionProperty::getLabelsEnd() const {
    return LabelSelectionIterator(value_, value_.end());
}

void PlotSelectionProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("values", value_);
    s.serialize("threeDimensional", threeDimensional_);
}

void PlotSelectionProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    s.deserialize("values", value_);
    s.deserialize("threeDimensional", threeDimensional_);
}

}   // namespace
