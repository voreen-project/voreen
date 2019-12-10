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

#include "plotcell.h"
#include "plotdata.h"
#include "aggregationfunction.h"

#include "voreen/core/io/serialization/serialization.h"

#include <algorithm>
#include <limits>

namespace voreen {

// PlotCellValue methods: ---------------------------------

PlotCellValue::PlotCellValue()
    : value_(0)
    , isValue_(false)
    , isTag_(false)
    , isHighlighted_(false)
{}

PlotCellValue::PlotCellValue(const PlotCellValue& rhs)
    : value_(rhs.value_)
    , isValue_(rhs.isValue_)
    , isTag_(rhs.isTag_)
    , isHighlighted_(rhs.isHighlighted_)
{
    if (isTag_)
        tag_ = new std::string(*rhs.tag_);
    else if (! isValue_) // no tag, no value -> null
        value_ = 0; // set value to 0 for safety reasons
}

PlotCellValue::PlotCellValue(plot_t value)
    : value_(value)
    , isValue_(true)
    , isTag_(false)
    , isHighlighted_(false)
{
}

PlotCellValue::PlotCellValue(const std::string& tag)
    : tag_(new std::string(tag))
    , isValue_(false)
    , isTag_(true)
    , isHighlighted_(false)
{
}

PlotCellValue::~PlotCellValue() {
    if (isTag())
        delete tag_;
}

PlotCellValue& PlotCellValue::operator=(PlotCellValue rhs) {
    std::swap(isValue_, rhs.isValue_);
    std::swap(isTag_, rhs.isTag_);
    std::swap(isHighlighted_, rhs.isHighlighted_);
    if (isValue_) {
        std::swap(value_, rhs.value_);
        return *this;
    }
    else if (isTag_) {
        std::swap(tag_, rhs.tag_); // no worries: old tag will deleted upon function exit (~rhs)
        return *this;
    }
    else { // isNull
        value_ = 0;
        return *this;
    }
}

bool PlotCellValue::operator==(const PlotCellValue& rhs) const {
    if (isValue()) {
        return rhs.isValue() && (getValue() == rhs.getValue())
            && (isHighlighted() == rhs.isHighlighted());
    }
    else {
        return rhs.isTag() && (getTag() == rhs.getTag())
            && (isHighlighted() == rhs.isHighlighted());
    }
}

bool PlotCellValue::operator<(const PlotCellValue& rhs) const {
    if (isValue()) {
        return rhs.isTag()
            || getValue() < rhs.getValue()
            || (getValue() == rhs.getValue() && isHighlighted() < rhs.isHighlighted());
    }
    else {
        return rhs.isTag()
            && (getTag() < rhs.getTag() || (getTag() == rhs.getTag() && isHighlighted() < rhs.isHighlighted()));
    }
}

bool PlotCellValue::operator>(const PlotCellValue& rhs) const {
    return !(*this < rhs || *this == rhs);
}

std::ostream& operator<<(std::ostream& ostr, const PlotCellValue& rhs) {
    std::stringstream ss;
    if (rhs.isValue())
        ss << " Value= " << rhs.getValue();
    else if (rhs.isTag())
        ss << " Tag = " << rhs.getTag();
    ostr << ss.str();
    return ostr;
}

plot_t PlotCellValue::getValue() const {
    return (isValue() ? value_ : std::numeric_limits<plot_t>::quiet_NaN());
}

std::string PlotCellValue::getTag() const {
    return (isTag() ? *tag_ : "");
}


void PlotCellValue::setValue(plot_t value) {
    if (isTag() && tag_)
        delete tag_;

    value_ = value;
    isValue_ = true;
    isTag_ = false;
}

void PlotCellValue::setTag(const std::string& tag) {
    if (isTag() && tag_)
        delete tag_;

    tag_ = new std::string(tag);
    isValue_ = false;
    isTag_ = true;
}

void PlotCellValue::setHighlighted(bool highlightedState) {
    isHighlighted_ = highlightedState;
}


void PlotCellValue::clear() {
    if (isTag() && tag_)
        delete tag_;

    tag_ = 0;
    isValue_ = false;
    isTag_ = false;
    isHighlighted_ = false;
}

void PlotCellValue::serialize(Serializer& s) const {
    s.serialize("isValue", isValue());
    s.serialize("isTag", isTag());
    s.serialize("isHighlighted", isHighlighted());

    if (isTag())
        s.serialize("tag", getTag());
    if (isValue())
        s.serialize("value", getValue());
}

void PlotCellValue::deserialize(Deserializer& d) {

    d.deserialize("isValue", isValue_);
    d.deserialize("isTag", isTag_);

    // now call the according constructors
    if (isValue_) {
        d.deserialize("value", value_);
//        deserializeSimpleTypes();
//        data = PlotCellValue(value);
    }
    else if (isTag_) {
        d.deserialize("tag", tag_);
//        deserializeSimpleTypes();
//        data = PlotCellValue(tag);
    }

    // finally set highlighted flag
//    data.setHighlighted(isHighlighted);
    d.deserialize("isHighlighted", isHighlighted_);
}

// PlotCellImplicit methods: ---------------------------------

PlotCellImplicit::PlotCellImplicit(AggregationFunction* func, const PlotData* data, int column)
    : plotData_(data)
    , column_(column)
{
    if (func != NULL) {
        aggregationFunction_ = func->clone();
    }
    else
        aggregationFunction_ = NULL;
}

PlotCellImplicit::~PlotCellImplicit() {
    // do not delete pointers as we do not own them and others might want to use them later
}

AggregationFunction* PlotCellImplicit::getAggregationFunction() const {
    return aggregationFunction_;
}

AggregationFunction* PlotCellImplicit::getAggregationFunctionClone() const {
    if (aggregationFunction_)
        return aggregationFunction_->clone();
    else
        return NULL;
}

plot_t PlotCellImplicit::getValue() const {
    if (!plotData_ || column_ < 0 || column_ >= plotData_->getColumnCount()) {
        tgtAssert(false, "PlotCellImplicit::getValue: column out of bounds.");
        return 0;
    }
    else if (aggregationFunction_)
        return plotData_->aggregate(column_, aggregationFunction_);
    else
        return 0;
}

bool PlotCellImplicit::isNull() const {
    return (plotData_ == 0 || column_ == 0);
}

void PlotCellImplicit::clear() {
    plotData_ = NULL;
    column_ = 0;
}

}
