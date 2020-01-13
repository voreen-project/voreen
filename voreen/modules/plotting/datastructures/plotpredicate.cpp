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

#include "plotpredicate.h"

#include <limits>
#include <sstream>
#include <iomanip>
#include <string>

namespace voreen {

namespace {
    void setSmartPrecision(std::stringstream& ss, size_t precision, plot_t value) {
        double base = log10(fabs(value));
        if (base > 0)
            ss << std::setprecision(static_cast<int>(ceil(base) + precision));
        else if (base < 0)
            ss << std::setprecision(precision);
    }
}

// PlotPredicateLess methods -------------------------------------------------------

PlotPredicateLess::PlotPredicateLess()
    : threshold_(PlotCellValue())
    {}

PlotPredicateLess::PlotPredicateLess(plot_t threshold)
    : threshold_(PlotCellValue(threshold))
    {}
PlotPredicateLess::PlotPredicateLess(std::string threshold)
    : threshold_(PlotCellValue(threshold))
    {}

PlotPredicateLess::PlotPredicateLess(const PlotCellValue& threshold)
    : threshold_(threshold.isTag() ? PlotCellValue(threshold.getTag()) : PlotCellValue(threshold.getValue()))
    {}

int PlotPredicateLess::getNumberOfThresholdValues() const {
    return 1;
}

void PlotPredicateLess::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (! thresholds.empty())
        threshold_ = thresholds[0];
}

std::vector<PlotCellValue> PlotPredicateLess::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(threshold_);
    return erg;
}

std::vector<std::string> PlotPredicateLess::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Upper Threshold");
    return erg;
}

bool PlotPredicateLess::check(const PlotCellValue& value) const {
    return ((value.isValue() && threshold_.isValue() && value.getValue() <  threshold_.getValue())
        || (value.isTag() && threshold_.isTag() && value.getTag() < threshold_.getTag()));
}

Interval<plot_t> PlotPredicateLess::getIntervalRepresentation() const {
    return Interval<plot_t>(-std::numeric_limits<plot_t>::max(), threshold_.getValue(), false, true);
}

PlotPredicate* PlotPredicateLess::clone() const {
    return new PlotPredicateLess(threshold_);
}

std::string PlotPredicateLess::toString() const {
    std::stringstream ss;
    ss << "< ";
    if (threshold_.isNull())
        ss << "NULL";
    else if (threshold_.isValue()) {
        setSmartPrecision(ss, 3, threshold_.getValue());
        ss << threshold_.getValue();
    }
    else
        ss << "\"" << threshold_.getTag() << "\"";
    return ss.str();
}

void PlotPredicateLess::serialize(Serializer& s) const {
    s.serialize("Threshold",threshold_);
}

void PlotPredicateLess::deserialize(Deserializer& s) {
    s.deserialize("Threshold",threshold_);
}

// PlotPredicateEqual methods -------------------------------------------------------

PlotPredicateEqual::PlotPredicateEqual()
    : threshold_(PlotCellValue())
    {}

PlotPredicateEqual::PlotPredicateEqual(plot_t threshold)
    : threshold_(PlotCellValue(threshold))
    {}

PlotPredicateEqual::PlotPredicateEqual(std::string threshold)
    : threshold_(PlotCellValue(threshold))
    {}

PlotPredicateEqual::PlotPredicateEqual(const PlotCellValue& threshold)
    : threshold_(threshold.isTag() ? PlotCellValue(threshold.getTag()) : PlotCellValue(threshold.getValue()))
    {}

int PlotPredicateEqual::getNumberOfThresholdValues() const {
    return 1;
}

void PlotPredicateEqual::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (! thresholds.empty())
        threshold_ = thresholds[0];
}

std::vector<PlotCellValue> PlotPredicateEqual::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(threshold_);
    return erg;
}

std::vector<std::string> PlotPredicateEqual::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Value");
    return erg;
}

bool PlotPredicateEqual::check(const PlotCellValue& value) const {
    return ((value.isValue() && threshold_.isValue() && value.getValue() == threshold_.getValue()) ||
        (value.isTag() && threshold_.isTag() && value.getTag() == threshold_.getTag()));
}

Interval<plot_t> PlotPredicateEqual::getIntervalRepresentation() const {
    return Interval<plot_t>(threshold_.getValue(), threshold_.getValue(), false, false);
}

PlotPredicate* PlotPredicateEqual::clone() const {
    return new PlotPredicateEqual(threshold_);
}

std::string PlotPredicateEqual::toString() const {
    std::stringstream ss;
    ss << "= ";
    if (threshold_.isNull())
        ss << "NULL";
    else if (threshold_.isValue()) {
        setSmartPrecision(ss, 3, threshold_.getValue());
        ss << threshold_.getValue();
    }
    else
        ss << "\"" << threshold_.getTag() << "\"";
    return ss.str();
}

void PlotPredicateEqual::serialize(Serializer& s) const {
    s.serialize("Threshold",threshold_);
}

void PlotPredicateEqual::deserialize(Deserializer& s) {
    s.deserialize("Threshold",threshold_);
}

// PlotPredicateGreater methods -------------------------------------------------------

PlotPredicateGreater::PlotPredicateGreater()
    : threshold_(PlotCellValue())
    {}

PlotPredicateGreater::PlotPredicateGreater(plot_t threshold)
    : threshold_(PlotCellValue(threshold))
    {}

PlotPredicateGreater::PlotPredicateGreater(std::string threshold)
    : threshold_(PlotCellValue(threshold))
    {}

PlotPredicateGreater::PlotPredicateGreater(const PlotCellValue& threshold)
    : threshold_(threshold.isTag() ? PlotCellValue(threshold.getTag()) : PlotCellValue(threshold.getValue()))
    {}

int PlotPredicateGreater::getNumberOfThresholdValues() const {
    return 1;
}

void PlotPredicateGreater::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (! thresholds.empty())
        threshold_ = thresholds[0];
}

std::vector<PlotCellValue> PlotPredicateGreater::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(threshold_);
    return erg;
}

std::vector<std::string> PlotPredicateGreater::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Upper Threshold");
    return erg;
}

bool PlotPredicateGreater::check(const PlotCellValue& value) const {
    return ((value.isValue() && threshold_.isValue() && value.getValue() > threshold_.getValue()) ||
        (value.isTag() && threshold_.isTag() && value.getTag() > threshold_.getTag()));
}

Interval<plot_t> PlotPredicateGreater::getIntervalRepresentation() const {
    return Interval<plot_t>(threshold_.getValue(), std::numeric_limits<plot_t>::max(), true, false);
}

PlotPredicate* PlotPredicateGreater::clone() const {
    return new PlotPredicateGreater(threshold_);
}


std::string PlotPredicateGreater::toString() const {
    std::stringstream ss;
    ss << "> ";
    if (threshold_.isNull())
        ss << "NULL";
    else if (threshold_.isValue()) {
        setSmartPrecision(ss, 3, threshold_.getValue());
        ss << threshold_.getValue();
    }
    else
        ss << "\"" << threshold_.getTag() << "\"";
    return ss.str();
}

void PlotPredicateGreater::serialize(Serializer& s) const {
    s.serialize("Threshold",threshold_);
}

void PlotPredicateGreater::deserialize(Deserializer& s) {
    s.deserialize("Threshold",threshold_);
}

// PlotPredicateBetween methods -------------------------------------------------------

PlotPredicateBetween::PlotPredicateBetween()
    : lowerThreshold_(PlotCellValue())
    , upperThreshold_(PlotCellValue())
    {}

PlotPredicateBetween::PlotPredicateBetween(plot_t lowerThreshold, plot_t upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
    {}

PlotPredicateBetween::PlotPredicateBetween(std::string lowerThreshold, std::string upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
    {}

PlotPredicateBetween::PlotPredicateBetween(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold)
    : lowerThreshold_(lowerThreshold.isTag() ? PlotCellValue(lowerThreshold.getTag()) : PlotCellValue(lowerThreshold.getValue()))
    , upperThreshold_(upperThreshold.isTag() ? PlotCellValue(upperThreshold.getTag()) : PlotCellValue(upperThreshold.getValue()))
    {}

int PlotPredicateBetween::getNumberOfThresholdValues() const {
    return 2;
}

void PlotPredicateBetween::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (thresholds.size() == 2) {
        lowerThreshold_ = thresholds[0];
        upperThreshold_ = thresholds[1];
    }
}

std::vector<PlotCellValue> PlotPredicateBetween::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(lowerThreshold_);
    erg.push_back(upperThreshold_);
    return erg;
}

std::vector<std::string> PlotPredicateBetween::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Lower Threshold");
    erg.push_back("Upper Threshold");
    return erg;
}

bool PlotPredicateBetween::check(const PlotCellValue& value) const {
    return ((value.isValue() && lowerThreshold_.isValue() && upperThreshold_.isValue() &&
        value.getValue() > lowerThreshold_.getValue() && value.getValue() < upperThreshold_.getValue()) ||
        (value.isTag() && lowerThreshold_.isTag() && upperThreshold_.isTag() &&
        value.getTag() > lowerThreshold_.getTag() && value.getTag() < upperThreshold_.getTag()));
}

Interval<plot_t> PlotPredicateBetween::getIntervalRepresentation() const {
    return Interval<plot_t>(lowerThreshold_.getValue(), upperThreshold_.getValue(), true, true);
}

PlotPredicate* PlotPredicateBetween::clone() const {
    return new PlotPredicateBetween(lowerThreshold_,upperThreshold_);
}

std::string PlotPredicateBetween::toString() const {
    std::stringstream ss;
    ss << "in (";
    if (lowerThreshold_.isNull())
        ss << "NULL";
    else if (lowerThreshold_.isValue()) {
        setSmartPrecision(ss, 3, lowerThreshold_.getValue());
        ss << lowerThreshold_.getValue();
    }
    else
        ss << "\"" << lowerThreshold_.getTag() << "\"";
    ss << ", ";
    if (upperThreshold_.isNull())
        ss << "NULL";
    else if (upperThreshold_.isValue()) {
        setSmartPrecision(ss, 3, upperThreshold_.getValue());
        ss << upperThreshold_.getValue();
    }
    else
        ss << "\"" << upperThreshold_.getTag() << "\"";
    ss << ")";
    return ss.str();
}

void PlotPredicateBetween::serialize(Serializer& s) const {
    s.serialize("lowerThreshold",lowerThreshold_);
    s.serialize("upperThreshold",upperThreshold_);
}

void PlotPredicateBetween::deserialize(Deserializer& s) {
    s.deserialize("lowerThreshold",lowerThreshold_);
    s.deserialize("upperThreshold",upperThreshold_);
}

// PlotPredicateNotBetween methods -------------------------------------------------------

PlotPredicateNotBetween::PlotPredicateNotBetween()
    : lowerThreshold_(PlotCellValue())
    , upperThreshold_(PlotCellValue())
    {}

PlotPredicateNotBetween::PlotPredicateNotBetween(plot_t lowerThreshold, plot_t upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
    {}

PlotPredicateNotBetween::PlotPredicateNotBetween(std::string lowerThreshold, std::string upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
    {}

PlotPredicateNotBetween::PlotPredicateNotBetween(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold)
    : lowerThreshold_(lowerThreshold.isTag() ? PlotCellValue(lowerThreshold.getTag()) : PlotCellValue(lowerThreshold.getValue()))
    , upperThreshold_(upperThreshold.isTag() ? PlotCellValue(upperThreshold.getTag()) : PlotCellValue(upperThreshold.getValue()))
    {}

int PlotPredicateNotBetween::getNumberOfThresholdValues() const {
    return 2;
}

void PlotPredicateNotBetween::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (thresholds.size() == 2) {
        lowerThreshold_ = thresholds[0];
        upperThreshold_ = thresholds[1];
    }
}

std::vector<PlotCellValue> PlotPredicateNotBetween::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(lowerThreshold_);
    erg.push_back(upperThreshold_);
    return erg;
}

std::vector<std::string> PlotPredicateNotBetween::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Lower Threshold");
    erg.push_back("Upper Threshold");
    return erg;
}

bool PlotPredicateNotBetween::check(const PlotCellValue& value) const {
    return ((value.isValue() && lowerThreshold_.isValue() && upperThreshold_.isValue() &&
        (value.getValue() <= lowerThreshold_.getValue() || value.getValue() >= upperThreshold_.getValue())) ||
        (value.isTag() && lowerThreshold_.isTag() && upperThreshold_.isTag() &&
        (value.getTag() <= lowerThreshold_.getTag() || value.getTag() >= upperThreshold_.getTag())));
}

Interval<plot_t> PlotPredicateNotBetween::getIntervalRepresentation() const {
    return Interval<plot_t>(upperThreshold_.getValue(),lowerThreshold_.getValue(), false, false);
}

PlotPredicate* PlotPredicateNotBetween::clone() const {
    return new PlotPredicateNotBetween(lowerThreshold_,upperThreshold_);
}

std::string PlotPredicateNotBetween::toString() const {
    std::stringstream ss;
    ss << "not in (";
    if (lowerThreshold_.isNull())
        ss << "NULL";
    else if (lowerThreshold_.isValue()) {
        setSmartPrecision(ss, 3, lowerThreshold_.getValue());
        ss << lowerThreshold_.getValue();
    }
    else
        ss << "\"" << lowerThreshold_.getTag() << "\"";
    ss << ", ";
    if (upperThreshold_.isNull())
        ss << "NULL";
    else if (upperThreshold_.isValue()) {
        setSmartPrecision(ss, 3, upperThreshold_.getValue());
        ss << upperThreshold_.getValue();
    }
    else
        ss << "\"" << upperThreshold_.getTag() << "\"";
    ss << ")";
    return ss.str();
}

void PlotPredicateNotBetween::serialize(Serializer& s) const {
    s.serialize("lowerThreshold",lowerThreshold_);
    s.serialize("upperThreshold",upperThreshold_);
}

void PlotPredicateNotBetween::deserialize(Deserializer& s) {
    s.deserialize("lowerThreshold",lowerThreshold_);
    s.deserialize("upperThreshold",upperThreshold_);
}

// PlotPredicateBetweenOrEqual methods -------------------------------------------------------

PlotPredicateBetweenOrEqual::PlotPredicateBetweenOrEqual()
    : lowerThreshold_(PlotCellValue())
    , upperThreshold_(PlotCellValue())
{}

PlotPredicateBetweenOrEqual::PlotPredicateBetweenOrEqual(plot_t lowerThreshold, plot_t upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
{}

PlotPredicateBetweenOrEqual::PlotPredicateBetweenOrEqual(std::string lowerThreshold, std::string upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
{}

PlotPredicateBetweenOrEqual::PlotPredicateBetweenOrEqual(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold)
    : lowerThreshold_(lowerThreshold.isTag() ? PlotCellValue(lowerThreshold.getTag()) : PlotCellValue(lowerThreshold.getValue()))
    , upperThreshold_(upperThreshold.isTag() ? PlotCellValue(upperThreshold.getTag()) : PlotCellValue(upperThreshold.getValue()))
{}

int PlotPredicateBetweenOrEqual::getNumberOfThresholdValues() const {
    return 2;
}

void PlotPredicateBetweenOrEqual::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (thresholds.size() == 2) {
        lowerThreshold_ = thresholds[0];
        upperThreshold_ = thresholds[1];
    }
}

std::vector<PlotCellValue> PlotPredicateBetweenOrEqual::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(lowerThreshold_);
    erg.push_back(upperThreshold_);
    return erg;
}

std::vector<std::string> PlotPredicateBetweenOrEqual::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Lower Threshold");
    erg.push_back("Upper Threshold");
    return erg;
}

bool PlotPredicateBetweenOrEqual::check(const PlotCellValue& value) const {
    if (value.isValue()) {
        return (   lowerThreshold_.isValue()
                && upperThreshold_.isValue()
                && value.getValue() >= lowerThreshold_.getValue()
                && value.getValue() <= upperThreshold_.getValue());
    }
    else if (value.isTag()) {
        return (   lowerThreshold_.isTag()
                && upperThreshold_.isTag()
                && value.getTag() >= lowerThreshold_.getTag()
                && value.getTag() <= upperThreshold_.getTag());
    }
    return false;
}

Interval<plot_t> PlotPredicateBetweenOrEqual::getIntervalRepresentation() const {
    return Interval<plot_t>(lowerThreshold_.getValue(), upperThreshold_.getValue(), false, false);
}

PlotPredicate* PlotPredicateBetweenOrEqual::clone() const {
    return new PlotPredicateBetweenOrEqual(lowerThreshold_,upperThreshold_);
}

std::string PlotPredicateBetweenOrEqual::toString() const {
    std::stringstream ss;
    ss << "in [";
    if (lowerThreshold_.isNull())
        ss << "NULL";
    else if (lowerThreshold_.isValue()) {
        setSmartPrecision(ss, 3, lowerThreshold_.getValue());
        ss << lowerThreshold_.getValue();
    }
    else
        ss << "\"" << lowerThreshold_.getTag() << "\"";
    ss << ", ";
    if (upperThreshold_.isNull())
        ss << "NULL";
    else if (upperThreshold_.isValue()) {
        setSmartPrecision(ss, 3, upperThreshold_.getValue());
        ss << upperThreshold_.getValue();
    }
    else
        ss << "\"" << upperThreshold_.getTag() << "\"";
    ss << "]";
    return ss.str();
}

void PlotPredicateBetweenOrEqual::serialize(Serializer& s) const {
    s.serialize("lowerThreshold",lowerThreshold_);
    s.serialize("upperThreshold",upperThreshold_);
}

void PlotPredicateBetweenOrEqual::deserialize(Deserializer& s) {
    s.deserialize("lowerThreshold",lowerThreshold_);
    s.deserialize("upperThreshold",upperThreshold_);
}

// PlotPredicateNotBetweenOrEqual methods -------------------------------------------------------

PlotPredicateNotBetweenOrEqual::PlotPredicateNotBetweenOrEqual()
    : lowerThreshold_(PlotCellValue())
    , upperThreshold_(PlotCellValue())
{}

PlotPredicateNotBetweenOrEqual::PlotPredicateNotBetweenOrEqual(plot_t lowerThreshold, plot_t upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
{}

PlotPredicateNotBetweenOrEqual::PlotPredicateNotBetweenOrEqual(std::string lowerThreshold, std::string upperThreshold)
    : lowerThreshold_(PlotCellValue(lowerThreshold))
    , upperThreshold_(PlotCellValue(upperThreshold))
{}

PlotPredicateNotBetweenOrEqual::PlotPredicateNotBetweenOrEqual(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold)
    : lowerThreshold_(lowerThreshold.isTag() ? PlotCellValue(lowerThreshold.getTag()) : PlotCellValue(lowerThreshold.getValue()))
    , upperThreshold_(upperThreshold.isTag() ? PlotCellValue(upperThreshold.getTag()) : PlotCellValue(upperThreshold.getValue()))
{}

int PlotPredicateNotBetweenOrEqual::getNumberOfThresholdValues() const {
    return 2;
}

void PlotPredicateNotBetweenOrEqual::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (thresholds.size() == 2) {
        lowerThreshold_ = thresholds[0];
        upperThreshold_ = thresholds[1];
    }
}

std::vector<PlotCellValue> PlotPredicateNotBetweenOrEqual::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(lowerThreshold_);
    erg.push_back(upperThreshold_);
    return erg;
}

std::vector<std::string> PlotPredicateNotBetweenOrEqual::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Lower Threshold");
    erg.push_back("Upper Threshold");
    return erg;
}

bool PlotPredicateNotBetweenOrEqual::check(const PlotCellValue& value) const {
    if (value.isValue()) {
        return (   lowerThreshold_.isValue()
                && upperThreshold_.isValue()
                && (value.getValue() < lowerThreshold_.getValue()
                || value.getValue() > upperThreshold_.getValue()));
    }
    else if (value.isTag()) {
        return (   lowerThreshold_.isTag()
                && upperThreshold_.isTag()
                && (value.getTag() < lowerThreshold_.getTag()
                || value.getTag() > upperThreshold_.getTag()));
    }
    return false;
}

Interval<plot_t> PlotPredicateNotBetweenOrEqual::getIntervalRepresentation() const {
    return Interval<plot_t>(upperThreshold_.getValue(),lowerThreshold_.getValue(), true, true);
}

PlotPredicate* PlotPredicateNotBetweenOrEqual::clone() const {
    return new PlotPredicateNotBetweenOrEqual(lowerThreshold_,upperThreshold_);
}

std::string PlotPredicateNotBetweenOrEqual::toString() const {
    std::stringstream ss;
    ss << "not in [";
    if (lowerThreshold_.isNull())
        ss << "NULL";
    else if (lowerThreshold_.isValue()) {
        setSmartPrecision(ss, 3, lowerThreshold_.getValue());
        ss << lowerThreshold_.getValue();
    }
    else
        ss << "\"" << lowerThreshold_.getTag() << "\"";
    ss << ", ";
    if (upperThreshold_.isNull())
        ss << "NULL";
    else if (upperThreshold_.isValue()) {
        setSmartPrecision(ss, 3, upperThreshold_.getValue());
        ss << upperThreshold_.getValue();
    }
    else
        ss << "\"" << upperThreshold_.getTag() << "\"";
    ss << "]";
    return ss.str();
}

void PlotPredicateNotBetweenOrEqual::serialize(Serializer& s) const {
    s.serialize("lowerThreshold",lowerThreshold_);
    s.serialize("upperThreshold",upperThreshold_);
}

void PlotPredicateNotBetweenOrEqual::deserialize(Deserializer& s) {
    s.deserialize("lowerThreshold",lowerThreshold_);
    s.deserialize("upperThreshold",upperThreshold_);
}

// PlotPredicateIsSubStr methods -------------------------------------------------------

PlotPredicateIsSubStr::PlotPredicateIsSubStr()
    : threshold_("")
    {}

PlotPredicateIsSubStr::PlotPredicateIsSubStr(std::string threshold)
    : threshold_(threshold)
    {}

bool PlotPredicateIsSubStr::check(const PlotCellValue& value) const {
    return (value.isTag() && (value.getTag().find(threshold_,0) != std::string::npos));
}


int PlotPredicateIsSubStr::getNumberOfThresholdValues() const {
    return 1;
}

void PlotPredicateIsSubStr::setThresholdValues(const std::vector<PlotCellValue>& thresholds) {
    if (! thresholds.empty())
        threshold_ = thresholds[0].getTag();
}

std::vector<PlotCellValue> PlotPredicateIsSubStr::getThresholdValues() const {
    std::vector<PlotCellValue> erg;
    erg.push_back(PlotCellValue(threshold_));
    return erg;
}

std::vector<std::string> PlotPredicateIsSubStr::getThresholdTitles() const {
    std::vector<std::string> erg;
    erg.push_back("Substring");
    return erg;
}

Interval<plot_t> PlotPredicateIsSubStr::getIntervalRepresentation() const {
    return Interval<plot_t>();
}

PlotPredicate* PlotPredicateIsSubStr::clone() const {
    return new PlotPredicateIsSubStr(threshold_);
}

std::string PlotPredicateIsSubStr::toString() const {
    std::stringstream ss;
    ss << "Substring: " << "\"" << threshold_ << "\"";
    return ss.str();
}

void PlotPredicateIsSubStr::serialize(Serializer& s) const {
    s.serialize("threshold",threshold_);
}

void PlotPredicateIsSubStr::deserialize(Deserializer& s) {
    s.deserialize("threshold",threshold_);
}

// PlotPredicateNotAlphaNumeric methods -------------------------------------------------------

bool PlotPredicateNotAlphaNumeric::check(const PlotCellValue& value) const {
    if (value.isValue())
        return false;
    if (value.isNull())
        return false;
    if (value.isTag()) {
        std::string s = value.getTag();
        std::stringstream Str(s);
        double num;
        if (Str >> num)
            return false;
        for (size_t i = 0; i < s.length(); ++i) {
            int c = s[i];
            if (c <= 0 || !isalnum(s[i])) {
                return true;
            }
        }
        return false;
    }

    return true;
}

Interval<plot_t> PlotPredicateNotAlphaNumeric::getIntervalRepresentation() const {
    return Interval<plot_t>();
}

PlotPredicate* PlotPredicateNotAlphaNumeric::clone() const {
    return new PlotPredicateNotAlphaNumeric();
}

int PlotPredicateNotAlphaNumeric::getNumberOfThresholdValues() const {
    return 0;
}

void PlotPredicateNotAlphaNumeric::setThresholdValues(const std::vector<PlotCellValue>& /*thresholds*/) {
}

std::vector<PlotCellValue> PlotPredicateNotAlphaNumeric::getThresholdValues() const {
    return std::vector<PlotCellValue>();
}

std::vector<std::string> PlotPredicateNotAlphaNumeric::getThresholdTitles() const {
    return std::vector<std::string>();
}

std::string PlotPredicateNotAlphaNumeric::toString() const {
    return "Not Alpha-Numeric";
}

void PlotPredicateNotAlphaNumeric::serialize(Serializer& /*s*/) const {
}

void PlotPredicateNotAlphaNumeric::deserialize(Deserializer& /*s*/) {
}

// PlotPredicateAlphaNumeric methods -------------------------------------------------------

bool PlotPredicateAlphaNumeric::check(const PlotCellValue& value) const {
    if (value.isValue())
        return true;
    if (value.isNull())
        return false;
    if (value.isTag()) {
        std::string s = value.getTag();
        std::stringstream Str(s);
        double num;
        if (Str >> num)
            return true;
        for (size_t i = 0; i < s.length(); ++i) {
            int c = s[i];
            if (c <= 0 || !isalnum(s[i])) {
                return false;
            }
        }
        return true;
    }

    return false;
}

Interval<plot_t> PlotPredicateAlphaNumeric::getIntervalRepresentation() const {
    return Interval<plot_t>();
}

PlotPredicate* PlotPredicateAlphaNumeric::clone() const {
    return new PlotPredicateAlphaNumeric();
}

int PlotPredicateAlphaNumeric::getNumberOfThresholdValues() const {
    return 0;
}

void PlotPredicateAlphaNumeric::setThresholdValues(const std::vector<PlotCellValue>& /*thresholds*/) {
}

std::vector<PlotCellValue> PlotPredicateAlphaNumeric::getThresholdValues() const {
    return std::vector<PlotCellValue>();
}

std::vector<std::string> PlotPredicateAlphaNumeric::getThresholdTitles() const {
    return std::vector<std::string>();
}

std::string PlotPredicateAlphaNumeric::toString() const {
    return "Alpha Numeric";
}

void PlotPredicateAlphaNumeric::serialize(Serializer& /*s*/) const {
}

void PlotPredicateAlphaNumeric::deserialize(Deserializer& /*s*/) {
}

// PlotPredicateEmpty methods -------------------------------------------------------

bool PlotPredicateEmpty::check(const PlotCellValue& value) const {
    return (value.isNull());
}

Interval<plot_t> PlotPredicateEmpty::getIntervalRepresentation() const {
    return Interval<plot_t>();
}

PlotPredicate* PlotPredicateEmpty::clone() const {
    return new PlotPredicateEmpty();
}

int PlotPredicateEmpty::getNumberOfThresholdValues() const {
    return 0;
}

void PlotPredicateEmpty::setThresholdValues(const std::vector<PlotCellValue>& /*thresholds*/) {
}

std::vector<PlotCellValue> PlotPredicateEmpty::getThresholdValues() const {
    return std::vector<PlotCellValue>();
}

std::vector<std::string> PlotPredicateEmpty::getThresholdTitles() const {
    return std::vector<std::string>();
}

std::string PlotPredicateEmpty::toString() const {
    return "Is Empty";
}

void PlotPredicateEmpty::serialize(Serializer& /*s*/) const {
}

void PlotPredicateEmpty::deserialize(Deserializer& /*s*/) {
}

// PlotPredicateNotEmpty methods -------------------------------------------------------

bool PlotPredicateNotEmpty::check(const PlotCellValue& value) const {
    return (!value.isNull());
}

Interval<plot_t> PlotPredicateNotEmpty::getIntervalRepresentation() const {
    return Interval<plot_t>(-std::numeric_limits<plot_t>::max(), std::numeric_limits<plot_t>::max(), false, false);
}

PlotPredicate* PlotPredicateNotEmpty::clone() const {
    return new PlotPredicateNotEmpty();
}

int PlotPredicateNotEmpty::getNumberOfThresholdValues() const {
    return 0;
}

void PlotPredicateNotEmpty::setThresholdValues(const std::vector<PlotCellValue>& /*thresholds*/) {
}

std::vector<PlotCellValue> PlotPredicateNotEmpty::getThresholdValues() const {
    return std::vector<PlotCellValue>();
}

std::vector<std::string> PlotPredicateNotEmpty::getThresholdTitles() const {
    return std::vector<std::string>();
}

std::string PlotPredicateNotEmpty::toString() const {
    return "Is Not Empty";
}

void PlotPredicateNotEmpty::serialize(Serializer& /*s*/) const {
}

void PlotPredicateNotEmpty::deserialize(Deserializer& /*s*/) {
}

}

