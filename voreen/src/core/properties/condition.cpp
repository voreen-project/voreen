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

#include "voreen/core/datastructures/callback/callback.h"
#include "voreen/core/properties/condition.h"
#include "voreen/core/properties/numericproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "tgt/logmanager.h"

namespace voreen {

const std::string Condition::loggerCat_("voreen.Condition");

Condition::~Condition() {
}


bool Condition::validate(std::string& errorMsg) const {

    bool valid = met();
    if (!valid) {
        errorMsg = "Property validation failed" + description();
        //LWARNING(msg);
    }
    return met();
}

// ----------------------------------------------------------------------------

Condition::ValidationFailed::ValidationFailed(const Condition* condition)
    : std::exception()
{
    conditionMsg_ = "ValidationFailed";
    if ((condition != 0) && (condition->description().empty() == false))
        conditionMsg_ += ": " + condition->description();
}

const char* Condition::ValidationFailed::what() const throw() {
    return conditionMsg_.c_str();
}

// ============================================================================

template<typename T>
const std::string NumericPropertyValidation<T>::loggerCat_("voreen.NumericPropertyValidation");

template<>
bool NumericPropertyValidation<float>::met() const throw() {
    // if min and max make sense do a real validation
    if (observed_->minValue_ <= observed_->maxValue_) {
        return ((observed_->value_ >= observed_->minValue_)
            && (observed_->value_ <= observed_->maxValue_));
    }
    else
        return true;
}

template<>
bool NumericPropertyValidation<int>::met() const throw() {
    // if min and max make sense do a real validation
    if (observed_->minValue_ <= observed_->maxValue_) {
        return ((observed_->value_ >= observed_->minValue_)
            && (observed_->value_ <= observed_->maxValue_));
    }
    else
        return true;
}

template<>
bool NumericPropertyValidation<double>::met() const throw() {
    // if min and max make sense do a real validation
    if (observed_->minValue_ <= observed_->maxValue_) {
        return ((observed_->value_ >= observed_->minValue_)
            && (observed_->value_ <= observed_->maxValue_));
    }
    else
        return true;
}
template<typename T>
bool NumericPropertyValidation<T>::met() const throw() {
    const T& min = observed_->minValue_;
    const T& max = observed_->maxValue_;

    // do min and max make sense at all?
    // If not just pretend the validation holds
    for (size_t i = 0; i < min.size; ++i) {
        if (min.elem[i] > max.elem[i])
            return true;
    }

    // If we have come this far we can do real validations:
    // test if any component is < min or > max
    const T& val = observed_->value_;
    for (size_t i = 0; i < val.size; ++i) {
        if ((val.elem[i] < min.elem[i]) || (val.elem[i] > max.elem[i]))
            return false;
    }
    return true;
}

template<class T>
std::string voreen::NumericPropertyValidation<T>::description() const{
    std::stringstream stream;
    stream << " for '" << observed_->getFullyQualifiedID() << "': "
           << observed_->value_ << " out of valid range [" << observed_->minValue_ << "," << observed_->maxValue_ << "]";
    return stream.str();
}

// ----------------------------------------------------------------------------

const std::string OptionPropertyValidation::loggerCat_("voreen.OptionPropertyValidation");

bool OptionPropertyValidation::met() const throw() {
    std::vector<std::string> allowedValues = observed_->getKeys();
    return (std::find(allowedValues.begin(), allowedValues.end(), observed_->get()) != allowedValues.end());
}

std::string voreen::OptionPropertyValidation::description() const {
    std::string keyStr;
    std::vector<std::string> allowedValues = observed_->getKeys();
    if (!allowedValues.empty()) {
        keyStr = allowedValues[0];
        for (size_t i=1; i<allowedValues.size(); i++)
            keyStr += ", " + allowedValues[i];
    }

    std::string desc = " for '" + observed_->getFullyQualifiedID() + "': " + "Unknown option key '" + observed_->get() + "'. "
        + "Accepted values: " + keyStr;
    return desc;
}

// ----------------------------------------------------------------------------
// explicit template instantiation to enable distribution of
// implementation of template class methods over .h and .cpp files
//
template class NumericPropertyValidation<float>;
template class NumericPropertyValidation<int>;
template class NumericPropertyValidation<double>;
template class NumericPropertyValidation<tgt::vec2>;
template class NumericPropertyValidation<tgt::vec3>;
template class NumericPropertyValidation<tgt::vec4>;
template class NumericPropertyValidation<tgt::ivec2>;
template class NumericPropertyValidation<tgt::ivec3>;
template class NumericPropertyValidation<tgt::ivec4>;
template class NumericPropertyValidation<tgt::dvec2>;
template class NumericPropertyValidation<tgt::dvec3>;
template class NumericPropertyValidation<tgt::dvec4>;
template class NumericPropertyValidation<tgt::mat2>;
template class NumericPropertyValidation<tgt::mat3>;
template class NumericPropertyValidation<tgt::mat4>;

}   // namespace
