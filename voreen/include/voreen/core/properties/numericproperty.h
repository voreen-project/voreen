/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_NUMERICPROPERTY_H
#define VRN_NUMERICPROPERTY_H

#include "voreen/core/properties/condition.h"
#include "voreen/core/properties/templateproperty.h"

namespace voreen {
/**
 * Template for properties that store numeric values.
 *
 * Besides the allowed value range, some GUI hints can be defined,
 * such as the stepping and number of decimals.
 */
template<typename T>
class NumericProperty : public TemplateProperty<T> {

template<typename> friend class NumericPropertyValidation;

public:
    typedef T ElemType;

    /**
     * Policy, how property boundaries [min and max value, range] should be handled.
     * Defualt value is STATIC.
     */
    enum BoundaryUpdatePolicy {
        CONSTANT,   ///< The boundaries can not be changed.
        STATIC,     ///< Default. Boundaries can be changed at runtime, but will not be serialized.
        DYNAMIC     ///< Boundaries can be changed at runtime and but will be serialized/deserialized.
    };

    NumericProperty(const std::string& id, const std::string& guiText, const T& value,
                    const T& minValue, const T& maxValue, const T& stepping,
                    int invalidationLevel=Processor::INVALID_RESULT, BoundaryUpdatePolicy boundPol = STATIC,
                    Property::LevelOfDetail lod = Property::LOD_DEFAULT);

    /**
     * Sets the minimum value the variable can take.
     * This is checked by a NumericPropertyValidation on each
     * set() call. Therefore, it is guaranteed that the property
     * never attains a value outside its valid range.
     */
    void setMinValue(const T& minValue);
    const T& getMinValue() const;

    /**
     * Sets the maximum value the variable can take.
     * This is checked by a NumericPropertyValidation on each
     * set() call. Therefore, it is guaranteed that the property
     * never attains a value outside its valid range.
     */
    void setMaxValue(const T& maxValue);
    const T& getMaxValue() const;

    /**
     * Sets the property's minimum increase/decrease,
     * which should be used in GUI widgets representing
     * the property, e.g., a spin box.
     *
     * @see increase, decrease
     */
    void setStepping(const T stepping);
    T getStepping() const;

    /**
     * Sets the number of decimals that should be
     * displayed by a GUI representation of the property.
     */
    void setNumDecimals(int numDecimals);
    int getNumDecimals() const;

    /**
     * If tracking is disabled, the property is not to be updated
     * during user interactions, e.g., while the user drags a slider.
     * Tracking is enabled by default.
     */
    void setTracking(bool tracking);
    bool hasTracking() const;

    /**
     * Increases the property's value by its stepping.
     * @see setStepping
     */
    void increase();

    /**
     * Decreases the property's value by its stepping.
     * @see setStepping
     */
    void decrease();

    /// @see Property::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Property::deserialize
    virtual void deserialize(Deserializer& s);

protected:
    using TemplateProperty<T>::value_;

    T minValue_;
    T maxValue_;
    T stepping_;
    bool tracking_;
    int numDecimals_;
    BoundaryUpdatePolicy boundaryPolicy_;
};

} // namespace

//HACK: resolving circular dependencie
#include "voreen/core/properties/link/linkevaluatoridgeneric.h"

namespace voreen {

template<typename T>
NumericProperty<T>::NumericProperty(const std::string& id, const std::string& guiText, const T& value,
                                    const T& minValue, const T& maxValue, const T& stepping,
                                    int invalidationLevel, BoundaryUpdatePolicy boundPol, Property::LevelOfDetail lod)
                                    : TemplateProperty<T>(id, guiText, value, invalidationLevel, lod),
                                    minValue_(minValue),
                                    maxValue_(maxValue),
                                    stepping_(stepping),
                                    tracking_(true),
                                    numDecimals_(2),
                                    boundaryPolicy_(boundPol)
{
    this->addValidation(NumericPropertyValidation<T>(this));
    std::string errorMsg;
    this->validate(value, errorMsg);
    if (!errorMsg.empty())
        LWARNINGC("voreen.NumericProperty", errorMsg);
}

template<typename T>
void NumericProperty<T>::setMaxValue( const T& maxValue ) {
    if(boundaryPolicy_ == CONSTANT) {
        LWARNINGC("voreen.NumericProperty", "Tried setting max value on numeric property with constant boundaries.");
        return;
    }

    if(maxValue_ != maxValue) {
        maxValue_ = maxValue;
        //check if min < max (so complex because of  T == vec3 etc.)
        T tempValue(tgt::clamp(minValue_,minValue_,maxValue_));
        if (tempValue != minValue_){
            LDEBUGC("voreen.NumericProperty",
                    "setMaxValue: maxValue < minValue. Setting minValue = min(minValue,maxValue).");
            this->setMinValue(tempValue);
            return; //invalidate has been called in setMinValue
        }
        //check, if value_ < maxValue_ (so complex because of  T == vec3 etc.)
        tempValue = tgt::clamp(value_,value_,maxValue_);
        if(value_ != tempValue)
            this->set(tempValue); //calls invalidate
        else
            this->invalidate();
    }
}

template<typename T>
const T& NumericProperty<T>::getMaxValue() const {
    return maxValue_;
}

template<typename T>
void NumericProperty<T>::setMinValue(const T& minValue) {
    if(boundaryPolicy_ == CONSTANT) {
        LWARNINGC("voreen.NumericProperty", "Tried setting min value on numeric property with constant boundaries.");
        return;
    }

    if(minValue_ != minValue) {
        minValue_ = minValue;
        //check if min < max (so complex because of  T == vec3 etc.)
        T tempValue(-tgt::clamp(-minValue_,-minValue_,-maxValue_));
        if (tempValue != maxValue_){
            LDEBUGC("voreen.NumericProperty",
                    "setMinValue: maxValue < minValue. Setting maxValue = max(minValue,maxValue).");
            this->setMaxValue(tempValue);
            return; //invalidate has been called in setMaxValue
        }
        //check, if value_ > minValue_ (so complex because of  T == vec3 etc.)
        tempValue = -tgt::clamp(-minValue_,-minValue_,-value_);
        if(value_ != tempValue)
            this->set(tempValue); //calls invalidate
        else
            this->invalidate();
    }
}

template<typename T>
const T& NumericProperty<T>::getMinValue() const {
    return minValue_;
}

template<typename T>
void NumericProperty<T>::setStepping(const T stepping) {
    stepping_ = stepping;
    this->updateWidgets();
}

template<typename T>
T NumericProperty<T>::getStepping() const {
    return stepping_;
}

template<typename T>
int NumericProperty<T>::getNumDecimals() const {
    return numDecimals_;
}

template<typename T>
void NumericProperty<T>::setNumDecimals(int numDecimals) {
    tgtAssert(numDecimals <= 64 && numDecimals >= 0, "Invalid number of decimals");
    numDecimals_ = numDecimals;
    this->updateWidgets();
}

template<typename T>
bool NumericProperty<T>::hasTracking() const {
    return tracking_;
}

template<typename T>
void NumericProperty<T>::setTracking(bool tracking) {
    tracking_ = tracking;
    this->updateWidgets();
}

template<typename T>
void NumericProperty<T>::increase() {
    this->set(value_ + stepping_);
}

template<typename T>
void NumericProperty<T>::decrease() {
    this->set(value_ - stepping_);
}

template<typename T>
void NumericProperty<T>::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("value", value_);

    if(boundaryPolicy_ == DYNAMIC) {
        s.serialize("minValue", minValue_);
        s.serialize("maxValue", maxValue_);
    }

    // serialize tracking mode, if it differs from default value
    if (!tracking_)
        s.serialize("tracking", tracking_);
}

template<typename T>
void NumericProperty<T>::deserialize(Deserializer& s) {
    Property::deserialize(s);

    // deserialize value
    T value;
    s.deserialize("value", value);

    try {
        if(boundaryPolicy_ == DYNAMIC) {
            s.deserialize("minValue", minValue_);
            s.deserialize("maxValue", maxValue_);
        }
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
    }

    try {
        this->set(value);
    }
    catch (Condition::ValidationFailed& e) {
        s.addError(e);
    }

    // deserialize tracking mode, if available
    try {
        s.deserialize("tracking", tracking_);
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
    }
}

} // namespace voreen

#endif // VRN_NUMERICPROPERTY_H
