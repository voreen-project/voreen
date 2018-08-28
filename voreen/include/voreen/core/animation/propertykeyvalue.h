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

#ifndef VRN_PROPERTYKEYVALUE_H
#define VRN_PROPERTYKEYVALUE_H

#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/animation/interpolationfunction.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"
#include "tgt/camera.h"

namespace voreen {

template <class T> class MultiPointInterpolationFunction;
template <class T> class TemplatePropertyTimelineState;

/**
 * Non-templated base class for generic PropertyKeyValues.
 */
class VRN_CORE_API PropertyKeyValueBase : public Serializable {

public:
    PropertyKeyValueBase() {}
    virtual float getTime() const = 0;
};

/**
 * Class for creating an animation.
 * The KeyValue stores one property value at one corresponding time.
 */
template <class T>
class PropertyKeyValue : public PropertyKeyValueBase {
public:
    /**
     * The key value has to be initialized with its value.
     */
    PropertyKeyValue(T value, float time);

    /**
     * Destructor.
     */
    virtual ~PropertyKeyValue();

    /**
     * Function to get the time.
     */
    virtual float getTime() const;

    /**
     * Function to get the value.
     */
    const T getValue() const;

    /**
     * Gets the interpolation function after the keyvalue.
     */
    const InterpolationFunction<T>* getFollowingInterpolationFunction() const;

    /**
     * Gets the interpolation function before the keyvalue.
     */
    const InterpolationFunction<T>* getForegoingInterpolationFunction() const;

    /**
     * Returns the variable smooth.
     * Smoothness is a flag of a key value enabling continuity of the respective interpolation function at that specific key value
     */
    bool getSmooth() const;

    /**
     * Returns true if smooth is true and the following and the foregoing interpolation function are of the same type of multipointinterpolation.
     */
    bool isSmoothed() const;

    /**
     * Returns a copy of the current keyvalue.
     * The value and the time will be copied, but not the linked interpolation functions.
     */
    PropertyKeyValue<T>* clone() const;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

protected:
    friend class TemplatePropertyTimelineState<T>;
    friend class TransFunc1DKeysPropertyTimelineState;
    friend class TransFunc2DPrimitivesPropertyTimelineState;
    friend class XmlDeserializer;

    /**
     * Default constructor.
     */
    PropertyKeyValue();

    /**
     * Function to set the value.
     */
    void setValue(T value);

    /**
     * Sets the inner variable smooth.
     * Explanation of this variable: see function isSmoothed().
     */
    void setSmooth(bool smooth);

    /**
     * Function to set the time.
     */
    void setTime(float time);

    /**
     * Sets the interpolation function before the key value.
     */
    void setForegoingInterpolationFunction(InterpolationFunction<T>* func);

    /**
     * Sets the interpolation function following the key value.
     */
    void setFollowingInterpolationFunction(InterpolationFunction<T>* func);

    float time_;
    T value_;
    InterpolationFunction<T>* before_;
    InterpolationFunction<T>* after_;
    bool smooth_;
};


template <class T>
PropertyKeyValue<T>::PropertyKeyValue()
    : time_(0)
    , before_(0)
    , after_(0)
    , smooth_(true)
{}

template <class T>
PropertyKeyValue<T>::PropertyKeyValue(T value, float time)
    : value_(value)
    , before_(0)
    , after_(0)
    , smooth_(true)
{
    //time_ = floor(time * 10000.f) / 10000.f;
    time_ = time;
}

template <class T>
PropertyKeyValue<T>::~PropertyKeyValue() {}

template <class T>
void PropertyKeyValue<T>::setTime(float time) {
    //time_ = floor(time * 10000.f) / 10000.f;
    time_ = time;
}

template <class T>
float PropertyKeyValue<T>::getTime() const {
    return time_;
}

template <class T>
void PropertyKeyValue<T>::setValue(T value) {
    value_ = value;
}

template <class T>
const T PropertyKeyValue<T>::getValue() const {
    return value_;
}

template <class T>
void PropertyKeyValue<T>::setFollowingInterpolationFunction(voreen::InterpolationFunction<T>* func) {
    after_ = func;
}

template <class T>
void PropertyKeyValue<T>::setForegoingInterpolationFunction(voreen::InterpolationFunction<T>* func) {
    before_ = func;
}

template <class T>
const voreen::InterpolationFunction<T>* PropertyKeyValue<T>::getFollowingInterpolationFunction() const {
    return after_;
}

template <class T>
const voreen::InterpolationFunction<T>* PropertyKeyValue<T>::getForegoingInterpolationFunction() const {
    return before_;
}

template <>
PropertyKeyValue<TransFunc1DKeys*>* PropertyKeyValue<TransFunc1DKeys*>::clone() const;

template <>
PropertyKeyValue<TransFunc2DPrimitives*>* PropertyKeyValue<TransFunc2DPrimitives*>::clone() const;

template <class T>
PropertyKeyValue<T>* PropertyKeyValue<T>::clone() const {
    PropertyKeyValue<T>* keyvalue = new PropertyKeyValue<T>();
    keyvalue->after_ = after_;
    keyvalue->before_ = before_;
    keyvalue->value_ = value_;
    keyvalue->time_ = time_;
    keyvalue->smooth_ = smooth_;
    return keyvalue;
}

template <class T>
void PropertyKeyValue<T>::setSmooth(bool smooth) {
    smooth_ = smooth;
}

template <class T>
bool PropertyKeyValue<T>::getSmooth() const {
    return smooth_;
}

template <class T>
bool PropertyKeyValue<T>::isSmoothed() const {
    if (!getSmooth())
        return false;

    MultiPointInterpolationFunction<T>* f1 = dynamic_cast<MultiPointInterpolationFunction<T>*>(before_);
    MultiPointInterpolationFunction<T>* f2 = dynamic_cast<MultiPointInterpolationFunction<T>*>(after_);
    if (f1 && f2 && (typeid(f1) == typeid(f2)))
        return true;
    else
        return false;
}

template <>
void PropertyKeyValue<tgt::Camera>::serialize(Serializer& s) const;

template <class T>
void PropertyKeyValue<T>::serialize(Serializer& s) const {
    s.serialize("time", time_);
    s.serialize("value", value_);
    s.serialize("smooth", smooth_);
    s.serialize("before", before_);
    s.serialize("after", after_);
}

template <>
void PropertyKeyValue<tgt::Camera>::deserialize(Deserializer& s);

template <class T>
void PropertyKeyValue<T>::deserialize(Deserializer& s) {
    s.deserialize("time", time_);
    s.deserialize("value", value_);
    s.deserialize("smooth", smooth_);
    s.deserialize("before", before_);
    s.deserialize("after", after_);
}

} // namespace voreen

#endif
