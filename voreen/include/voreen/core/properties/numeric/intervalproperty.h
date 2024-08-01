/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_INTINTERVALPROPERTY_H
#define VRN_INTINTERVALPROPERTY_H

#include "voreen/core/properties/numericproperty.h"
#include "tgt/vector.h"
#include "tgt/interval.h"

namespace voreen {

template<typename TYPE>
class IntervalProperty :public TemplateProperty<tgt::Vector2<TYPE> >, public NumericPropertyBase {
public:
    IntervalProperty(const std::string& id, const std::string& guiText,
        TYPE value = 0, TYPE minValue = 0, TYPE maxValue = 100,
        TYPE minRange=0, TYPE maxRange=INT_MAX, int invalidationLevel=Processor::INVALID_RESULT,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    IntervalProperty(const std::string& id, const std::string& guiText,
        typename tgt::Vector2<TYPE>  value, TYPE minValue = 0, TYPE maxValue = 100,
        TYPE minRange=0, TYPE maxRange=INT_MAX, int invalidationLevel=Processor::INVALID_RESULT,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    IntervalProperty();

    void set(const TYPE &intervalmin);

    void set(const typename tgt::Vector2<TYPE>& v);

    const tgt::Vector2<TYPE>& get() const;

    /**
     * The minimal value for the interval
     */
    TYPE getMinValue() const;
    /**
     * The maximal value for the interval
     */
    TYPE getMaxValue() const;
    /**
     * The minimal size of th interval range
     */
    TYPE getMinRange() const;
    /**
     * The maximal size of the interval range
     */
    TYPE getMaxRange() const;


    /**
     * Set the minimium value of the interval
     */
    void setMinValue(TYPE value);

    /**
     * Set the maximium value of the interval
     */
    void setMaxValue(TYPE value);

    /**
     * Set the minimium range of the interval
     */
    void setMinRange(TYPE value);

    /**
     * Set the maximium range of the interval
     */
    void setMaxRange(TYPE value);

    /**
     * Sets the property's minimum increase/decrease,
     * which should be used in GUI widgets representing
     * the property, e.g., a spin box.
     */
    void setStepping(TYPE value);
    TYPE getStepping() const;

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
     * Get the internal tgt::Interval
     */
    tgt::Interval<TYPE> getInterval() const;

    /**
     * Set the internal tgt::Interval
     */
    void setInterval(tgt::Interval<TYPE> interval);

    virtual void reset();

    /// @see Property::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Property::deserialize
    virtual void deserialize(Deserializer& s);

private:
    tgt::Interval<TYPE> interval_;
    tgt::Interval<TYPE> defaultValue_;
    TYPE stepping_;
    int numDecimals_;
    bool tracking_;
};

class IntIntervalProperty :public IntervalProperty<int>{
public:
    IntIntervalProperty(const std::string& id, const std::string& guiText,
        int value = 0, int minValue = 0, int maxValue = 100,
        int minRange=0, int maxRange=INT_MAX, int invalidationLevel=Processor::INVALID_RESULT,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT)
        :IntervalProperty<int>(id, guiText, value, minValue, maxValue,
                          minRange, maxRange, invalidationLevel, lod){
    }
    IntIntervalProperty(const std::string& id, const std::string& guiText,
        tgt::Vector2<int>  value, int minValue = 0, int maxValue = 100,
        int minRange=0, int maxRange=INT_MAX, int invalidationLevel=Processor::INVALID_RESULT,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT)
        :IntervalProperty<int>(id, guiText, value, minValue, maxValue,
                          minRange, maxRange, invalidationLevel, lod){
    }
    IntIntervalProperty(){}

    virtual std::string getClassName() const{
        return "IntIntervalProperty";
    }
    virtual std::string getTypeDescription() const{
        return "Integer interval";
    }
    virtual Property* create() const{
        return new IntIntervalProperty();
    }
};

class FloatIntervalProperty :public IntervalProperty<float>{
public:
    FloatIntervalProperty(const std::string& id, const std::string& guiText,
        float value = 0, float minValue = 0, float maxValue = 100,
        float minRange=0, float maxRange=FLT_MAX, int invalidationLevel=Processor::INVALID_RESULT,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT)
        :IntervalProperty<float>(id, guiText, value, minValue, maxValue,
                          minRange, maxRange, invalidationLevel, lod){
        setStepping(0.01f);
    }
    FloatIntervalProperty(const std::string& id, const std::string& guiText,
        tgt::Vector2<float>  value, float minValue = 0, float maxValue = 100,
        float minRange=0, float maxRange=FLT_MAX, int invalidationLevel=Processor::INVALID_RESULT,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT)
        :IntervalProperty<float>(id, guiText, value, minValue, maxValue,
                          minRange, maxRange, invalidationLevel, lod){
        setStepping(0.01f);
    }
    FloatIntervalProperty(){}

    virtual std::string getClassName() const{
        return "FloatIntervalProperty";
    }
    virtual std::string getTypeDescription() const{
        return "Float interval";
    }
    virtual Property* create() const{
        return new FloatIntervalProperty();
    }

    /*
     * Adapt stepping and numDecimals to the current range (max - min),
     * so that numSignificantDecimals decimals will be adjustable.
     */
    void adaptDecimalsToRange(size_t numSignificantDecimals) {
        const float range = getMaxValue() - getMinValue();
        if(range <= 0) {
            return;
        }
        // We can only omit digits behind the period, digits has to be >= 0.
        int digits = std::max(0, tgt::iround(-log10(range)) + static_cast<int>(numSignificantDecimals));
        setNumDecimals(digits);
        setStepping(pow(10.0f, -digits));

        // Hack: If we set min or max to a very specific value before,
        // the GUI representation might be slightly off, because it is
        // affected by the stepping, but the actual min/max values are not.
        // So now that we have set the stepping appropriately, we force
        // an update by changing min/max around and back to their original
        // values.
        float min = getMinValue();
        float max = getMaxValue();
        setMinValue(std::numeric_limits<float>::lowest());
        setMaxValue(std::numeric_limits<float>::max());
        setMinValue(min);
        setMaxValue(max);
    }

};

template<typename TYPE>
class IntervalPropertyValidation : public Condition{
public:
    IntervalPropertyValidation(IntervalProperty<TYPE>* prop);
    virtual Condition* clone() const;
    virtual bool met() const throw();
    virtual std::string description() const;
private:
    IntervalProperty<TYPE> *prop_;
};



}   // namespace

 /**********************************************************
  *                                                        *
  *                    Implementation                      *
  *                                                        *
  **********************************************************/
#include "tgt/assert.h"
#include <sstream>
namespace voreen {

template<typename TYPE>
IntervalProperty<TYPE>::IntervalProperty(const std::string& id, const std::string& guiText,
                    TYPE value, TYPE minValue, TYPE maxValue, TYPE minRange, TYPE maxRange,
                    int invalidationLevel, Property::LevelOfDetail lod)
    : TemplateProperty<tgt::Vector2<TYPE> >(id, guiText, tgt::Vector2<TYPE>(value), invalidationLevel, lod)
    , interval_(tgt::Vector2<TYPE>(value), minValue, maxValue, minRange, maxRange)
    , defaultValue_(interval_)
    , stepping_(1)
    , numDecimals_(2)
{
    //this->addValidation(IntervalPropertyValidation<TYPE>(this));
}

template<typename TYPE>
IntervalProperty<TYPE>::IntervalProperty(const std::string& id, const std::string& guiText,
                    typename tgt::Vector2<TYPE>  value, TYPE minValue, TYPE maxValue, TYPE minRange, TYPE maxRange,
                    int invalidationLevel, Property::LevelOfDetail lod)
    : TemplateProperty<tgt::Vector2<TYPE> >(id, guiText, value, invalidationLevel, lod)
    , interval_(value, minValue, maxValue, minRange, maxRange)
    , defaultValue_(interval_)
    , stepping_(1)
    , numDecimals_(2)
{
//    this->addValidation(IntervalPropertyValidation<TYPE>(this));
}

template<typename TYPE>
IntervalProperty<TYPE>::IntervalProperty()
    : TemplateProperty<tgt::Vector2<TYPE> > ()
    , defaultValue_(interval_)
    , stepping_(1)
    , numDecimals_(2)
{
    //this->addValidation(IntervalPropertyValidation<TYPE>(this));
}

template<typename TYPE>
void IntervalProperty<TYPE>::set(const TYPE &intervalmin){
    if (tgt::Vector2<TYPE>(intervalmin, intervalmin+interval_.getMinRange()) != interval_.get()) {
        interval_.set(tgt::Vector2<TYPE> (intervalmin, intervalmin+interval_.getMinRange()));
        TemplateProperty<tgt::Vector2<TYPE> >::set(interval_.get());
    }
}

template<typename TYPE>
void IntervalProperty<TYPE>::set(const tgt::Vector2<TYPE>& val){
    if (val == interval_.get())
        return;
    if (interval_.isValid(val)){
        interval_.set(val);
        TemplateProperty<tgt::Vector2<TYPE> >::set(interval_.get());
    }else{
        LWARNINGC("voreen.IntervalProperty", "Set invalid value.");
    }
}

template<typename TYPE>
const tgt::Vector2<TYPE>& IntervalProperty<TYPE>::get() const{
    const_cast<IntervalProperty*>(this)->TemplateProperty<tgt::Vector2<TYPE> >::set(interval_.get()); // SHOULD do nothing
    return TemplateProperty<tgt::Vector2<TYPE> >::get();
}


template<typename TYPE>
TYPE IntervalProperty<TYPE>::getMinValue() const{
    return interval_.getMinValue();
}

template<typename TYPE>
TYPE IntervalProperty<TYPE>::getMaxValue() const{
    return interval_.getMaxValue();
}

template<typename TYPE>
TYPE IntervalProperty<TYPE>::getMinRange() const{
    return interval_.getMinRange();
}

template<typename TYPE>
TYPE IntervalProperty<TYPE>::getMaxRange() const{
    return interval_.getMaxRange();
}

template<typename TYPE>
void IntervalProperty<TYPE>::setMinValue(TYPE value){
    if (interval_.getMinValue() != value) {
        interval_.setMinValue(value);
        this->invalidate();
    }
}

template<typename TYPE>
void IntervalProperty<TYPE>::setMaxValue(TYPE value){
    if (interval_.getMaxValue() != value) {
        interval_.setMaxValue(value);
        this->invalidate();
    }
}

template<typename TYPE>
void IntervalProperty<TYPE>::setMinRange(TYPE value){
    if (interval_.getMinRange() != value) {
        interval_.setMinRange(value);
        this->invalidate();
    }
}

template<typename TYPE>
void IntervalProperty<TYPE>::setMaxRange(TYPE value){
    if (interval_.getMaxRange() != value) {
        interval_.setMaxRange(value);
        this->invalidate();
    }
}

template<typename TYPE>
tgt::Interval<TYPE> IntervalProperty<TYPE>::getInterval() const{
    return interval_;
}

template<typename TYPE>
void IntervalProperty<TYPE>::setInterval(tgt::Interval<TYPE> interval){
    if (interval != interval_) {
        interval_ = interval;
        this->invalidate();
    }
}

template<typename TYPE>
void IntervalProperty<TYPE>::reset(){
    interval_ = defaultValue_;
    this->invalidate();
}

template<typename TYPE>
void IntervalProperty<TYPE>::serialize(Serializer& s) const{
    tgt::Vector2<TYPE> t = interval_.get();
    s.serialize("value", t);
    s.serialize("minValue", this->getMinValue());
    s.serialize("maxValue", this->getMaxValue());
    s.serialize("minRange", this->getMinRange());
    s.serialize("maxRange", this->getMaxRange());
    s.serialize("stepping", this->getStepping());
    s.serialize("tracking", this->hasTracking());
}

template<typename TYPE>
void IntervalProperty<TYPE>::deserialize(Deserializer& s){
    typename tgt::Vector2<TYPE> t;
    s.optionalDeserialize("value", t, tgt::Vector2<TYPE>::zero);

    TYPE minValue, maxValue, minRange, maxRange, stepping;
    bool tracking;

    s.optionalDeserialize("minValue", minValue, static_cast<TYPE>(0));
    s.optionalDeserialize("maxValue", maxValue, static_cast<TYPE>(1));
    s.optionalDeserialize("minRange", minRange, static_cast<TYPE>(0));
    s.optionalDeserialize("maxRange", maxRange, std::numeric_limits<TYPE>::max());
    s.optionalDeserialize("stepping", stepping, static_cast<TYPE>(1));
    s.optionalDeserialize("tracking", tracking, tracking_);

    this->setMinValue(minValue);
    this->setMaxValue(maxValue);
    this->setMinRange(minRange);
    this->setMaxRange(maxRange);
    this->setStepping(stepping);
    this->setTracking(tracking);

    interval_.set(t);
}

template<typename TYPE>
void IntervalProperty<TYPE>::setStepping(TYPE stepping){
    if (stepping != stepping_){
        stepping_ = stepping;
        this->invalidate();
    }

}

template<typename TYPE>
TYPE IntervalProperty<TYPE>::getStepping() const{
    return stepping_;
}

template<typename TYPE>
void IntervalProperty<TYPE>::setNumDecimals(int numDecimals) {
    tgtAssert(numDecimals <= 64 && numDecimals >= 0, "Invalid number of decimals");
    numDecimals_ = numDecimals;
    this->updateWidgets();
}

template<typename TYPE>
int IntervalProperty<TYPE>::getNumDecimals() const {
    return numDecimals_;
}

template<typename TYPE>
bool IntervalProperty<TYPE>::hasTracking() const {
    return tracking_;
}

template<typename TYPE>
void IntervalProperty<TYPE>::setTracking(bool tracking) {
    tracking_ = tracking;
    this->updateWidgets();
}


template<typename TYPE>
IntervalPropertyValidation<TYPE>::IntervalPropertyValidation(IntervalProperty<TYPE>* prop)
    :prop_(prop)
{
}

template<typename TYPE>
Condition* IntervalPropertyValidation<TYPE>::clone() const{
    return new IntervalPropertyValidation(*this);
}

template<typename TYPE>
bool IntervalPropertyValidation<TYPE>::met() const throw(){
    TYPE min      = prop_->get().x;
    TYPE max      = prop_->get().y;
    TYPE minval   = prop_->getMinValue();
    TYPE maxval   = prop_->getMaxValue();
    TYPE minrange = prop_->getMinRange();
    TYPE maxrange = prop_->getMaxRange();

    if (minrange < 0)
        return false;
    if (maxrange < 0)
        return false;
    if (minval > maxval)
        return false;
    if (min < minval)
        return false;
    if (min > maxval)
        return false;
    if (max < minval)
        return false;
    if (max > maxval)
        return false;
    if (max-min < minrange)
        return false;
    if (max-min > maxrange)
        return false;
    return true;
}

template<typename TYPE>
std::string IntervalPropertyValidation<TYPE>::description() const{
    std::stringstream s;
    TYPE min      = prop_->get().x;
    TYPE max      = prop_->get().y;
    TYPE minval   = prop_->getMinValue();
    TYPE maxval   = prop_->getMaxValue();
    TYPE minrange = prop_->getMinRange();
    TYPE maxrange = prop_->getMaxRange();


    if (minrange < 0)
        s<< "minRange should be > 0; is ["<<minrange<<"].";
    if (maxrange < 0)
        s<< "maxRange should be > 0; is ["<<maxrange<<"].";
    if (minval > maxval)
        s << "minValue is greater then maxValue.";
    if (min < minval)
        s << "min["<<min<<"] should not be less then minValue["<<minval<<"].";
    if (min > maxval)
        s << "min["<<min<<"] should not be greater then maxValue["<<maxval<<"].";
    if (max < minval)
        s << "max["<<max<<"] should not be less then minValue["<<minval<<"].";
    if (max > maxval)
        s << "max["<<max<<"] should not be greater then maxValue["<<maxval<<"].";
    if (max-min < minrange)
        s << "interval["<< min <<","<< max << "] should not be smaller then minRange["<<minrange<<"].";
    if (max-min > maxrange)
        s << "interval["<< min <<","<< max << "] should not be bigger then maxRange["<<maxrange<<"].";
    return s.str();
}
}   // namespace

#endif
