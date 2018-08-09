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

#ifndef VRN_TEMPLATEPROPERTYTIMELINESTATE_H
#define VRN_TEMPLATEPROPERTYTIMELINESTATE_H

#include <vector>
#include <map>
#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/animation/propertykeyvalue.h"
#include "voreen/core/animation/interpolationfunction.h"
#include "voreen/core/animation/interpolationfunctionfactory.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"
#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

class Animation;
class ShaderSource;

//FIXME: move into class!
enum ChangeTimeOfKeyValueReturn {
    KV_NOT_FOUND = 0,
    KV_EQUAL_TO_OLD = 1,
    KV_TIME_CHANGED = 2,
    KV_EXISTS_AT_NEW_TIME = 3,
    KV_MOVED_TO_DIFFERENT_INTERVAL = 4,
    KV_TIME_AFTER_DURATION = 5
};

//FIXME: move into class!
enum DeleteKeyValueReturn {
    KV_NOT_THERE = 0,
    KV_IS_THE_ONLY_ONE = 1,
    KV_DELETED = 2
};

template <class T>
class TemplatePropertyTimeline;

/**
 * This class organizes all animationsettings of one property of the type T.
 */
template <class T>
class VRN_CORE_API TemplatePropertyTimelineState : public Serializable {
public:
    virtual ~TemplatePropertyTimelineState();

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

protected:
    friend class TemplatePropertyTimeline<T>;
    friend class XmlDeserializer;

    /**
     * Constructor.
     * Parameter is the initial value.
     */
    TemplatePropertyTimelineState(PropertyKeyValue<T>* kv);

    /**
     * Returns a map of all keyvalues of the timeline.
     */
    virtual const std::map<float,PropertyKeyValue<T>*>& getKeyValues() const;

    /**
     * Adds a new keyvalue at a certain time.
     * If there already exists a Keyvalue at the given time, this function does nothing.
     *
     * The interplationfunction at the time of the new keyvalue will be deletetd and replaced by two default
     * interpolations before and after the new keyvalue.
     * Returns the pointer to the new keyvalue.
     * If there already existed a keyvalue at the given time a nullpointer is returned.
     */
    virtual const PropertyKeyValue<T>* newKeyValue(float time);

    const InterpolationFunction<T>* getInterpolationFunctionAt(float time);

    /**
     * Changes the value of a keyvalue:
     * Returns true if the value was changed.
     * Returns false if the given keyvalue couldn't be found in the timeline.
     */
    virtual bool changeValueOfKeyValue(T value, const PropertyKeyValue<T>* keyvalue);

    /**
     * Changes the smooth-variable of a keyvalue:
     * Returns true, if the value was changed.
     * Returns false, if the given keyvalue couldn't be found in the timeline.
     */
    virtual bool changeSmoothnessOfKeyValue(bool smooth, const PropertyKeyValue<T>* keyvalue);

    /**
     * Changes the time of a keyvalue:
     * If there is already a keyvalue at the given time, this function deletes the given keyvalue
     * and sets the value of the already existing keyvalue to the deleted value.
     * If the order of the keyvalues changes by setting the new time, the function reacts as if
     * the given keyvalue was deleted and a new one was added at the given time with the same value as the deleted one.
     *
     * returns KV_NOT_FOUND                   if the given keyvalue couldn't be found in the timeline
     * returns KV_EQUAL_TO_OLD                if the new time is equal to the old one
     * returns KV_TIME_CHANGED                if the new time is in the same intervall between the same other keyvalues
     * returns KV_EXISTS_AT_NEW_TIME          if there already exists a new keyvalue at the new time
     * returns KV_MOVED_TO_DIFFERENT_INTERVAL if the new time is in another interval between two other keyvalues.
     */
    virtual ChangeTimeOfKeyValueReturn changeTimeOfKeyValue(float time, const PropertyKeyValue<T>* keyvalue);

    /**
     * Changes the time of a keyvalue and shifts following keyvalues by the same offset.
     *
     * returns KV_NOT_FOUND                   if the given keyvalue couldn't be found in the timeline
     * returns KV_EQUAL_TO_OLD                if the new time is equal to the old one
     */
    virtual ChangeTimeOfKeyValueReturn shiftKeyValue(float time, const PropertyKeyValue<T>* keyvalue);

    /**
     * Deletes the given keyvalue:
     * The two interpolationfunctions before and after the deleted keyvalue will be deleted too and replaced by one default interpolation.
     * returns KV_NOT_FOUND       if the given keyvalue couldn't be found in the timeline
     * returns KV_IS_THE_ONLY_ONE if teh given keyvalue is the only value in the timeline and mustn't be deleted
     * returns KV_DELETED         if the keyvalue was deleted successfully
    */
    virtual DeleteKeyValueReturn deleteKeyValue(const PropertyKeyValue<T>* keyvalue);

    /**
     * Sets the interpolationfunction between the given keyvalue and the successor.
     */
    virtual void setInterpolationFunctionAfter(InterpolationFunction<T>* func, PropertyKeyValue<T>* keyvalue);

    /**
     * Sets the interpolationfunction btween the given keyvalue and its predecessor.
     */
    virtual void setInterpolationFunctionBefore(InterpolationFunction<T>* func, PropertyKeyValue<T>* keyvalue);

    /**
     * Deletes all keyvalues after time of duration.
     */
    virtual void setDuration(float duration);

    /**
     * Returns a clone of the timeline.
     */
    virtual TemplatePropertyTimelineState<T>* clone() const;

    /**
     * Calculates the animated value to a certain point of time.
     */
    virtual const T getPropertyAt(float time);

     /**
      * Default constructor for clone-method.
      */
    TemplatePropertyTimelineState();

    static float round(float x);

    /**
     * Map of all the keyvalues in this timelinestate.
     */
    std::map<float,PropertyKeyValue<T>*> values_;
};

/**
 * Special implementation of the timeline for transfunc-property,
 * because of the fact that the template-variable is a pointer.
 * For a description of the functions see further above in superclass.
 */
class VRN_CORE_API TransFunc1DKeysPropertyTimelineState : public TemplatePropertyTimelineState<TransFunc1DKeys*> {
public:
    TransFunc1DKeysPropertyTimelineState(PropertyKeyValue<TransFunc1DKeys*>* kv);
    virtual ~TransFunc1DKeysPropertyTimelineState();
protected:
    friend class TemplatePropertyTimeline<TransFunc1DKeys*>;
    friend class XmlDeserializer;
    TransFunc1DKeysPropertyTimelineState();
    const PropertyKeyValue<TransFunc1DKeys*>* newKeyValue(float time);
    bool changeValueOfKeyValue(TransFunc1DKeys* value, const PropertyKeyValue<TransFunc1DKeys*>* keyvalue);
    DeleteKeyValueReturn deleteKeyValue(const PropertyKeyValue<TransFunc1DKeys*>* keyvalue);
    /// TransFuncBase* const is const T (not const TransFuncBase*)
    virtual TransFunc1DKeys* const getPropertyAt(float time);
    virtual TemplatePropertyTimelineState<TransFunc1DKeys*>* clone() const;
};

class VRN_CORE_API TransFunc2DPrimitivesPropertyTimelineState : public TemplatePropertyTimelineState<TransFunc2DPrimitives*> {
public:
    TransFunc2DPrimitivesPropertyTimelineState(PropertyKeyValue<TransFunc2DPrimitives*>* kv);
    virtual ~TransFunc2DPrimitivesPropertyTimelineState();
protected:
    friend class TemplatePropertyTimeline<TransFunc2DPrimitives*>;
    friend class XmlDeserializer;
    TransFunc2DPrimitivesPropertyTimelineState();
    const PropertyKeyValue<TransFunc2DPrimitives*>* newKeyValue(float time);
    bool changeValueOfKeyValue(TransFunc2DPrimitives* value, const PropertyKeyValue<TransFunc2DPrimitives*>* keyvalue);
    DeleteKeyValueReturn deleteKeyValue(const PropertyKeyValue<TransFunc2DPrimitives*>* keyvalue);
    /// TransFuncBase* const is const T (not const TransFuncBase*)
    virtual TransFunc2DPrimitives* const getPropertyAt(float time);
    virtual TemplatePropertyTimelineState<TransFunc2DPrimitives*>* clone() const;
};

template <class T>
TemplatePropertyTimelineState<T>::TemplatePropertyTimelineState() {
}

template <class T>
TemplatePropertyTimelineState<T>::~TemplatePropertyTimelineState() {
    // delete all saved keyvalues and the interpolationfunctions between them
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    for (it = values_.begin(); it != values_.end(); ++it) {
        if ((*it).second->getFollowingInterpolationFunction())
            delete ((*it).second->getFollowingInterpolationFunction());

        delete ((*it).second);
    }
    values_.clear();
}

template <class T>
float TemplatePropertyTimelineState<T>::round(float x) {
    //return floor(x * 10000.f) / 10000.f;
    return x;
}

template <class T>
TemplatePropertyTimelineState<T>::TemplatePropertyTimelineState(PropertyKeyValue<T>* kv) {
    // creates the class with one given standard-keyvalue
    values_.insert(std::pair<float,PropertyKeyValue<T>*>(kv->getTime(),kv));
}

template <class T>
const std::map<float,PropertyKeyValue<T>*>& TemplatePropertyTimelineState<T>::getKeyValues() const {
    return values_;
}

template <class T>
const PropertyKeyValue<T>* TemplatePropertyTimelineState<T>::newKeyValue(float time) {
    time = round(time);
    // if there already exists a keyvalue at the given point of time -> do nothing
    if (values_.find(time) != values_.end())
        return 0;

    // calculates the value of the property at the given time to set this value into the new keyvalue
    T value = getPropertyAt(time);
    const InterpolationFunction<T>* oldFunc = getInterpolationFunctionAt(time);

    // create the new keyvalue
    PropertyKeyValue<T>* kv = new PropertyKeyValue<T>(value,time);

    // insert the new keyvalue into the map
    values_.insert(std::pair<float,PropertyKeyValue<T>*>(time,kv));

    // adapt the interpolationfunctions between the keyvalues
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(time);

    // if new value is the first value:
    if (it == values_.begin()) {
        // if new value is the only value:
        it++;
        if (it != values_.end()) {
            InterpolationFunctionFactory fac;
            InterpolationFunction<T>* func = fac.getDefaultFunction<T>()->create();
            //InterpolationFunction<T>* func = new InterpolationFunction<T>();
            (*it).second->setForegoingInterpolationFunction(func);
            it--;
            (*it).second->setFollowingInterpolationFunction(func);
        }
    }
    else {
        it++;
        // if the new value is the last one
        if (it == values_.end()) {
            it--;
            InterpolationFunctionFactory fac;
            InterpolationFunction<T>* func = fac.getDefaultFunction<T>()->create();
            //InterpolationFunction<T>* func = new InterpolationFunction<T>();
            (*it).second->setForegoingInterpolationFunction(func);
            it--;
            (*it).second->setFollowingInterpolationFunction(func);
            it++;
        }
        else {
            // if the value is somewhere in the middle of the timeline
            InterpolationFunction<T>* func1 = oldFunc->create();
            InterpolationFunction<T>* func2 = oldFunc->create();
            (*it).second->setForegoingInterpolationFunction(func2);
            it--;
            (*it).second->setFollowingInterpolationFunction(func2);
            (*it).second->setForegoingInterpolationFunction(func1);
            it--;
            delete (*it).second->getFollowingInterpolationFunction();
            (*it).second->setFollowingInterpolationFunction(func1);
            it++;
        }
    }
    // returns the pointer to the new keyvalue
    return kv;
}

template <class T>
bool TemplatePropertyTimelineState<T>::changeValueOfKeyValue(T value, const PropertyKeyValue<T>* keyvalue) {
    const float time = keyvalue->getTime();
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(time);
    if (it != values_.end()) {
        (*it).second->setValue(value);
        return true;
    }
    return false;
}

template <class T>
bool TemplatePropertyTimelineState<T>::changeSmoothnessOfKeyValue(bool smooth, const PropertyKeyValue<T>* keyvalue) {
    const float time = keyvalue->getTime();
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(time);
    if (it != values_.end()) {
        (*it).second->setSmooth(smooth);
        return true;
    }
    return false;
}

template <class T>
ChangeTimeOfKeyValueReturn TemplatePropertyTimelineState<T>::changeTimeOfKeyValue(float time, const PropertyKeyValue<T>* keyvalue) {
    time = round(time);
    const float oldTime = keyvalue->getTime();
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(oldTime);

    // if the given keyvalue is not in the timeline
    if (it == values_.end())
        return KV_NOT_FOUND;

    // if the times are equal -> nothing to do
    if (time == oldTime)
        return KV_EQUAL_TO_OLD;

    // set the new time of the keyvalue
    PropertyKeyValue<T>* keyvalueToChange = (*it).second;
    keyvalueToChange->setTime(time);

    // adapt all the interpolationfunctions between the keyvalues
    if (time > oldTime) {
        it++;
        if ((it == values_.end()) || ((*it).first > time )) {
            values_.erase(oldTime);
            values_.insert(std::pair<float,PropertyKeyValue<T>*>(time,keyvalueToChange));
            return KV_TIME_CHANGED;
        }
        it--;
    }
    if (time < oldTime) {
        if ((it == values_.begin()) || ((*(--it)).first<time)) {
            values_.erase(oldTime);
            values_.insert(std::pair<float,PropertyKeyValue<T>*>(time,keyvalueToChange));
            return KV_TIME_CHANGED;
        }
        else {
            it++;
        }
    }
    // keyvalue has to be insertet at another point
    // delete keyvalue from the old point:
    // if value is the first one
    if (it == values_.begin()) {
        (*it).second->setFollowingInterpolationFunction(0);
        it++;
        delete ((*it).second->getForegoingInterpolationFunction());
        (*it).second->setForegoingInterpolationFunction(0);
        values_.erase(oldTime);
    }
    else {
        // if value is the last one
        it++;
        if (it == values_.end()) {
            it--;
            (*it).second->setForegoingInterpolationFunction(0);
            it--;
            delete ((*it).second->getFollowingInterpolationFunction());
            (*it).second->setFollowingInterpolationFunction(0);
            values_.erase(oldTime);
        }
        else {
            // if value is in the middle
            InterpolationFunctionFactory fac;
            InterpolationFunction<T>* func = fac.getDefaultFunction<T>()->create();
            //InterpolationFunction<T>* func = new InterpolationFunction<T>();
            (*it).second->setForegoingInterpolationFunction(func);
            it--;
            delete (*it).second->getFollowingInterpolationFunction();
            delete (*it).second->getForegoingInterpolationFunction();
            (*it).second->setFollowingInterpolationFunction(0);
            (*it).second->setForegoingInterpolationFunction(0);
            it--;
            (*it).second->setFollowingInterpolationFunction(func);
            values_.erase(oldTime);
        }
    }
        // insert keyvalue at the new point:
            // if there already exists a keyvalue at the new time
    if (values_.find(time) != values_.end()) {
        it = values_.find(time);
        (*it).second->setValue(keyvalueToChange->getValue());
        delete keyvalueToChange;
        return KV_EXISTS_AT_NEW_TIME;
    }

    values_.insert(std::pair<float,PropertyKeyValue<T>*>(time,keyvalueToChange));

    it = values_.find(time);

    // if new value is the first value:
    if (it == values_.begin()) {
        // only do something if there are multiple values
        it++;
        if (it != values_.end()) {
            InterpolationFunctionFactory fac;
            InterpolationFunction<T>* func = fac.getDefaultFunction<T>()->create();
            //InterpolationFunction<T>* func = new InterpolationFunction<T>()->create();
            (*it).second->setForegoingInterpolationFunction(func);
            it--;
            (*it).second->setFollowingInterpolationFunction(func);
        }
    }
    else {
        it++;
        // if the new value is the last one
        if (it == values_.end()) {
            it--;
            InterpolationFunctionFactory fac;
            InterpolationFunction<T>* func = fac.getDefaultFunction<T>()->create();
            //InterpolationFunction<T>* func = new InterpolationFunction<T>();
            (*it).second->setForegoingInterpolationFunction(func);
            it--;
            (*it).second->setFollowingInterpolationFunction(func);
        }
        else {
            // if new value is in the middle
            InterpolationFunctionFactory fac;
            InterpolationFunction<T>* func1 = fac.getDefaultFunction<T>()->create();
            InterpolationFunction<T>* func2 = fac.getDefaultFunction<T>()->create();
            //InterpolationFunction<T>* func1 = new InterpolationFunction<T>();
            //InterpolationFunction<T>* func2 = new InterpolationFunction<T>();
            (*it).second->setForegoingInterpolationFunction(func2);
            it--;
            (*it).second->setFollowingInterpolationFunction(func2);
            (*it).second->setForegoingInterpolationFunction(func1);
            it--;
            delete (*it).second->getFollowingInterpolationFunction();
            (*it).second->setFollowingInterpolationFunction(func1);
        }
    }
    return KV_MOVED_TO_DIFFERENT_INTERVAL;
}

template <class T>
ChangeTimeOfKeyValueReturn TemplatePropertyTimelineState<T>::shiftKeyValue(float time, const PropertyKeyValue<T>* keyvalue) {
    const float oldTime = keyvalue->getTime();
    const float diff = time - oldTime;
    time = round(time);

    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(oldTime);

    // if the given keyvalue is not in the timeline
    if (it == values_.end()) {
        return KV_NOT_FOUND;
    }

    // if the times are equal -> nothing to do
    if (time == oldTime) {
        return KV_EQUAL_TO_OLD;
    }

    // set the new time of the keyvalue
    PropertyKeyValue<T>* keyvalueToChange = (*it).second;

    // adapt all the interpolationfunctions between the keyvalues
    if (time > oldTime) {
        typename std::map<float,PropertyKeyValue<T>*>::reverse_iterator it = values_.rbegin();

        bool end = false;
        while(!end) {
            float tempTime = (*it).first;
            PropertyKeyValue<T>* tempKey = (*it).second;

            values_.erase(tempTime);
            tempTime = round(tempTime + diff);
            tempKey->setTime(tempTime);
            values_.insert(std::pair<float, PropertyKeyValue<T>*>(tempTime, tempKey));

            if(tempKey == keyvalueToChange)
                end = true;

            it++;
        }
    }
    else {
        //typename std::map<float,PropertyKeyValue<T>*>::iterator it = values_.begin();

        //if (((*it).second->getTime() + diff) < 0.0f)
            //return KV_TIME_AFTER_DURATION;

        while(it != values_.end()) {
            float tempTime = (*it).first;
            PropertyKeyValue<T>* tempKey = (*it).second;

            values_.erase(tempTime);
            tempTime = round(tempTime + diff);
            tempKey->setTime(tempTime);
            values_.insert(std::pair<float, PropertyKeyValue<T>*>(tempTime, tempKey));

            it++;
        }
    }

    return KV_TIME_CHANGED;
}

template <class T>
DeleteKeyValueReturn TemplatePropertyTimelineState<T>::deleteKeyValue(const PropertyKeyValue<T>* keyvalue){
    const float time = keyvalue->getTime();
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(time);

    // if wrong parameter do nothing
    if (it == values_.end()) {
        return KV_NOT_THERE;
    }

    // if keyvalue ist the only value do nothing
    if (it == values_.begin()) {
        it++;
        if (it == values_.end())
            return KV_IS_THE_ONLY_ONE;
        it--;
    }

    // if value is the first one
    if (it == values_.begin()) {
        it++;
        delete ((*it).second->getForegoingInterpolationFunction());
        (*it).second->setForegoingInterpolationFunction(0);
        values_.erase(time);
        return KV_DELETED;
    }
    // if value is the last one
    it++;
    if (it == values_.end()){
        it--;
        it--;
        delete ((*it).second->getFollowingInterpolationFunction());
        (*it).second->setFollowingInterpolationFunction(0);
        values_.erase(time);
        return KV_DELETED;
    }
    // if value is in the middle
    InterpolationFunctionFactory fac;
    InterpolationFunction<T>* func = fac.getDefaultFunction<T>()->create();
    //InterpolationFunction<T>* func = new InterpolationFunction<T>();
    (*it).second->setForegoingInterpolationFunction(func);
    it--;
    delete (*it).second->getFollowingInterpolationFunction();
    delete (*it).second->getForegoingInterpolationFunction();
    it--;
    (*it).second->setFollowingInterpolationFunction(func);
    values_.erase(time);
    return KV_DELETED;
}

template <class T>
void TemplatePropertyTimelineState<T>::setInterpolationFunctionBefore(InterpolationFunction<T>* func,PropertyKeyValue<T>* keyvalue) {
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(keyvalue->getTime());
    if (it == values_.end())
        return;

    if (it == values_.begin())
        // there cannot be an interpolationfunction before the first keyvalue
        return;

    (*it).second->setForegoingInterpolationFunction(func);
    it--;
    (*it).second->setFollowingInterpolationFunction(func);
}

template <class T>
void TemplatePropertyTimelineState<T>::setInterpolationFunctionAfter(InterpolationFunction<T>* func,PropertyKeyValue<T>* keyvalue) {
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(keyvalue->getTime());
    if (it == values_.end())
        return;
    it++;
    if (it == values_.end())
        // there cannot be an interpolationfunction after the last keyvalue
        return;

    (*it).second->setForegoingInterpolationFunction(func);
    it--;
    (*it).second->setFollowingInterpolationFunction(func);
}

template <class T>
const InterpolationFunction<T>* TemplatePropertyTimelineState<T>::getInterpolationFunctionAt(float time) {
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(time);

    // if the given point of time is exactly on a keyvalue
    if (it != values_.end())
        return 0;

    // if the given point of time is earlier than the first keyalue
    it = values_.upper_bound(time);
    if (it == values_.begin())
        return 0;

    // if the given point of time is later than the last keyalue
    if (it == values_.end())
        return 0;

    // if the given point of time lies between two keyvalues -> call the corresponding interpolationfunction
    it--;
    return (*it).second->getFollowingInterpolationFunction();
}

template <class T>
const T TemplatePropertyTimelineState<T>::getPropertyAt(float time) {
    typename std::map<float,PropertyKeyValue<T>*>::iterator it;
    it = values_.find(time);

    // if the given point of time is exactly on a keyvalue
    if (it != values_.end())
        return ((*it).second->getValue());

    // if the given point of time is earlier than the first keyalue
    it = values_.upper_bound(time);
    if (it == values_.begin())
        return ((*it).second->getValue());

    // if the given point of time is later than the last keyalue
    if (it == values_.end()) {
        it--;
        return ((*it).second->getValue());
    }
    // if the given point of time lies between two keyvalues -> call the corresponding interpolationfunction
    typename std::map<float,PropertyKeyValue<T>*>::iterator it2;
    it2 = it;
    it--;

    const InterpolationFunction<T>* func = (*it).second->getFollowingInterpolationFunction();
    const MultiPointInterpolationFunction<T>* multifunc = dynamic_cast<const MultiPointInterpolationFunction<T>*>(func);
    if (multifunc) {
        // call a function with multiple points
        // create vector of the used keyvalues
        std::vector<PropertyKeyValue<T>*> keys;
        // search for the first value in the multi-point interval
        while ((it!=values_.begin()) && (it->second->isSmoothed()))
            it--;

        do {
            keys.push_back((*it).second->clone());
            it++;
        }
        while ((it != values_.end()) && it->second->isSmoothed());

        if (it != values_.end())
            keys.push_back((*it).second->clone());

        // interpolate value
        T returnvalue = multifunc->interpolate(keys,time);

        // delete all copied keys
        typename std::vector<PropertyKeyValue<T>*>::const_iterator delIt;
        for (delIt = keys.begin(); delIt != keys.end(); ++delIt)
            delete (*delIt);
        keys.clear();

        // return
        return returnvalue;
    }
    else
        return func->interpolate( it->second->getValue(), it2->second->getValue(), (time- (it->first)) / ((it2->first) - (it->first)));
}

template <class T>
TemplatePropertyTimelineState<T>* TemplatePropertyTimelineState<T>::clone() const {
    TemplatePropertyTimelineState<T>* timeline = new TemplatePropertyTimelineState<T>();

    // copy all keyvalues
    typename std::map<float,PropertyKeyValue<T>*>::const_iterator it;
    for (it = values_.begin(); it != values_.end(); ++it)
        timeline->values_.insert(std::pair<float,PropertyKeyValue<T>*>(it->first, it->second->clone()));

    // copy all interpolationfunctions between the keyvalues
    typename std::map<float,PropertyKeyValue<T>*>::const_iterator it2;
    it2 = timeline->values_.begin();
    for (it = values_.begin(); it != values_.end();) {
        const InterpolationFunction<T>* func;
        InterpolationFunction<T>* func2;
        func = (*it).second->getFollowingInterpolationFunction();
        if (!func) {
            it2++;
            it++;
            continue;
        }
        func2 = func->create();
        (*it2).second->setFollowingInterpolationFunction(func2);
        it2++;
        it++;
        (*it2).second->setForegoingInterpolationFunction(func2);
    }

    return timeline;
}

template <class T>
void TemplatePropertyTimelineState<T>::setDuration(float duration) {
    std::map<float,PropertyKeyValue<T>*> newValues;

    typename std::map<float,PropertyKeyValue<T>*>::const_iterator it = values_.begin();

    // if first keyvalue is too late -> set it to time=0
    if ((*it).first > duration) {
        (*it).second->setTime(0);
        this->changeTimeOfKeyValue(0,(*it).second);
    }

    it = values_.begin();

    // copy all keyvalues with time<=duration
    while((it != values_.end()) && (it->first <= duration)) {
        newValues.insert(std::pair<float,PropertyKeyValue<T>*>(it->first, it->second->clone()));
        it++;
    }

    // copy all interpolationfunctions
    typename std::map<float,PropertyKeyValue<T>*>::const_iterator it2;
    it2 = values_.begin();
    for (it = newValues.begin(); it != newValues.end();) {
        const InterpolationFunction<T>* func;
        InterpolationFunction<T>* func2;

        func = it2->second->getFollowingInterpolationFunction();
        if (!func) {
            it2++;
            it++;
            continue;
        }
        func2 = func->create();
        it++;
        it2++;
        if (it != newValues.end()) {
            it--;
            it->second->setFollowingInterpolationFunction(func2);
            it++;
            it->second->setForegoingInterpolationFunction(func2);
        }
        else {
            delete func2;
            it--;
            it->second->setFollowingInterpolationFunction(0);
            it++;
        }
    }

    // delete all old keyvalues
    for (it = values_.begin(); it != values_.end(); ++it) {
        if (it->second->getFollowingInterpolationFunction())
            delete (it->second->getFollowingInterpolationFunction());

        delete (it->second);
    }
    values_.clear();

    values_ = newValues;
}

template <class T>
void TemplatePropertyTimelineState<T>::serialize(Serializer& s) const {
    s.serialize("values", values_);
}

template <class T>
void TemplatePropertyTimelineState<T>::deserialize(Deserializer& s) {
    s.deserialize("values", values_);
}

} // namespace voreen

#endif
