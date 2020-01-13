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

#ifndef VRN_TEMPLATEPROPERTYTIMELINE_H
#define VRN_TEMPLATEPROPERTYTIMELINE_H

#include <list>
#include <deque>
#include <map>
#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include "voreen/core/animation/propertykeyvalue.h"
#include "voreen/core/animation/interpolationfunction.h"
#include "voreen/core/animation/propertytimeline.h"
#include "voreen/core/animation/templatepropertytimelinestate.h"
#include "voreen/core/animation/timelineobserver.h"
#include "voreen/core/animation/animation.h"

#include "voreen/core/io/serialization/serialization.h"
#include "tgt/bounds.h"

namespace voreen {

class Animation;
class TransFunc1DKeys;
class TransFunc2DPrimitives;

/**
 * This class organizes all animationsettings of one property of the type T.
 */
template <class T>
class TemplatePropertyTimeline : public PropertyTimeline, public Observable<TimelineObserver> {
public:
    /**
     * Parameter is the property object which should be animated.
     */
    TemplatePropertyTimeline(TemplateProperty<T>* prop);

    virtual ~TemplatePropertyTimeline();

    /// returns true if the timeline is empty
    bool isEmpty();

    virtual bool isCompatibleWith(const Property* p) const;

    /**
     * Returns a sorted vector of all keyvalues of the timeline.
     */
    virtual const std::map<float,PropertyKeyValue<T>*> getKeyValues() const;

    /**
     * Returns a pointer to the corresponding property for the gui.
     * With this pointer it is possible to determine the subtype of the property.
     */
    virtual const TemplateProperty<T>* getCorrespondingProperty() const;

    /**
     * Adds a new keyvalue at a certain time:
     * If there already exists a keyvalue at the given time, this function does nothing.
     *
     * The interplationfunction at the time of the new keyvalue will be deleted and replaced by two default
     * interpolations before and after the new keyvalue.
     * Returns the pointer to the new keyvalue.
     * If there already existed a keyvalue at the given time a nullpointer is returned.
     * If the given time is after duration time a nullpointer is returned.
     */
    virtual const PropertyKeyValue<T>* newKeyValue(float time);

    /**
     * Changes the value of a keyvalue:
     * returns true, if the value was changed,
     * returns false, if the given keyvalue couldn't be found in the timeline.
     */
    virtual bool changeValueOfKeyValue(T value, const PropertyKeyValue<T>* keyvalue);

    /**
     * Changes the smooth-variable of a keyvalue:
     * returns true, if the value was changed,
     * returns false, if the given keyvalue couldn't be found in the timeline.
     */
    virtual bool changeSmoothnessOfKeyValue(bool smooth, const PropertyKeyValue<T>* keyvalue);

    /**
     * Implements method of abstract class PropertyTimeline.
     * Checks if the key value is suitable for this template class and if so, calls changeTimeOfKeyValue.
     */
    virtual void changeTimeOfKeyValueAbstract(float time, const PropertyKeyValueBase* key);

    /**
     * Changes the time of a keyvalue:
     * If there is already a keyvalue at the given time, this function deletes the given keyvalue
     * and sets the value of the already existing keyvalue to the deleted value.
     * If the order of the keyvalues changes by setting the new time, the function reacts as if
     * the given keyvalue was deleted and a new one was added at the given time with the same value as the deleted one.
     * Returns true if the given time is valid and the change is done.
     *
     * returns KV_NOT_FOUND                   if the given keyvalue couldn't be found in the timeline
     * returns KV_EQUAL_TO_OLD                if the new time is equal to the old one
     * returns KV_TIME_CHANGED                if the new time is in the same interval between the same other keyvalues
     * returns KV_EXISTS_AT_NEW_TIME          if there already exists a new keyvalue at the new time
     * returns KV_MOVED_TO_DIFFERENT_INTERVAL if the new time is in another interval between two other keyvalues.
     */
    virtual ChangeTimeOfKeyValueReturn changeTimeOfKeyValue(float time, const PropertyKeyValue<T>* keyvalue);

    /**
     * Implements method of abstract class PropertyTimeline.
     * Checks if the key value is suitable for this template class and if so, calls shiftKeyValue.
     */
    virtual void shiftKeyValueAbstract(float time, const PropertyKeyValueBase* key);

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
     * returns KV_IS_THE_ONLY_ONE if teh given keyvalue is the only value in the timeline and mustn#t be deleted
     * returns KV_DELETED         if the keyvalue was succesful deleted
     */
    virtual DeleteKeyValueReturn deleteKeyValue(const PropertyKeyValue<T>* keyvalue);

    /**
     * Implements method of abstract class PropertyTimeline.
     * Checks if the key value is suitable for this template class and if so, calls deleteKeyValue.
     */
    virtual DeleteKeyValueReturn removeKeyValueAbstract(const PropertyKeyValueBase* key);

    /**
     * Sets the interpolationfunction between the given keyvalue and its successor.
     */
    virtual void setInterpolationFunctionAfter(InterpolationFunction<T>* func, PropertyKeyValue<T>* keyvalue);
    /**
     * Sets the interpolationfunction between the given keyvalue and its predecessor.
     */
    virtual void setInterpolationFunctionBefore(InterpolationFunction<T>* func, PropertyKeyValue<T>* keyvalue);

    /**
     * This function deletes all keyvalues of the timeline and takes the current setting of the processornetwork
     * as the standard value at time 0.
     */
    virtual void resetTimeline();

    /**
     * This function returns true if the corresponding property is linked with another one.
     */
    virtual bool propertyIsLinked() const;

    /**
     * Sets the current setting as a keyvalue at the given time if it differs from the animated value.
     */
    virtual void setCurrentSettingAsKeyvalue(float time, bool forceKeyValue);

    /**
     * Returns a copy of the current timelinestate.
     */
    virtual TemplatePropertyTimelineState<T>* getCurrentTimelineState() const ;

    /**
     * This function registers the given animation class as an utils/observer.
     * Everytime there is a change in the timeline the observer will be contacted.
     */
    virtual void registerUndoObserver(Animation* observer);

    /**
     * Returns the name of the animated property.
     */
    virtual std::string getPropertyName() const ;

    /**
     * Sets the interactionmode of the corresponding property.
     */
    virtual void setInteractionMode(bool interactionmode,void* source);

    /**
     * Calculates the animated value to a certain point of time.
     */
    const T getPropertyAt(float time);

    /**
     * Creates a new undo-state.
     */
    virtual void setNewUndoState();

    /**
     *  Sets a bool if this property should be animated on rendering.
     */
    virtual void setActiveOnRendering(bool activeOnRendering);

    /**
     * Returns if the property is animated on rendering.
     */
    virtual bool getActiveOnRendering() const;

    /**
     * Returns true if the timeline was changed.
     */
    virtual bool isChanged() const ;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    virtual Property* getProperty() const;

    virtual void setProperty(Property* p);

    virtual AbstractSerializable* create() const = 0;

    virtual std::string getClassName() const = 0;

    /**
     * Returns a sorted vector of the timesteps for all key values.
     */
    virtual std::vector<float> getTimesteps() const;

    /**
     * Returns a sorted vector of all keyvalues of the timeline, but with the abstract key value class (for generic use with propertytimeline class)
     */
    virtual const std::map<float,PropertyKeyValueBase*> getKeyValuesAbstract() const;

    virtual const PropertyKeyValueBase* newKeyValueAbstract(float time);

    virtual bool setCurrentSettingToKeyValue(PropertyKeyValueBase* key);

protected:
    friend class PropertyTimelineFactory;
    friend class XmlDeserializer;
    friend class TimelineObserver;
    friend class Animation;

    /**
     * Default constructor.
     */
    TemplatePropertyTimeline();

    /**
     * This function calculates the propertyvalue at the given time and changes the corresponding property object to that value.
     */
    virtual void renderAt(float time);

    /**
     * Called from animation-class.
     * All keyvalues after the durationtime will be deleted.
     */
    virtual void setDuration(float duration);

    /**
     * Pointer to the corresponding property.
     */
    TemplateProperty<T>* property_;

    bool activeOnRendering_;

    bool timelineChanged_;

    /**
     * Calculates the animated value to a certain point of time.
     */
    T privateGetPropertyAt(float time);

    /**
     * This function provides an undo-functionality for this timeline.
     * It should only be called by the undoobserver which manages all undos and redos of the whole animation.
     * If the undo-function is called another timelinestate is set as current, this means that all pointers to the keyvalues of this timeline get invalid.
     */
    virtual void undo();

    /**
     * This function provides a redo-functionality for this timeline.
     * It should only be called by the undoObserver which manages all undos and redos of the whole animation.
     * If the undo-function is called another timelinestate is set as current, this means that all pointers to the keyvalues of this timeline get invalid.
     */
    virtual void redo();

    /*
     * This function deletes all redo-states (stack lastundos_).
     */
    virtual void clearRedoStates();

    /**
     * This function deletes all states (stack lastundos_ and lastchanges_)
     */
    virtual void clearAllStates();

    /**
     * This function removes the oldest undo-state.
     * If the number of possible undo-steps is decreased, this function adapts the related stack lastchanges_.
     */
    virtual void removeOldestUndoState();

    /**
     * Pointer to the current timelinestate.
     */
    TemplatePropertyTimelineState<T>* timeline_;

    /**
     * Deques of all undo and redo timelinestates.
     */
    std::deque<TemplatePropertyTimelineState<T>*> lastChanges_;
    std::deque<TemplatePropertyTimelineState<T>*> lastUndos_;

    /**
     * Pointer to the animationclass as undoutils/observer.
     */
    Animation* undoObserver_;
};

template <class T>
TemplatePropertyTimeline<T>::TemplatePropertyTimeline(TemplateProperty<T>* prop)
: property_(prop)
, activeOnRendering_(true)
, timelineChanged_(false)
, undoObserver_(0)
{
    duration_ = 60.f * 15.f;

    timeline_ = new TemplatePropertyTimelineState<T>(new PropertyKeyValue<T>(property_->get(),0));
}

template <>
TemplatePropertyTimeline<TransFunc1DKeys*>::TemplatePropertyTimeline(TemplateProperty<TransFunc1DKeys*>* prop);
template <>
TemplatePropertyTimeline<TransFunc2DPrimitives*>::TemplatePropertyTimeline(TemplateProperty<TransFunc2DPrimitives*>* prop);

template <class T>
TemplatePropertyTimeline<T>::TemplatePropertyTimeline()
: property_(0)
, activeOnRendering_(true)
, timelineChanged_(false)
, undoObserver_(0)
{
    duration_ = 60.f * 15.f;

    timeline_ = new TemplatePropertyTimelineState<T>();
}

template <>
TemplatePropertyTimeline<TransFunc1DKeys*>::TemplatePropertyTimeline();
template <>
TemplatePropertyTimeline<TransFunc2DPrimitives*>::TemplatePropertyTimeline();

template <class T>
TemplatePropertyTimeline<T>::~TemplatePropertyTimeline() {
    delete timeline_;

    TemplatePropertyTimelineState<T>* temp;
    while (lastChanges_.size()) {
        temp = lastChanges_.back();
        lastChanges_.pop_back();
        delete temp;
    }
    while (lastUndos_.size()) {
        temp = lastUndos_.back();
        lastUndos_.pop_back();
        delete temp;
    }
}

template <class T>
bool TemplatePropertyTimeline<T>::isCompatibleWith(const Property* p) const {
    return dynamic_cast<const TemplateProperty<T>*>(p);
}

template <class T>
bool TemplatePropertyTimeline<T>::isEmpty() {
    if (timeline_->getKeyValues().size() > 1)
        return false;
    else
        return true;
}

template <class T>
void TemplatePropertyTimeline<T>::resetTimeline() {
    timelineChanged_ = true;
    lastChanges_.push_back(timeline_);
    undoObserver_->animationChanged(this);

    timeline_ = new TemplatePropertyTimelineState<T>(new PropertyKeyValue<T>(property_->get(),0));

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();
}

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::resetTimeline();
template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::resetTimeline();

template <class T>
const T TemplatePropertyTimeline<T>::getPropertyAt(float time) {
    return privateGetPropertyAt(time);
}

template <class T>
T TemplatePropertyTimeline<T>::privateGetPropertyAt(float time) {
    return timeline_->getPropertyAt(time);
}

template <>
TransFunc1DKeys* TemplatePropertyTimeline<TransFunc1DKeys*>::privateGetPropertyAt(float time);
template <>
TransFunc2DPrimitives* TemplatePropertyTimeline<TransFunc2DPrimitives*>::privateGetPropertyAt(float time);

template <class T>
const std::map<float,PropertyKeyValue<T>*> TemplatePropertyTimeline<T>::getKeyValues() const{
    return timeline_->getKeyValues();
}

template <class T>
std::vector<float> TemplatePropertyTimeline<T>::getTimesteps() const {
    //get the key values
    const std::map<float, PropertyKeyValue<T>*> keys = timeline_->getKeyValues();
    std::vector<float> timesteps;
    //iterate over map and put keys (ie. float values) into vector
    for (typename std::map<float,PropertyKeyValue<T>*>::const_iterator i = keys.begin(); i != keys.end(); ++i)
        timesteps.push_back(i->first);
    return timesteps;
}

template <class T>
const std::map<float,PropertyKeyValueBase*> TemplatePropertyTimeline<T>::getKeyValuesAbstract() const {
    //get the key values
    const std::map<float, PropertyKeyValue<T>*> keys = timeline_->getKeyValues();
    //insert into new map with abstract type
    std::map<float, PropertyKeyValueBase*> abstractKeys;
    for (typename std::map<float,PropertyKeyValue<T>*>::const_iterator i = keys.begin(); i != keys.end(); ++i) {
        float time = i->first;
        PropertyKeyValueBase* key = i->second;
        abstractKeys.insert(std::make_pair(time, key));
    }
    return abstractKeys;
}

template <class T>
const PropertyKeyValueBase* TemplatePropertyTimeline<T>::newKeyValueAbstract(float time) {
    return newKeyValue(time);
}

template <class T>
const PropertyKeyValue<T>* TemplatePropertyTimeline<T>::newKeyValue(float time) {
    //time = floor(time * 10000.f) / 10000.f;

    if (time > duration_) {
        timelineChanged_ = false;
        return 0;
    }

    lastChanges_.push_back(timeline_->clone());
    undoObserver_->animationChanged(this);

    const PropertyKeyValue<T>* kv = timeline_->newKeyValue(time);

    //setInterpolationFunctionAfter(, kv); //TODO
    //setInterpolationFunctionBefore(, kv);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();

    return kv;
}

template <class T>
bool TemplatePropertyTimeline<T>::setCurrentSettingToKeyValue(PropertyKeyValueBase* key) {
    timelineChanged_ = true;

    PropertyKeyValue<T>* tKey = dynamic_cast<PropertyKeyValue<T>*>(key);
    if (!tKey)
        return false;

    T value = property_->get();
    std::string errorMsg;
    if (!(property_->isValidValue(value, errorMsg))) {
        LWARNINGC("voreen.TemplatePropertyTimeline", "Invalid property value");
        return false;
    }
    T animatedValue = tKey->getValue();
    if (value != animatedValue)
        return changeValueOfKeyValue(value,tKey);
    else
        return true;
}

/*template <>
void TemplatePropertyTimeline<TransFuncBase*>::setCurrentSettingToKeyvalue(PropertyKeyValueBase* key);*/

template <class T>
bool TemplatePropertyTimeline<T>::changeValueOfKeyValue(T value, const PropertyKeyValue<T>* keyvalue) {
    timelineChanged_ = true;
    std::string errorMsg;
    if (!(property_->isValidValue(value, errorMsg))) {
        LWARNINGC("voreen.TemplatePropertyTimeline", "Invalid property value");
        return false;
    }
    bool temp = timeline_->changeValueOfKeyValue(value, keyvalue);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();

    return temp;
}

template <class T>
bool TemplatePropertyTimeline<T>::changeSmoothnessOfKeyValue(bool smooth, const PropertyKeyValue<T>* keyvalue) {
    timelineChanged_ = true;
    bool temp = timeline_->changeSmoothnessOfKeyValue(smooth, keyvalue);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();

    return temp;
}

template <class T>
void TemplatePropertyTimeline<T>::changeTimeOfKeyValueAbstract(float time, const PropertyKeyValueBase* key) {
    if (const PropertyKeyValue<T>* k = dynamic_cast<const PropertyKeyValue<T>*>(key))
        changeTimeOfKeyValue(time, k);
}

template <class T>
ChangeTimeOfKeyValueReturn TemplatePropertyTimeline<T>::changeTimeOfKeyValue(float time, const PropertyKeyValue<T>* keyvalue) {
    if (time > duration_) {
        timelineChanged_ = false;
        return KV_TIME_AFTER_DURATION;
    }

    //time = floor(time * 10000.f) / 10000.f;
    timelineChanged_ = true;

    ChangeTimeOfKeyValueReturn temp = timeline_->changeTimeOfKeyValue(time,keyvalue);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();

    return temp;
}

template <class T>
void TemplatePropertyTimeline<T>::shiftKeyValueAbstract(float time, const PropertyKeyValueBase* key) {
    if (const PropertyKeyValue<T>* k = dynamic_cast<const PropertyKeyValue<T>*>(key))
        shiftKeyValue(time, k);
}

template <class T>
ChangeTimeOfKeyValueReturn TemplatePropertyTimeline<T>::shiftKeyValue(float time, const PropertyKeyValue<T>* keyvalue) {
    //time = floor(time * 10000.f) / 10000.f;
    const float oldTime = keyvalue->getTime();
    const float diff = time - oldTime;
    timelineChanged_ = true;

    const std::map<float,PropertyKeyValue<T>*>& values = timeline_->getKeyValues();

    {
        typename std::map<float,PropertyKeyValue<T>*>::const_iterator it = values.find(oldTime);

        // if the given keyvalue is not in the timeline
        if (it == values.end())
            return KV_NOT_FOUND;

        // if the times are equal -> nothing to do
        if (time == oldTime)
            return KV_EQUAL_TO_OLD;

        if(time > oldTime) {
            typename std::map<float,PropertyKeyValue<T>*>::const_reverse_iterator last = values.rbegin();

            if(((*last).second->getTime() + diff) > duration_)
                return KV_TIME_AFTER_DURATION;
        }
        else {
            if(it != values.begin()) {
                it--; // check if movement would colide with predecessor

                if(((*it).second->getTime() > time) )
                    return KV_TIME_AFTER_DURATION;
            }
        }
    }

    ChangeTimeOfKeyValueReturn temp = timeline_->shiftKeyValue(time,keyvalue);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();

    return temp;
}

template <class T>
DeleteKeyValueReturn TemplatePropertyTimeline<T>::removeKeyValueAbstract(const PropertyKeyValueBase* key) {
    if (const PropertyKeyValue<T>* k = dynamic_cast<const PropertyKeyValue<T>*>(key))
        return deleteKeyValue(k);
    else
        return KV_NOT_THERE;
}

template <class T>
DeleteKeyValueReturn TemplatePropertyTimeline<T>::deleteKeyValue(const PropertyKeyValue<T>* keyvalue) {
    timelineChanged_ = true;
    lastChanges_.push_back(timeline_->clone());
    undoObserver_->animationChanged(this);

    DeleteKeyValueReturn temp = timeline_->deleteKeyValue(keyvalue);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();

    return temp;
}

template <class T>
void TemplatePropertyTimeline<T>::setInterpolationFunctionBefore(InterpolationFunction<T>* func,PropertyKeyValue<T>* keyvalue) {
    timelineChanged_ = true;
    lastChanges_.push_back(timeline_->clone());
    undoObserver_->animationChanged(this);

    timeline_->setInterpolationFunctionBefore(func,keyvalue);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();
}

template <class T>
void TemplatePropertyTimeline<T>::setInterpolationFunctionAfter(InterpolationFunction<T>* func,PropertyKeyValue<T>* keyvalue) {
    timelineChanged_ = true;
    lastChanges_.push_back(timeline_->clone());
    undoObserver_->animationChanged(this);

    timeline_->setInterpolationFunctionAfter(func,keyvalue);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();
}

template <class T>
void TemplatePropertyTimeline<T>::renderAt(float time) {
    if (!activeOnRendering_)
        return;

    property_->set(getPropertyAt(time));
}

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::renderAt(float time);
template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::renderAt(float time);

template <class T>
bool TemplatePropertyTimeline<T>::propertyIsLinked() const {
    if (property_->getLinks().size() > 0)
        return true;
    else
        return false;
}

template <class T>
const TemplateProperty<T>* TemplatePropertyTimeline<T>::getCorrespondingProperty() const {
    return property_;
}

template <class T>
void TemplatePropertyTimeline<T>::registerUndoObserver(Animation* observer) {
    undoObserver_ = observer;
}

template <class T>
std::string TemplatePropertyTimeline<T>::getPropertyName() const {
    return property_->getGuiName();
}

template <class T>
void TemplatePropertyTimeline<T>::undo() {
    lastUndos_.push_back(timeline_);
    timeline_ = lastChanges_.back();
    lastChanges_.pop_back();

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();
}

template <class T>
void TemplatePropertyTimeline<T>::redo() {
    lastChanges_.push_back(timeline_);
    timeline_ = lastUndos_.back();
    lastUndos_.pop_back();

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();
}

template <class T>
void TemplatePropertyTimeline<T>::clearRedoStates() {
    TemplatePropertyTimelineState<T>* temp;
    while (lastUndos_.size() > 0) {
        temp = lastUndos_.back();
        lastUndos_.pop_back();
        delete temp;
    }
}

template <class T>
void TemplatePropertyTimeline<T>::clearAllStates() {
    TemplatePropertyTimelineState<T>* temp;
    while (lastUndos_.size() > 0) {
        temp = lastUndos_.back();
        lastUndos_.pop_back();
        delete temp;
    }
    while (lastChanges_.size() > 0) {
        temp = lastChanges_.back();
        lastChanges_.pop_back();
        delete temp;
    }
}

template <class T>
void TemplatePropertyTimeline<T>::removeOldestUndoState() {
    TemplatePropertyTimelineState<T>* temp;
    temp = lastChanges_.front();
    lastChanges_.pop_front();
    delete temp;
}

template <class T>
void TemplatePropertyTimeline<T>::setInteractionMode(bool interactionmode,void* source) {
    QualityMode.requestQualityMode((interactionmode ? VoreenQualityMode::RQ_INTERACTIVE : VoreenQualityMode::RQ_DEFAULT),source);
}

template <class T>
void TemplatePropertyTimeline<T>::setNewUndoState() {
    lastChanges_.push_back(timeline_->clone());
    undoObserver_->animationChanged(this);
}

template <class T>
bool TemplatePropertyTimeline<T>::getActiveOnRendering() const {
    return activeOnRendering_;
}

template <class T>
void TemplatePropertyTimeline<T>::setActiveOnRendering(bool activeOnRendering) {
    activeOnRendering_ = activeOnRendering;
}

template <class T>
bool TemplatePropertyTimeline<T>::isChanged() const {
    return timelineChanged_;
}

template <class T>
TemplatePropertyTimelineState<T>* TemplatePropertyTimeline<T>::getCurrentTimelineState() const {
    return timeline_->clone();
}

template <class T>
void TemplatePropertyTimeline<T>::setDuration(float duration) {
    duration_ = duration;

    timeline_->setDuration(duration);

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it)
        (*it)->timelineChanged();
}

template <class T>
void TemplatePropertyTimeline<T>::setCurrentSettingAsKeyvalue(float time, bool forceKeyValue) {
    //time = floor(time*10000)/10000;
    T value = property_->get();
    T animatedValue = getPropertyAt(time);
    if (value != animatedValue) {
        const PropertyKeyValue<T>* kv = timeline_->newKeyValue(time);
        if (!kv) {
            kv = new PropertyKeyValue<T>(value,time);
        }
        this->changeValueOfKeyValue(value,kv);
    }
    else {
        if (forceKeyValue){
            newKeyValue(time);
        }
    }
}

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::setCurrentSettingAsKeyvalue(float time, bool forceKeyValue);
template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::setCurrentSettingAsKeyvalue(float time, bool forceKeyValue);

template <class T>
void TemplatePropertyTimeline<T>::serialize(Serializer& s) const {
    s.serialize("activeOnRendering", activeOnRendering_);
    s.serialize("propertyOwner", property_->getOwner());
    s.serialize("propertyId", property_->getID());
    s.serialize("duration", duration_);
    s.serialize("timeline", timeline_);
    s.serialize("tlchanged", timelineChanged_);
}

template <class T>
void TemplatePropertyTimeline<T>::deserialize(Deserializer& s) {
    s.deserialize("activeOnRendering", activeOnRendering_);
    PropertyOwner* propertyOwner = 0;
    s.deserialize("propertyOwner", propertyOwner);
    std::string propertyId;
    s.deserialize("propertyId", propertyId);
    property_ = dynamic_cast<TemplateProperty<T>*>(propertyOwner->getProperty(propertyId));
    if (!property_)
        LWARNINGC("TemplatePropertyTimeline", "Property defined in animation timeline does not exist: "
        << propertyOwner->getID() << "::" << propertyId);

    s.deserialize("duration", duration_);
    s.deserialize("timeline", timeline_);
    s.deserialize("tlchanged", timelineChanged_);
}

template <class T>
Property* TemplatePropertyTimeline<T>::getProperty() const {
    return property_;
}

template <class T>
void TemplatePropertyTimeline<T>::setProperty(Property* p) {
    TemplateProperty<T>* tp = dynamic_cast<TemplateProperty<T>*>(p);
    if(tp) {
        property_ = tp;
        if(timeline_->getKeyValues().size() == 0) {
            delete timeline_;
            timeline_ = 0;
            timeline_ = new TemplatePropertyTimelineState<T>(new PropertyKeyValue<T>(property_->get(),0));
        }
    }
    else {
        LERRORC("voreen.TemplatePropertyTimeline", "Property type mismatch!");
    }
}

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::setProperty(Property* p);
template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::setProperty(Property* p);

class VRN_CORE_API PropertyTimelineFloat : public TemplatePropertyTimeline<float> {
public:
    PropertyTimelineFloat() : TemplatePropertyTimeline<float>() {}
    PropertyTimelineFloat(TemplateProperty<float>* prop) : TemplatePropertyTimeline<float>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineFloat(); }

    virtual std::string getClassName() const { return "PropertyTimelineFloat"; }
};

class VRN_CORE_API PropertyTimelineInt : public TemplatePropertyTimeline<int> {
public:
    PropertyTimelineInt() : TemplatePropertyTimeline<int>() {}
    PropertyTimelineInt(TemplateProperty<int>* prop) : TemplatePropertyTimeline<int>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineInt(); }

    virtual std::string getClassName() const { return "PropertyTimelineInt"; }
};

class VRN_CORE_API PropertyTimelineBool : public TemplatePropertyTimeline<bool> {
public:
    PropertyTimelineBool()  : TemplatePropertyTimeline<bool>() {}
    PropertyTimelineBool(TemplateProperty<bool>* prop) : TemplatePropertyTimeline<bool>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineBool(); }

    virtual std::string getClassName() const { return "PropertyTimelineBool"; }
};

class VRN_CORE_API PropertyTimelineIntBoundingbox : public TemplatePropertyTimeline<tgt::IntBounds> {
public:
    PropertyTimelineIntBoundingbox()  : TemplatePropertyTimeline<tgt::IntBounds>() {}
    PropertyTimelineIntBoundingbox(TemplateProperty<tgt::IntBounds>* prop) : TemplatePropertyTimeline<tgt::IntBounds>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineIntBoundingbox(); }

    virtual std::string getClassName() const { return "PropertyTimelineIntBoundingbox"; }
};

class VRN_CORE_API PropertyTimelineFloatBoundingbox : public TemplatePropertyTimeline<tgt::Bounds> {
public:
    PropertyTimelineFloatBoundingbox()  : TemplatePropertyTimeline<tgt::Bounds>() {}
    PropertyTimelineFloatBoundingbox(TemplateProperty<tgt::Bounds>* prop) : TemplatePropertyTimeline<tgt::Bounds>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineFloatBoundingbox(); }

    virtual std::string getClassName() const { return "PropertyTimelineFloatBoundingbox"; }
};



class VRN_CORE_API PropertyTimelineIVec2 : public TemplatePropertyTimeline<tgt::ivec2> {
public:
    PropertyTimelineIVec2()  : TemplatePropertyTimeline<tgt::ivec2>() {}
    PropertyTimelineIVec2(TemplateProperty<tgt::ivec2>* prop) : TemplatePropertyTimeline<tgt::ivec2>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineIVec2(); }

    virtual std::string getClassName() const { return "PropertyTimelineIVec2"; }
};

class VRN_CORE_API PropertyTimelineIVec3 : public TemplatePropertyTimeline<tgt::ivec3> {
public:
    PropertyTimelineIVec3()  : TemplatePropertyTimeline<tgt::ivec3>() {}
    PropertyTimelineIVec3(TemplateProperty<tgt::ivec3>* prop) : TemplatePropertyTimeline<tgt::ivec3>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineIVec3(); }

    virtual std::string getClassName() const { return "PropertyTimelineIVec3"; }
};

class VRN_CORE_API PropertyTimelineIVec4 : public TemplatePropertyTimeline<tgt::ivec4> {
public:
    PropertyTimelineIVec4()  : TemplatePropertyTimeline<tgt::ivec4>() {}
    PropertyTimelineIVec4(TemplateProperty<tgt::ivec4>* prop) : TemplatePropertyTimeline<tgt::ivec4>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineIVec4(); }

    virtual std::string getClassName() const { return "PropertyTimelineIVec4"; }
};

class VRN_CORE_API PropertyTimelineVec2 : public TemplatePropertyTimeline<tgt::vec2> {
public:
    PropertyTimelineVec2()  : TemplatePropertyTimeline<tgt::vec2>() {}
    PropertyTimelineVec2(TemplateProperty<tgt::vec2>* prop) : TemplatePropertyTimeline<tgt::vec2>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineVec2(); }

    virtual std::string getClassName() const { return "PropertyTimelineVec2"; }
};

class VRN_CORE_API PropertyTimelineVec3 : public TemplatePropertyTimeline<tgt::vec3> {
public:
    PropertyTimelineVec3()  : TemplatePropertyTimeline<tgt::vec3>() {}
    PropertyTimelineVec3(TemplateProperty<tgt::vec3>* prop) : TemplatePropertyTimeline<tgt::vec3>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineVec3(); }

    virtual std::string getClassName() const { return "PropertyTimelineVec3"; }
};

class VRN_CORE_API PropertyTimelineVec4 : public TemplatePropertyTimeline<tgt::vec4> {
public:
    PropertyTimelineVec4()  : TemplatePropertyTimeline<tgt::vec4>() {}
    PropertyTimelineVec4(TemplateProperty<tgt::vec4>* prop) : TemplatePropertyTimeline<tgt::vec4>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineVec4(); }

    virtual std::string getClassName() const { return "PropertyTimelineVec4"; }
};

class VRN_CORE_API PropertyTimelineMat2 : public TemplatePropertyTimeline<tgt::mat2> {
public:
    PropertyTimelineMat2()  : TemplatePropertyTimeline<tgt::mat2>() {}
    PropertyTimelineMat2(TemplateProperty<tgt::mat2>* prop) : TemplatePropertyTimeline<tgt::mat2>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineMat2(); }

    virtual std::string getClassName() const { return "PropertyTimelineMat2"; }
};

class VRN_CORE_API PropertyTimelineMat3 : public TemplatePropertyTimeline<tgt::mat3> {
public:
    PropertyTimelineMat3()  : TemplatePropertyTimeline<tgt::mat3>() {}
    PropertyTimelineMat3(TemplateProperty<tgt::mat3>* prop) : TemplatePropertyTimeline<tgt::mat3>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineMat3(); }

    virtual std::string getClassName() const { return "PropertyTimelineMat3"; }
};

class VRN_CORE_API PropertyTimelineMat4 : public TemplatePropertyTimeline<tgt::mat4> {
public:
    PropertyTimelineMat4()  : TemplatePropertyTimeline<tgt::mat4>() {}
    PropertyTimelineMat4(TemplateProperty<tgt::mat4>* prop) : TemplatePropertyTimeline<tgt::mat4>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineMat4(); }

    virtual std::string getClassName() const { return "PropertyTimelineMat4"; }
};

class VRN_CORE_API PropertyTimelineCamera : public TemplatePropertyTimeline<tgt::Camera> {
public:
    PropertyTimelineCamera()  : TemplatePropertyTimeline<tgt::Camera>() {}
    PropertyTimelineCamera(TemplateProperty<tgt::Camera>* prop) : TemplatePropertyTimeline<tgt::Camera>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineCamera(); }

    virtual std::string getClassName() const { return "PropertyTimelineCamera"; }
};

class VRN_CORE_API PropertyTimelineString : public TemplatePropertyTimeline<std::string> {
public:
    PropertyTimelineString()  : TemplatePropertyTimeline<std::string>() {}
    PropertyTimelineString(TemplateProperty<std::string>* prop) : TemplatePropertyTimeline<std::string>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineString(); }

    virtual std::string getClassName() const { return "PropertyTimelineString"; }
};

class VRN_CORE_API PropertyTimelineShaderSource : public TemplatePropertyTimeline<ShaderSource> {
public:
    PropertyTimelineShaderSource()  : TemplatePropertyTimeline<ShaderSource>() {}
    PropertyTimelineShaderSource(TemplateProperty<ShaderSource>* prop) : TemplatePropertyTimeline<ShaderSource>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineShaderSource(); }

    virtual std::string getClassName() const { return "PropertyTimelineShaderSource"; }
};

class VRN_CORE_API PropertyTimelineTransFunc1DKeys : public TemplatePropertyTimeline<TransFunc1DKeys*> {
public:
    PropertyTimelineTransFunc1DKeys()  : TemplatePropertyTimeline<TransFunc1DKeys*>() {}
    PropertyTimelineTransFunc1DKeys(TemplateProperty<TransFunc1DKeys*>* prop) : TemplatePropertyTimeline<TransFunc1DKeys*>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineTransFunc1DKeys(); }

    virtual std::string getClassName() const { return "PropertyTimelineTransFunc1DKeys"; }
};

class VRN_CORE_API PropertyTimelineTransFunc2DPrimitives : public TemplatePropertyTimeline<TransFunc2DPrimitives*> {
public:
    PropertyTimelineTransFunc2DPrimitives()  : TemplatePropertyTimeline<TransFunc2DPrimitives*>() {}
    PropertyTimelineTransFunc2DPrimitives(TemplateProperty<TransFunc2DPrimitives*>* prop) : TemplatePropertyTimeline<TransFunc2DPrimitives*>(prop) {}
    virtual AbstractSerializable* create() const { return new PropertyTimelineTransFunc2DPrimitives(); }

    virtual std::string getClassName() const { return "PropertyTimelineTransFunc"; }
};

} // namespace voreen

#endif
