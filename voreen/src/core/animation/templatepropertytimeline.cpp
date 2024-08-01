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

#include "voreen/core/animation/templatepropertytimeline.h"
#include "voreen/core/animation/animation.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"
#include "voreen/core/properties/cameraproperty.h"
#include "tgt/camera.h"

using tgt::Camera;

namespace voreen {

////////////////// special case TransFunc1DKeys*
template <>
TemplatePropertyTimeline<TransFunc1DKeys*>::TemplatePropertyTimeline(TemplateProperty<TransFunc1DKeys*>* prop)
    : property_(prop)
    , activeOnRendering_(true)
    , timelineChanged_(false)
    , undoObserver_(0)
{
    duration_ = 60.f * 15.f;
    TransFunc1DKeys* func = property_->get();
    if (func)
        func = func->clone();
    else {
        LWARNINGC("", "warn 1");
        func = new TransFunc1DKeys();
        property_->set(func->clone());
    }
    timeline_ = new TransFunc1DKeysPropertyTimelineState(new PropertyKeyValue<TransFunc1DKeys*>(func,0));
}

template <>
TemplatePropertyTimeline<TransFunc1DKeys*>::TemplatePropertyTimeline()
    : property_(0)
    , activeOnRendering_(true)
    , timelineChanged_(false)
    , undoObserver_(0)
{
    duration_ = 60.f * 15.f;
    timeline_ = new TransFunc1DKeysPropertyTimelineState();
}

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::setProperty(Property* p) {
    TemplateProperty<TransFunc1DKeys*>* tp = dynamic_cast<TemplateProperty<TransFunc1DKeys*>*>(p);
    if(tp) {
        property_ = tp;
        if(timeline_->getKeyValues().size() == 0) {
            delete timeline_;
            timeline_ = 0;

            TransFunc1DKeys* func = property_->get();
            if (func)
                func = func->clone();
            else {
                LWARNINGC("", "warn 2");
                func = new TransFunc1DKeys();
                property_->set(func->clone());
            }
            timeline_ = new TransFunc1DKeysPropertyTimelineState(new PropertyKeyValue<TransFunc1DKeys*>(func,0));
        }
    }
    else {
        LERRORC("voreen.TemplatePropertyTimeline", "Property type mismatch!");
    }
}

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::resetTimeline() {
    timelineChanged_ = true;
    lastChanges_.push_back(timeline_);
    undoObserver_->animationChanged(this);

    TransFunc1DKeys* func = property_->get();
    if (func)
        func = func->clone();
    else {
        LWARNINGC("", "warn 3");
        func = new TransFunc1DKeys();
        property_->set(func->clone());
    }
    timeline_ = new TransFunc1DKeysPropertyTimelineState(new PropertyKeyValue<TransFunc1DKeys*>(func,0));

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it) {
        (*it)->timelineChanged();
    }
}

template <>
TransFunc1DKeys* TemplatePropertyTimeline<TransFunc1DKeys*>::privateGetPropertyAt(float time) {
    TemplatePropertyTimelineState<TransFunc1DKeys*>* tl = dynamic_cast<TemplatePropertyTimelineState<TransFunc1DKeys*>* >(timeline_);
    if (tl) {
        TransFunc1DKeys* func = const_cast<TransFunc1DKeys*>(tl->getPropertyAt(time));
        return func;
    }
    else {
        LERRORC("TemplatePropertyTimeline<TransFunc1DKeys*>", "no timeline state");
        return 0;
    }
}

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::renderAt(float time) {
    if (!activeOnRendering_)
        return;

    property_->set(getPropertyAt(time));
}

/*template <>
TemplatePropertyTimeline<TransFunc1DKeys*>::setCurrentSettingToKeyValue(PropertyKeyValueBase* key) {

    PropertyKeyValue<TransFunc1DKeys*>* tKey = dynamic_cast<PropertyKeyValue<TransFunc1DKeys*>*>(key);
    if (!tKey)
        return false;

    TransFunc1DKeys* value = property_->get();
    TransFunc1DKeys* animatedValue = tKey->getValue();

    timelineChanged_ = true;

    if (value != animatedValue)
        return changeValueOfKeyValue(value,tKey);
    else
        return true;
}*/

template <>
void TemplatePropertyTimeline<TransFunc1DKeys*>::setCurrentSettingAsKeyvalue(float time, bool forceKeyValue) {
    //time = floor(time*10000)/10000;
    TransFunc1DKeys* value = property_->get();
    TransFunc1DKeys* animatedValue = getPropertyAt(time);


    if (value && animatedValue) {
        if ((*value) != (*animatedValue)) {
            const PropertyKeyValue<TransFunc1DKeys*>* kv = timeline_->newKeyValue(time);
            if (!kv) {
                kv = new PropertyKeyValue<TransFunc1DKeys*>(value,time);
            }
            this->changeValueOfKeyValue(value,kv);
            return;
        }
        else {
            if (forceKeyValue){
                newKeyValue(time);
            }
        }
        return;
    }

    const PropertyKeyValue<TransFunc1DKeys*>* kv = timeline_->newKeyValue(time);
    if (!kv) {
        kv = new PropertyKeyValue<TransFunc1DKeys*>(value,time);
    }
    this->changeValueOfKeyValue(value,kv);
}

////////////////// special case TransFunc2DPrimitives*
template <>
TemplatePropertyTimeline<TransFunc2DPrimitives*>::TemplatePropertyTimeline(TemplateProperty<TransFunc2DPrimitives*>* prop)
    : property_(prop)
    , activeOnRendering_(true)
    , timelineChanged_(false)
    , undoObserver_(0)
{
    duration_ = 60.f * 15.f;
    TransFunc2DPrimitives* func = property_->get();
    if (func)
        func = func->clone();
    else {
        LWARNINGC("", "warn 1");
        func = new TransFunc2DPrimitives();
        property_->set(func->clone());
    }
    timeline_ = new TransFunc2DPrimitivesPropertyTimelineState(new PropertyKeyValue<TransFunc2DPrimitives*>(func,0));
}

template <>
TemplatePropertyTimeline<TransFunc2DPrimitives*>::TemplatePropertyTimeline()
    : property_(0)
    , activeOnRendering_(true)
    , timelineChanged_(false)
    , undoObserver_(0)
{
    duration_ = 60.f * 15.f;
    timeline_ = new TransFunc2DPrimitivesPropertyTimelineState();
}

template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::setProperty(Property* p) {
    TemplateProperty<TransFunc2DPrimitives*>* tp = dynamic_cast<TemplateProperty<TransFunc2DPrimitives*>*>(p);
    if(tp) {
        property_ = tp;
        if(timeline_->getKeyValues().size() == 0) {
            delete timeline_;
            timeline_ = 0;

            TransFunc2DPrimitives* func = property_->get();
            if (func)
                func = func->clone();
            else {
                LWARNINGC("", "warn 2");
                func = new TransFunc2DPrimitives();
                property_->set(func->clone());
            }
            timeline_ = new TransFunc2DPrimitivesPropertyTimelineState(new PropertyKeyValue<TransFunc2DPrimitives*>(func,0));
        }
    }
    else {
        LERRORC("voreen.TemplatePropertyTimeline", "Property type mismatch!");
    }
}

template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::resetTimeline() {
    timelineChanged_ = true;
    lastChanges_.push_back(timeline_);
    undoObserver_->animationChanged(this);

    TransFunc2DPrimitives* func = property_->get();
    if (func)
        func = func->clone();
    else {
        LWARNINGC("", "warn 3");
        func = new TransFunc2DPrimitives();
        property_->set(func->clone());
    }
    timeline_ = new TransFunc2DPrimitivesPropertyTimelineState(new PropertyKeyValue<TransFunc2DPrimitives*>(func,0));

    const std::vector<TimelineObserver*> timelineObservers = getObservers();
    std::vector<TimelineObserver*>::const_iterator it;
    for (it = timelineObservers.begin(); it != timelineObservers.end(); ++it) {
        (*it)->timelineChanged();
    }
}

template <>
TransFunc2DPrimitives* TemplatePropertyTimeline<TransFunc2DPrimitives*>::privateGetPropertyAt(float time) {
    TemplatePropertyTimelineState<TransFunc2DPrimitives*>* tl = dynamic_cast<TemplatePropertyTimelineState<TransFunc2DPrimitives*>* >(timeline_);
    if (tl) {
        TransFunc2DPrimitives* func = const_cast<TransFunc2DPrimitives*>(tl->getPropertyAt(time));
        return func;
    }
    else {
        LERRORC("TemplatePropertyTimeline<TransFunc2DPrimitives*>", "no timeline state");
        return 0;
    }
}

template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::renderAt(float time) {
    if (!activeOnRendering_)
        return;

    property_->set(getPropertyAt(time));
}

/*template <>
TemplatePropertyTimeline<TransFunc2DPrimitives*>::setCurrentSettingToKeyValue(PropertyKeyValueBase* key) {

    PropertyKeyValue<TransFunc2DPrimitives*>* tKey = dynamic_cast<PropertyKeyValue<TransFunc2DPrimitives*>*>(key);
    if (!tKey)
        return false;

    TransFunc2DPrimitives* value = property_->get();
    TransFunc2DPrimitives* animatedValue = tKey->getValue();

    timelineChanged_ = true;

    if (value != animatedValue)
        return changeValueOfKeyValue(value,tKey);
    else
        return true;
}*/

template <>
void TemplatePropertyTimeline<TransFunc2DPrimitives*>::setCurrentSettingAsKeyvalue(float time, bool forceKeyValue) {
    //time = floor(time*10000)/10000;
    TransFunc2DPrimitives* value = property_->get();
    TransFunc2DPrimitives* animatedValue = getPropertyAt(time);

    if (value && animatedValue) {
        if (*value != *animatedValue) {
            const PropertyKeyValue<TransFunc2DPrimitives*>* kv = timeline_->newKeyValue(time);
            if (!kv) {
                kv = new PropertyKeyValue<TransFunc2DPrimitives*>(value,time);
            }
            this->changeValueOfKeyValue(value,kv);
            return;
        }
        else {
            if (forceKeyValue){
                newKeyValue(time);
            }
        }
        return;
    }
    const PropertyKeyValue<TransFunc2DPrimitives*>* kv = timeline_->newKeyValue(time);
    if (!kv) {
        kv = new PropertyKeyValue<TransFunc2DPrimitives*>(value,time);
    }
    this->changeValueOfKeyValue(value,kv);
}

} // namespace voreen
