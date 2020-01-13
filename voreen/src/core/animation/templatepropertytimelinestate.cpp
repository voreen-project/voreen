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

#include "voreen/core/animation/templatepropertytimelinestate.h"
#include "voreen/core/animation/animation.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/animation/interpolation/camerainterpolationfunctions.h"
#include "voreen/core/animation/interpolation/transfunc1dkeysinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/transfunc2dprimitivesinterpolationfunctions.h"
#include "tgt/camera.h"

using tgt::Camera;

namespace voreen {

/////////////////// Special implementation for TransFuncBase*-Property

TransFunc1DKeysPropertyTimelineState::TransFunc1DKeysPropertyTimelineState(PropertyKeyValue<TransFunc1DKeys*>* kv) {
    values_.insert(std::pair<float,PropertyKeyValue<TransFunc1DKeys*>*>(kv->getTime(),kv));
}

TransFunc1DKeysPropertyTimelineState::TransFunc1DKeysPropertyTimelineState() {}


TransFunc1DKeysPropertyTimelineState::~TransFunc1DKeysPropertyTimelineState() {
    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::iterator it;
    for (it = values_.begin(); it != values_.end(); ++it) {
        if (it->second->getFollowingInterpolationFunction())
            delete (it->second->getFollowingInterpolationFunction());

        delete (it->second->getValue());
        delete (it->second);
    }
    values_.clear();
}

const PropertyKeyValue<TransFunc1DKeys*>* TransFunc1DKeysPropertyTimelineState::newKeyValue(float time) {
    time = round(time);
    if (values_.find(time) != values_.end())
        return 0;

    TransFunc1DKeys* value = const_cast<TransFunc1DKeys*>(getPropertyAt(time));

    PropertyKeyValue<TransFunc1DKeys*>* kv = new PropertyKeyValue<TransFunc1DKeys*>(value,time);

    values_.insert(std::pair<float,PropertyKeyValue<TransFunc1DKeys*>*>(time,kv));

    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::iterator it;
    it = values_.find(time);

    // if new value is the first value:
    if (it == values_.begin()) {
        // only do something if there are multiple values
        it++;
        if (it != values_.end()) {
            InterpolationFunction<TransFunc1DKeys*>* func = new TransFunc1DKeysInterpolationFunction();
            it->second->setForegoingInterpolationFunction(func);
            it--;
            it->second->setFollowingInterpolationFunction(func);
        }
    }
    else {
        it++;
        // if the new value is the last one
        if (it == values_.end()) {
            it--;
            InterpolationFunction<TransFunc1DKeys*>* func = new TransFunc1DKeysInterpolationFunction();
            it->second->setForegoingInterpolationFunction(func);
            it--;
            it->second->setFollowingInterpolationFunction(func);
            it++;
        }
        else {
            InterpolationFunction<TransFunc1DKeys*>* func1 = new TransFunc1DKeysInterpolationFunction();
            InterpolationFunction<TransFunc1DKeys*>* func2 = new TransFunc1DKeysInterpolationFunction();
            it->second->setForegoingInterpolationFunction(func2);
            it--;
            it->second->setFollowingInterpolationFunction(func2);
            it->second->setForegoingInterpolationFunction(func1);
            it--;
            delete (*it).second->getFollowingInterpolationFunction();
            it->second->setFollowingInterpolationFunction(func1);
            it++;
        }
    }

    return kv;
}

bool TransFunc1DKeysPropertyTimelineState::changeValueOfKeyValue(TransFunc1DKeys* value, const PropertyKeyValue<TransFunc1DKeys*>* keyvalue) {
    const float time = keyvalue->getTime();
    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::iterator it;
    it = values_.find(time);
    if (it != values_.end()) {
        if (value != it->second->getValue()) {
            delete (it->second->getValue());
            it->second->setValue(static_cast<TransFunc1DKeys*>(value->clone()));
        }
        return true;
    }
    return false;
}

DeleteKeyValueReturn TransFunc1DKeysPropertyTimelineState::deleteKeyValue(const PropertyKeyValue<TransFunc1DKeys*>* keyvalue) {
    const float time = keyvalue->getTime();
    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::iterator it;
    it = values_.find(time);

    // if wrong parameter do nothing
    if (it == values_.end())
        return KV_NOT_THERE;

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
        delete (it->second->getForegoingInterpolationFunction());
        it->second->setForegoingInterpolationFunction(0);
        it--;
        delete (it->second->getValue());
        values_.erase(time);
        return KV_DELETED;
    }
    // if value is the last one
    it++;
    if (it == values_.end()) {
        it--;
        it--;
        delete (it->second->getFollowingInterpolationFunction());
        it->second->setFollowingInterpolationFunction(0);
        it++;
        delete (it->second->getValue());
        values_.erase(time);
        return KV_DELETED;
    }
    // if value is in the middle
    InterpolationFunction<TransFunc1DKeys*>* func = new TransFunc1DKeysInterpolationFunction();
    it->second->setForegoingInterpolationFunction(func);
    it--;
    delete it->second->getFollowingInterpolationFunction();
    delete it->second->getValue();
    delete it->second->getForegoingInterpolationFunction();
    it--;
    it->second->setFollowingInterpolationFunction(func);
    values_.erase(time);
    return KV_DELETED;
}

TransFunc1DKeys* const TransFunc1DKeysPropertyTimelineState::getPropertyAt(float time) {
    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::iterator it;
    it = values_.find(time);
    if (it != values_.end()) {
        if (it->second->getValue()) {
            return static_cast<TransFunc1DKeys*>((*it).second->getValue()->clone());
        }
        else {
            LERRORC("TransFunc1DKeysPropertyTimelineState", "Keyvalue contains no value");
            return 0;
        }
    }

    it = values_.upper_bound(time);
    if (it == values_.begin())
        return static_cast<TransFunc1DKeys*>(it->second->getValue()->clone());

    if (it == values_.end()) {
        it--;
        return static_cast<TransFunc1DKeys*>(it->second->getValue()->clone());
    }
    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::iterator it2;
    it2 = it;
    it--;

    const InterpolationFunction<TransFunc1DKeys*>* func = (*it).second->getFollowingInterpolationFunction();
    const MultiPointInterpolationFunction<TransFunc1DKeys*>* multifunc = dynamic_cast<const MultiPointInterpolationFunction<TransFunc1DKeys*>*>(func);
    if (multifunc) {
        // call a function with multiple points
        // create vector of the used keyvalues
        std::vector<PropertyKeyValue<TransFunc1DKeys*>*> keys;
        // search for the first value in the multi-point interval
        while ((it!=values_.begin()) && (it->second->isSmoothed()))
            it--;

        do {
            keys.push_back((*it).second->clone());
            it++;
        } while ((it != values_.end()) && (it->second->isSmoothed()));

        if (it != values_.end())
            keys.push_back(it->second->clone());

        // interpolate value
        TransFunc1DKeys* returnvalue = multifunc->interpolate(keys,time);

        // delete all copied keys
        std::vector<PropertyKeyValue<TransFunc1DKeys*>*>::const_iterator delIt;
        for (delIt = keys.begin(); delIt != keys.end(); ++delIt)
            delete (*delIt);

        keys.clear();

        // return
        return returnvalue;
    }
    else {
        return func->interpolate(
            it->second->getValue(),
            it2->second->getValue(),
            (time - (it->first))/((it2->first) - (it->first)));
    }
}

TemplatePropertyTimelineState<TransFunc1DKeys*>* TransFunc1DKeysPropertyTimelineState::clone() const {
    TransFunc1DKeysPropertyTimelineState* timeline = new TransFunc1DKeysPropertyTimelineState();

    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::const_iterator it;
    for (it = values_.begin(); it != values_.end(); ++it) {
        timeline->values_.insert(std::pair<float,PropertyKeyValue<TransFunc1DKeys*>*>(it->first, it->second->clone()));
    }

    std::map<float,PropertyKeyValue<TransFunc1DKeys*>*>::const_iterator it2;
    it2 = timeline->values_.begin();
    for (it = values_.begin(); it != values_.end(); ) {
        const InterpolationFunction<TransFunc1DKeys*>* func;
        InterpolationFunction<TransFunc1DKeys*>* func2;
        func = it->second->getFollowingInterpolationFunction();
        if (!func) {
            it2++;
            it++;
            continue;
        }
        func2 = func->create();
        it2->second->setFollowingInterpolationFunction(func2);
        it2++;
        it++;
        it2->second->setForegoingInterpolationFunction(func2);
    }

    return timeline;
}

//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
TransFunc2DPrimitivesPropertyTimelineState::TransFunc2DPrimitivesPropertyTimelineState(PropertyKeyValue<TransFunc2DPrimitives*>* kv) {
    values_.insert(std::pair<float,PropertyKeyValue<TransFunc2DPrimitives*>*>(kv->getTime(),kv));
}

TransFunc2DPrimitivesPropertyTimelineState::TransFunc2DPrimitivesPropertyTimelineState() {}


TransFunc2DPrimitivesPropertyTimelineState::~TransFunc2DPrimitivesPropertyTimelineState() {
    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::iterator it;
    for (it = values_.begin(); it != values_.end(); ++it) {
        if (it->second->getFollowingInterpolationFunction())
            delete (it->second->getFollowingInterpolationFunction());

        delete (it->second->getValue());
        delete (it->second);
    }
    values_.clear();
}

const PropertyKeyValue<TransFunc2DPrimitives*>* TransFunc2DPrimitivesPropertyTimelineState::newKeyValue(float time) {
    time = round(time);
    if (values_.find(time) != values_.end())
        return 0;

    TransFunc2DPrimitives* value = const_cast<TransFunc2DPrimitives*>(getPropertyAt(time));

    PropertyKeyValue<TransFunc2DPrimitives*>* kv = new PropertyKeyValue<TransFunc2DPrimitives*>(value,time);

    values_.insert(std::pair<float,PropertyKeyValue<TransFunc2DPrimitives*>*>(time,kv));

    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::iterator it;
    it = values_.find(time);

    // if new value is the first value:
    if (it == values_.begin()) {
        // only do something if there are multiple values
        it++;
        if (it != values_.end()) {
            InterpolationFunction<TransFunc2DPrimitives*>* func = new TransFunc2DPrimitivesStartEndInterpolationFunction();
            it->second->setForegoingInterpolationFunction(func);
            it--;
            it->second->setFollowingInterpolationFunction(func);
        }
    }
    else {
        it++;
        // if the new value is the last one
        if (it == values_.end()) {
            it--;
            InterpolationFunction<TransFunc2DPrimitives*>* func = new TransFunc2DPrimitivesStartEndInterpolationFunction();
            it->second->setForegoingInterpolationFunction(func);
            it--;
            it->second->setFollowingInterpolationFunction(func);
            it++;
        }
        else {
            InterpolationFunction<TransFunc2DPrimitives*>* func1 = new TransFunc2DPrimitivesStartEndInterpolationFunction();
            InterpolationFunction<TransFunc2DPrimitives*>* func2 = new TransFunc2DPrimitivesStartEndInterpolationFunction();
            it->second->setForegoingInterpolationFunction(func2);
            it--;
            it->second->setFollowingInterpolationFunction(func2);
            it->second->setForegoingInterpolationFunction(func1);
            it--;
            delete (*it).second->getFollowingInterpolationFunction();
            it->second->setFollowingInterpolationFunction(func1);
            it++;
        }
    }

    return kv;
}

bool TransFunc2DPrimitivesPropertyTimelineState::changeValueOfKeyValue(TransFunc2DPrimitives* value, const PropertyKeyValue<TransFunc2DPrimitives*>* keyvalue) {
    const float time = keyvalue->getTime();
    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::iterator it;
    it = values_.find(time);
    if (it != values_.end()) {
        if (value != it->second->getValue()) {
            delete (it->second->getValue());
            it->second->setValue(static_cast<TransFunc2DPrimitives*>(value->clone()));
        }
        return true;
    }
    return false;
}

DeleteKeyValueReturn TransFunc2DPrimitivesPropertyTimelineState::deleteKeyValue(const PropertyKeyValue<TransFunc2DPrimitives*>* keyvalue) {
    const float time = keyvalue->getTime();
    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::iterator it;
    it = values_.find(time);

    // if wrong parameter do nothing
    if (it == values_.end())
        return KV_NOT_THERE;

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
        delete (it->second->getForegoingInterpolationFunction());
        it->second->setForegoingInterpolationFunction(0);
        it--;
        delete (it->second->getValue());
        values_.erase(time);
        return KV_DELETED;
    }
    // if value is the last one
    it++;
    if (it == values_.end()) {
        it--;
        it--;
        delete (it->second->getFollowingInterpolationFunction());
        it->second->setFollowingInterpolationFunction(0);
        it++;
        delete (it->second->getValue());
        values_.erase(time);
        return KV_DELETED;
    }
    // if value is in the middle
    InterpolationFunction<TransFunc2DPrimitives*>* func = new TransFunc2DPrimitivesStartEndInterpolationFunction();
    it->second->setForegoingInterpolationFunction(func);
    it--;
    delete it->second->getFollowingInterpolationFunction();
    delete it->second->getValue();
    delete it->second->getForegoingInterpolationFunction();
    it--;
    it->second->setFollowingInterpolationFunction(func);
    values_.erase(time);
    return KV_DELETED;
}

TransFunc2DPrimitives* const TransFunc2DPrimitivesPropertyTimelineState::getPropertyAt(float time) {
    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::iterator it;
    it = values_.find(time);
    if (it != values_.end()) {
        if (it->second->getValue()) {
            return static_cast<TransFunc2DPrimitives*>((*it).second->getValue()->clone());
        }
        else {
            LERRORC("TransFunc2DPrimitivesPropertyTimelineState", "Keyvalue contains no value");
            return 0;
        }
    }

    it = values_.upper_bound(time);
    if (it == values_.begin())
        return static_cast<TransFunc2DPrimitives*>(it->second->getValue()->clone());

    if (it == values_.end()) {
        it--;
        return static_cast<TransFunc2DPrimitives*>(it->second->getValue()->clone());
    }
    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::iterator it2;
    it2 = it;
    it--;

    const InterpolationFunction<TransFunc2DPrimitives*>* func = (*it).second->getFollowingInterpolationFunction();
    const MultiPointInterpolationFunction<TransFunc2DPrimitives*>* multifunc = dynamic_cast<const MultiPointInterpolationFunction<TransFunc2DPrimitives*>*>(func);
    if (multifunc) {
        // call a function with multiple points
        // create vector of the used keyvalues
        std::vector<PropertyKeyValue<TransFunc2DPrimitives*>*> keys;
        // search for the first value in the multi-point interval
        while ((it!=values_.begin()) && (it->second->isSmoothed()))
            it--;

        do {
            keys.push_back((*it).second->clone());
            it++;
        } while ((it != values_.end()) && (it->second->isSmoothed()));

        if (it != values_.end())
            keys.push_back(it->second->clone());

        // interpolate value
        TransFunc2DPrimitives* returnvalue = multifunc->interpolate(keys,time);

        // delete all copied keys
        std::vector<PropertyKeyValue<TransFunc2DPrimitives*>*>::const_iterator delIt;
        for (delIt = keys.begin(); delIt != keys.end(); ++delIt)
            delete (*delIt);

        keys.clear();

        // return
        return returnvalue;
    }
    else {
        return func->interpolate(
            it->second->getValue(),
            it2->second->getValue(),
            (time - (it->first))/((it2->first) - (it->first)));
    }
}

TemplatePropertyTimelineState<TransFunc2DPrimitives*>* TransFunc2DPrimitivesPropertyTimelineState::clone() const {
    TransFunc2DPrimitivesPropertyTimelineState* timeline = new TransFunc2DPrimitivesPropertyTimelineState();

    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::const_iterator it;
    for (it = values_.begin(); it != values_.end(); ++it) {
        timeline->values_.insert(std::pair<float,PropertyKeyValue<TransFunc2DPrimitives*>*>(it->first, it->second->clone()));
    }

    std::map<float,PropertyKeyValue<TransFunc2DPrimitives*>*>::const_iterator it2;
    it2 = timeline->values_.begin();
    for (it = values_.begin(); it != values_.end(); ) {
        const InterpolationFunction<TransFunc2DPrimitives*>* func;
        InterpolationFunction<TransFunc2DPrimitives*>* func2;
        func = it->second->getFollowingInterpolationFunction();
        if (!func) {
            it2++;
            it++;
            continue;
        }
        func2 = func->create();
        it2->second->setFollowingInterpolationFunction(func2);
        it2++;
        it++;
        it2->second->setForegoingInterpolationFunction(func2);
    }

    return timeline;
}


} // namespace voreen
