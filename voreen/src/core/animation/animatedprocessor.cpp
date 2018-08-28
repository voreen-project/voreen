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

#include "voreen/core/animation/animatedprocessor.h"
#include "voreen/core/animation/templatepropertytimeline.h"
#include "voreen/core/animation/serializationfactories.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/properties/propertyvector.h"

#include "tgt/camera.h"

namespace voreen {

const std::string AnimatedProcessor::loggerCat_("voreen.AnimatedProcessor");

AnimatedProcessor::AnimatedProcessor(Processor* processor)
    : processor_(processor)
{
}

PropertyTimeline* AnimatedProcessor::addTimeline(Property* prop) {
    PropertyTimeline* tl = PropertyTimelineFactory::getInstance()->createTimeline(prop);

    if (tl) {
        properties_.push_back(tl);
        return tl;
    }

    //PropertyVector* myprop14 = dynamic_cast<PropertyVector*>(prop);
    //if (myprop14 !=0) {
        //const std::vector<Property*>& tempproperties = myprop14->getProperties();
        //for (std::vector<Property*>::const_iterator it = tempproperties.begin(); it != tempproperties.end(); it++) {
            //addTemplatePropertyTimeline(*it);
        //}
        //return;
    //}

    return 0;
}

void AnimatedProcessor::removeTimeline(Property* prop) {
    std::vector<PropertyTimeline*>::iterator it;
    for (it = properties_.begin(); it!=properties_.end(); it++) {
        if((*it)->getProperty() == prop) {
            delete *it;
            properties_.erase(it);
            return;
        }
    }
}

AnimatedProcessor::AnimatedProcessor() {}

AnimatedProcessor::~AnimatedProcessor() {
    std::vector<PropertyTimeline*>::const_iterator it;
    for (it = properties_.begin(); it!=properties_.end(); ++it)
        delete (*it);

    properties_.clear();
}

bool AnimatedProcessor::isPropertyAnimated(const Property* p) {
    std::vector<PropertyTimeline*>::const_iterator it;
    for (it = properties_.begin(); it!=properties_.end(); ++it) {
        if((*it)->getProperty() == p)
            return true;
    }
    return false;
}

const std::vector<PropertyTimeline*>& AnimatedProcessor::getPropertyTimelines() const {
    return properties_;
}

void AnimatedProcessor::renderAt(float time) {
    //time = floor(time * 10000.f) / 10000.f;
    std::vector<PropertyTimeline*>::const_iterator it;
    for (it = properties_.begin(); it!=properties_.end(); ++it) {
        try {
            (*it)->renderAt(time);
        }
        catch (const std::exception& e) {
            LERROR("renderAt(): " << e.what() << " (" << getProcessorName() << ")");
        }
    }
}

const std::string AnimatedProcessor::getProcessorName() const {
    return processor_->getID();
}

const Processor* AnimatedProcessor::getCorrespondingProcessor() const {
    return processor_;
}

Processor* AnimatedProcessor::getCorrespondingProcessor() {
    return processor_;
}

void AnimatedProcessor::serialize(Serializer& s) const {
    s.serialize("processor", processor_);
    s.serialize("properties", properties_, "Property");
}

void AnimatedProcessor::deserialize(Deserializer& s) {
    s.deserialize("processor", processor_);
    s.deserialize("properties", properties_, "Property");

    // We need to remove timelines for properties that do not exist (any more)
    std::vector<PropertyTimeline*>::iterator it = properties_.begin();
    while (it != properties_.end()) {
        if ((*it)->getProperty() == 0) {
            delete *it;
            it = properties_.erase(it);
        } else {
            ++it;
        }
    }
}

} // namespace voreen
