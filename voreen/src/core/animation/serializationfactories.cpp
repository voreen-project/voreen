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

#include "voreen/core/animation/serializationfactories.h"

#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/animation/templatepropertytimeline.h"
#include "voreen/core/animation/propertykeyvalue.h"

#include "voreen/core/animation/interpolation/boolinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/camerainterpolationfunctions.h"
#include "voreen/core/animation/interpolation/floatinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/intinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/vec2interpolationfunctions.h"
#include "voreen/core/animation/interpolation/vec3interpolationfunctions.h"
#include "voreen/core/animation/interpolation/vec4interpolationfunctions.h"
#include "voreen/core/animation/interpolation/ivec2interpolationfunctions.h"
#include "voreen/core/animation/interpolation/ivec3interpolationfunctions.h"
#include "voreen/core/animation/interpolation/ivec4interpolationfunctions.h"
#include "voreen/core/animation/interpolation/mat2interpolationfunctions.h"
#include "voreen/core/animation/interpolation/mat3interpolationfunctions.h"
#include "voreen/core/animation/interpolation/mat4interpolationfunctions.h"
#include "voreen/core/animation/interpolation/transfunc1dkeysinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/transfunc2dprimitivesinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/shadersourceinterpolationfunctions.h"
#include "voreen/core/animation/interpolation/stringinterpolationfunctions.h"

namespace voreen {

PropertyTimelineFactory* PropertyTimelineFactory::instance_ = 0;

PropertyTimelineFactory::PropertyTimelineFactory() {
    registerType(new PropertyTimelineFloat);
    registerType(new PropertyTimelineInt);
    registerType(new PropertyTimelineBool);
    registerType(new PropertyTimelineIntBoundingbox);
    registerType(new PropertyTimelineFloatBoundingbox);
    registerType(new PropertyTimelineIVec2);
    registerType(new PropertyTimelineIVec3);
    registerType(new PropertyTimelineIVec4);
    registerType(new PropertyTimelineVec2);
    registerType(new PropertyTimelineVec3);
    registerType(new PropertyTimelineVec4);
    registerType(new PropertyTimelineMat2);
    registerType(new PropertyTimelineMat3);
    registerType(new PropertyTimelineMat4);
    registerType(new PropertyTimelineCamera);
    registerType(new PropertyTimelineString);
    registerType(new PropertyTimelineShaderSource);
    registerType(new PropertyTimelineTransFunc1DKeys);
    registerType(new PropertyTimelineTransFunc2DPrimitives);
}

PropertyTimeline* PropertyTimelineFactory::createTimeline(Property* p) const {
    std::map<std::string, const PropertyTimeline*>::const_iterator it;
    for (it = typeMap_.begin(); it != typeMap_.end(); ++it) {
        if(it->second->isCompatibleWith(p)) {
            AbstractSerializable* as = it->second->create();
            if(as) {
                PropertyTimeline* tl = dynamic_cast<PropertyTimeline*>(as);
                if(tl) {
                    tl->setProperty(p);
                    return tl;
                }
            }
        }
    }

    return 0;
}

bool PropertyTimelineFactory::canPropertyBeAnimated(const Property* p) const {
    std::map<std::string, const PropertyTimeline*>::const_iterator it;
    for (it = typeMap_.begin(); it != typeMap_.end(); ++it) {
        if(it->second->isCompatibleWith(p))
            return true;
    }

    return false;
}

PropertyTimelineFactory* PropertyTimelineFactory::getInstance() {
    if (!instance_)
        instance_ = new PropertyTimelineFactory();

    return instance_;
}

} // namespace voreen
