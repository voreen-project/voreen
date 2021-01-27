/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_NETWORKCONFIGURATOR_H
#define VRN_NETWORKCONFIGURATOR_H

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/volumeurlproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "tgt/vector.h"

#include <string>

namespace voreen {

class ProcessorNetwork;
class Processor;

/**
 * Configures the associated network via property value assignments,
 * which can be passed as strings.
 */
class VRN_CORE_API NetworkConfigurator {

public:
    NetworkConfigurator(ProcessorNetwork* network);

    /// Expects <ProcessorName>.<PropertyName>=<value>
    void setPropertyValue(const std::string& assignmentString);

    /// Expects the value to be set as a string.
    void setPropertyValue(const std::string& processorName, const std::string& propertyName, const std::string& valueStr);

    /// Expects <ProcessorName>.<PropertyName>
    void triggerButtonProperty(const std::string& qualifiedPropertyName);

    /// Expects the value to set as a string.
    void triggerButtonProperty(const std::string& processorName, const std::string& propertyName);

    void setPropertyValue(Processor* processor, BoolProperty* property, bool value);

    void setPropertyValue(Processor* processor, IntProperty* property, int value);

    void setPropertyValue(Processor* processor, FloatProperty* property, float value);

    void setPropertyValue(Processor* processor, IntVec2Property* property, const tgt::ivec2& value);

    void setPropertyValue(Processor* processor, IntVec3Property* property, const tgt::ivec3& value);

    void setPropertyValue(Processor* processor, IntVec4Property* property, const tgt::ivec4& value);

    void setPropertyValue(Processor* processor, FloatVec2Property* property, const tgt::vec2& value);

    void setPropertyValue(Processor* processor, FloatVec3Property* property, const tgt::vec3& value);

    void setPropertyValue(Processor* processor, FloatVec4Property* property, const tgt::vec4& value);

    void setPropertyValue(Processor* processor, IntIntervalProperty* property, const tgt::ivec2& value);

    void setPropertyValue(Processor* processor, FloatIntervalProperty* property, const tgt::vec2& value);

    void setPropertyValue(Processor* processor, StringProperty* property, const std::string& value);

    void setPropertyValue(Processor* processor, VolumeURLProperty* property, const std::string& value);

    void setPropertyValue(Processor* processor, OptionPropertyBase* property, const std::string& value);

private:
    void getProcessorAndProperty(const std::string& processorName, const std::string& propertyName,
        Processor*& processor, Property*& property);

    void parsePropertyName(const std::string& inputStr, std::string& processorName, std::string& propertyName);

    void parsePropertyAssignment(const std::string& inputStr, std::string& processorName, std::string& propertyName, std::string& value);

    ProcessorNetwork* network_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_NETWORKCONFIGURATOR_H
