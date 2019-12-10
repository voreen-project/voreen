/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/network/networkconfigurator.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/utils/exception.h"
#include "voreen/core/utils/stringutils.h"

#include "voreen/core/properties/buttonproperty.h"

#include "tgt/logmanager.h"
#include "tgt/exception.h"
#include "tgt/assert.h"

namespace voreen {

const std::string NetworkConfigurator::loggerCat_("voreentool.NetworkConfigurator");

NetworkConfigurator::NetworkConfigurator(ProcessorNetwork* network) :
    network_(network)
{
    tgtAssert(network_, "null pointer passed");
}

void NetworkConfigurator::setPropertyValue(const std::string& assignmentString) {
    std::string procName, propName, valueStr;
    parsePropertyAssignment(assignmentString, procName, propName, valueStr);
    setPropertyValue(procName, propName, valueStr);
}

void NetworkConfigurator::setPropertyValue(const std::string& processorName, const std::string& propertyName, const std::string& valueStr) {
    Processor* processor = 0;
    Property* property = 0;
    getProcessorAndProperty(processorName, propertyName, processor, property);
    tgtAssert(processor && property, "null pointer returned (exception expected)");

    if (dynamic_cast<BoolProperty*>(property)) {
        bool value;
        if (toLower(valueStr) == "false" || valueStr == "0")
            value = false;
        else if (toLower(valueStr) == "true" || valueStr == "1")
            value = true;
        else
            throw VoreenException("failed to convert '" + valueStr + "' to bool");
       setPropertyValue(processor, static_cast<BoolProperty*>(property), value);
    }
    else if (dynamic_cast<IntProperty*>(property)) {
        int value;
        try {
            value = genericFromString<int>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to int");
        }
        setPropertyValue(processor, static_cast<IntProperty*>(property), value);
    }
    else if (dynamic_cast<FloatProperty*>(property)) {
        float value;
        try {
            value = genericFromString<float>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to float");
        }
        setPropertyValue(processor, static_cast<FloatProperty*>(property), value);
    }
    else if (dynamic_cast<IntVec2Property*>(property)) {
        tgt::ivec2 value;
        try {
            value = genericFromString<tgt::ivec2>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to ivec2 (expected format: [x y])");
        }
        setPropertyValue(processor, static_cast<IntVec2Property*>(property), value);
    }
    else if (dynamic_cast<IntVec3Property*>(property)) {
        tgt::ivec3 value;
        try {
            value = genericFromString<tgt::ivec3>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to ivec3 (expected format: [x y z])");
        }
        setPropertyValue(processor, static_cast<IntVec3Property*>(property), value);
    }
    else if (dynamic_cast<IntVec4Property*>(property)) {
        tgt::ivec4 value;
        try {
            value = genericFromString<tgt::ivec4>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to ivec4 (expected format: [x y z w])");
        }
        setPropertyValue(processor, static_cast<IntVec4Property*>(property), value);
    }
    else if (dynamic_cast<FloatVec2Property*>(property)) {
        tgt::vec2 value;
        try {
            value = genericFromString<tgt::vec2>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to vec2 (expected format: [x y])");
        }
        setPropertyValue(processor, static_cast<FloatVec2Property*>(property), value);
    }
    else if (dynamic_cast<FloatVec3Property*>(property)) {
        tgt::vec3 value;
        try {
            value = genericFromString<tgt::vec3>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to vec3 (expected format: [x y z])");
        }
        setPropertyValue(processor, static_cast<FloatVec3Property*>(property), value);
    }
    else if (dynamic_cast<FloatVec4Property*>(property)) {
        tgt::vec4 value;
        try {
            value = genericFromString<tgt::vec4>(valueStr);
        }
        catch (VoreenException& /*e*/) {
            throw VoreenException("failed to convert '" + valueStr + "' to vec4 (expected format: [x y z w])");
        }
        setPropertyValue(processor, static_cast<FloatVec4Property*>(property), value);
    }
    else if (dynamic_cast<VolumeURLProperty*>(property)) {
        setPropertyValue(processor, static_cast<VolumeURLProperty*>(property), valueStr);
    }
    else if (dynamic_cast<OptionPropertyBase*>(property)) {
        setPropertyValue(processor, static_cast<OptionPropertyBase*>(property), valueStr);
    }
    else if (dynamic_cast<StringProperty*>(property)) {
        setPropertyValue(processor, static_cast<StringProperty*>(property), valueStr);
    }
    else if (dynamic_cast<ButtonProperty*>(property)) {
        throw VoreenException("'" + property->getFullyQualifiedID() + "' is a ButtonProperty and therefore cannot be assigned a value");
    }
    else {
        throw VoreenException(property->getFullyQualifiedID() + " is of unsupported type: " + property->getClassName());
    }

}

void NetworkConfigurator::triggerButtonProperty(const std::string& qualifiedPropertyName) {
    std::string procName, propName;
    parsePropertyName(qualifiedPropertyName, procName, propName);
    triggerButtonProperty(procName, propName);
}

void NetworkConfigurator::triggerButtonProperty(const std::string& processorName, const std::string& propertyName) {
    Processor* processor = 0;
    Property* property = 0;
    getProcessorAndProperty(processorName, propertyName, processor, property);
    tgtAssert(processor && property, "null pointer returned (exception expected)");

    ButtonProperty* buttonProperty = dynamic_cast<ButtonProperty*>(property);
    if (buttonProperty)
        buttonProperty->clicked();
    else
        throw VoreenException("'" + property->getFullyQualifiedID() + "' is not a ButtonProperty");
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, BoolProperty* property, bool value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, IntProperty* property, int value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, FloatProperty* property, float value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, IntVec2Property* property, const tgt::ivec2& value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, IntVec3Property* property, const tgt::ivec3& value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, IntVec4Property* property, const tgt::ivec4& value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, FloatVec2Property* property, const tgt::vec2& value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, FloatVec3Property* property, const tgt::vec3& value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, FloatVec4Property* property, const tgt::vec4& value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, StringProperty* property, const std::string& value) {
    tgtAssert(property, "null pointer passed");
    property->set(value);
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, VolumeURLProperty* property, const std::string& value) {
    tgtAssert(property, "null pointer passed");

    property->setURL(value);
    try {
        LINFO(property->getFullyQualifiedID() << ": loading volume " << value << "");
        property->loadVolume();
    }
    catch (std::exception& e) {
        throw VoreenException(property->getFullyQualifiedID() + ": " + e.what());
    }
}

void NetworkConfigurator::setPropertyValue(Processor* /*processor*/, OptionPropertyBase* property, const std::string& value) {
    tgtAssert(property, "null pointer passed");
    property->select(value);
}

// private functions

void NetworkConfigurator::getProcessorAndProperty(const std::string& processorName, const std::string& propertyName,
    Processor*& processor, Property*& property) {
    tgtAssert(network_, "no processor network");

    processor = network_->getProcessorByName(processorName);
    if (!processor)
        throw VoreenException("network contains no processor with name '" + processorName + "'");

    property = processor->getProperty(propertyName);
    if (!property)
        property = processor->getPropertyByName(propertyName);
    if (!property)
        throw VoreenException("processor '" + processorName + "' has no property with id or name '" + propertyName + "'");
}

void NetworkConfigurator::parsePropertyName(const std::string& inputStr, std::string& processorName, std::string& propertyName) {
    const std::string expectedFormat = "ProcessorName.PropertyName";

    std::vector<std::string> tokens = strSplit(inputStr, '.');
    if (tokens.size() <= 1)
        throw VoreenException("invalid property name '" + inputStr + "' (expected format: " + expectedFormat + ")");

    processorName = tokens.at(0);
    tokens.erase(tokens.begin());
    propertyName = strJoin(tokens, ".");

    if (processorName.empty() || propertyName.empty())
        throw VoreenException("invalid property name '" + inputStr + "' (expected format: " + expectedFormat + ")");
}

void NetworkConfigurator::parsePropertyAssignment(const std::string& inputStr, std::string& processorName, std::string& propertyName, std::string& value) {
    const std::string expectedFormat = "ProcessorName.PropertyName=value";

    std::vector<std::string> tokens = strSplit(inputStr, '=');
    if (tokens.size() != 2)
        throw VoreenException("invalid property value assignment (expected format: " + expectedFormat + ")");

    try {
        parsePropertyName(tokens.at(0), processorName, propertyName);
    }
    catch (const VoreenException& /*e*/) {
        throw VoreenException("invalid property value assignment (expected format: " + expectedFormat + ")");
    }

    value = tokens.at(1);
    if (value.empty())
        throw VoreenException("invalid property value assignment (expected format: " + expectedFormat + ")");
}

} // namespace
