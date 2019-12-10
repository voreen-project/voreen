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

#include "voreen/core/io/serialization/xmlserializerbase.h"
#include <limits>

namespace voreen {

const std::string XmlSerializerBase::loggerCat_ = "voreen.XmlSerializerBase";

XmlSerializerBase::TemporaryNodeChanger::TemporaryNodeChanger(XmlSerializerBase& serializer, TiXmlNode* node)
    : serializer_(serializer)
    , storedNode_(serializer_.node_)
{
    serializer_.node_ = node;
}

XmlSerializerBase::TemporaryNodeChanger::~TemporaryNodeChanger() {
    serializer_.node_ = storedNode_;
}

XmlSerializerBase::XmlSerializerBase()
    : SerializerBase()
    , node_(0)
    , useAttributes_(false)
{
}


std::string XmlSerializerBase::getTypeString(const std::type_info& type) {
    for (FactoryListType::iterator it = factories_.begin(); it != factories_.end(); ++it) {
        std::string typeString = (*it)->getSerializableTypeString(type);
        if (!typeString.empty())
            return typeString;
    }

    return "";
}

bool XmlSerializerBase::isPrimitiveType(const std::type_info& type) const {
    if (type == typeid(bool)
        || type == typeid(char)
        || type == typeid(signed char)
        || type == typeid(unsigned char)
        || type == typeid(signed short)
        || type == typeid(unsigned short)
        || type == typeid(signed int)
        || type == typeid(unsigned int)
        || type == typeid(signed long)
        || type == typeid(unsigned long)
        || type == typeid(float)
        || type == typeid(double)
        || type == typeid(long double)
        || type == typeid(std::string)
        || type == typeid(uint8_t)
        || type == typeid(int8_t)
        || type == typeid(uint16_t)
        || type == typeid(int16_t)
        || type == typeid(uint32_t)
        || type == typeid(int32_t)
        || type == typeid(int64_t)
        || type == typeid(uint64_t))
    {
        return true;
    }

    return false;
}

bool XmlSerializerBase::isPrimitivePointerType(const std::type_info& type) const {
    if (type == typeid(bool*)
        || type == typeid(char*)
        || type == typeid(signed char*)
        || type == typeid(unsigned char*)
        || type == typeid(signed short*)
        || type == typeid(unsigned short*)
        || type == typeid(signed int*)
        || type == typeid(unsigned int*)
        || type == typeid(signed long*)
        || type == typeid(unsigned long*)
        || type == typeid(float*)
        || type == typeid(double*)
        || type == typeid(long double*)
        || type == typeid(std::string*)
        || type == typeid(uint8_t*)
        || type == typeid(int8_t*)
        || type == typeid(uint16_t*)
        || type == typeid(int16_t*)
        || type == typeid(uint32_t*)
        || type == typeid(int32_t*)
        || type == typeid(int64_t*)
        || type == typeid(uint64_t*))
    {
        return true;
    }

    return false;
}

void XmlSerializerBase::registerFactory(SerializableFactory* factory) {
    factories_.push_back(factory);
}

void XmlSerializerBase::registerFactories(const std::vector<SerializableFactory*>& factories) {
    for (size_t i=0; i<factories.size(); i++)
        registerFactory(factories.at(i));
}

void XmlSerializerBase::setUseAttributes(const bool& useAttributes) {
    useAttributes_ = useAttributes;
}

bool XmlSerializerBase::getUseAttributes() const {
    return useAttributes_;
}

std::string XmlSerializerBase::convertDataToString(const float& data) {
    std::stringstream stream;
    stream.precision(std::numeric_limits<float>::max_digits10);
    stream.setf(std::ios::fixed);

    stream << data;
    std::string s = stream.str();

    // remove trailing zeros and decimal point
    size_t point_pos = s.find('.');
    if (point_pos != std::string::npos) {
        size_t new_size = s.size();
        for (size_t i = s.size() - 1; i >= point_pos; i--) {
            if (s[i] == '0' || s[i] == '.')
                new_size--;
            else
                break;
        }
        s.resize(new_size);
    }

    return s;
}

std::string XmlSerializerBase::convertDataToString(const double& data) {
    std::stringstream stream;
    stream.precision(std::numeric_limits<double>::max_digits10);
    stream.setf(std::ios::fixed);

    stream << data;
    std::string s = stream.str();

    // remove trailing zeros and decimal point
    size_t point_pos = s.find('.');
    if (point_pos != std::string::npos) {
        size_t new_size = s.size();
        for (size_t i = s.size() - 1; i >= point_pos; i--) {
            if (s[i] == '0' || s[i] == '.')
                new_size--;
            else
                break;
        }
        s.resize(new_size);
    }

    return s;
}



} // namespace
