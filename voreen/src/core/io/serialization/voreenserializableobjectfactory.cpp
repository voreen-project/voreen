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

#include "voreen/core/io/serialization/voreenserializableobjectfactory.h"

#include "tgt/assert.h"

namespace voreen {

const std::string VoreenSerializableObjectFactory::loggerCat_("voreen.VoreenSerializableObjectFactory");

VoreenSerializableObjectFactory::~VoreenSerializableObjectFactory() {
    // delete all mapped type instances
    std::map<std::string, const VoreenSerializableObject*>::const_iterator it;
    for (it = typeMap_.begin(); it != typeMap_.end(); ++it) {
        delete it->second;
    }

    types_.clear();
    typeNames_.clear();
    typeMap_.clear();
}

std::string VoreenSerializableObjectFactory::getSerializableTypeString(const std::type_info& type) const {
    std::map<std::string, const VoreenSerializableObject*>::const_iterator it;
    for (it = typeMap_.begin(); it != typeMap_.end(); ++it) {
        if (type == typeid(*(it->second)))
            return it->first;
    }
    return "";
}

VoreenSerializableObject* VoreenSerializableObjectFactory::createSerializableType(const std::string& typeString) const {
    const VoreenSerializableObject* instance = getSerializableType(typeString);
    if (instance)
        return instance->create();
    else
        return 0;
}

const VoreenSerializableObject* VoreenSerializableObjectFactory::getSerializableType(const std::string& typeString) const {
    std::map<std::string, const VoreenSerializableObject*>::const_iterator it = typeMap_.find(typeString);
    if (it == typeMap_.end())
        return 0;
    else {
        return it->second;
    }
}

void VoreenSerializableObjectFactory::registerSerializableType(const VoreenSerializableObject* type)  {
    tgtAssert(type, "null pointer passed");
    tgtAssert(!type->getClassName().empty(), "passed type has empty classname");
    registerSerializableType(type->getClassName(), type);
}

void VoreenSerializableObjectFactory::registerSerializableType(const std::string& typeName, const VoreenSerializableObject* type) {
    tgtAssert(!typeName.empty(), "empty string passed as typename");
    tgtAssert(type, "null pointer passed");
    if (isSerializableTypeKnown(typeName)) {
        LERROR("type name '" << typeName << "' already registered");
    }
    else {
        typeNames_.push_back(typeName);
        types_.push_back(type);
        typeMap_.insert(std::make_pair(typeName, type));
    }
}

bool VoreenSerializableObjectFactory::isSerializableTypeKnown(const std::string& typeName) const {
    return (typeMap_.find(typeName) != typeMap_.end());
}

const std::vector<const VoreenSerializableObject*> VoreenSerializableObjectFactory::getSerializableTypes() const {
    return types_;
}

const std::vector<std::string>& VoreenSerializableObjectFactory::getSerializableTypeNames() const {
    return typeNames_;
}

} // namespace
