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

#ifndef VRN_RESOURCEFACTORY_H
#define VRN_RESOURCEFACTORY_H

#include "voreen/core/io/serialization/serializablefactory.h"
#include "voreen/core/voreencoreapi.h"
#include "tgt/assert.h"
#include "tgt/logmanager.h"

#include <typeinfo>
#include <string>
#include <vector>
#include <map>

namespace voreen {

class VRN_CORE_API Serializable;

/**
 * Implementation of SerializeFactory that handles types based on their class name.
 * The types to be handled are populated to the factory via the registerTypes() methods.
 *
 * @note The resource type T must provide the following methods:
 *       - std::string getClassName(), returning a unique name for the class
 *       - Serializable* create(), returning a new instance of the class
 */
template<class T>
class ResourceFactory : public SerializableFactory {
public:
    /// Deletes the registered instances.
    virtual ~ResourceFactory();

    /**
     * Returns a type string corresponding to the given @c type_info object,
     * in case this type has been registered at this factory.
     *
     * @see SerializableFactory::getTypeString
     *
     * @returns either the string corresponding to the given @c type_info object,
     *          or an empty string, if the type has not been registered at this factory.
     */
    virtual std::string getSerializableTypeString(const std::type_info& type) const;

    /**
     * Returns a new instance of the class corresponding to the given @c typeString.
     *
     * @see SerializableFactory::createType
     *
     * @returns either the new instance or @c 0, if no class with this type string
     *  has been registered at this factory.
     */
    virtual Serializable* createSerializableType(const std::string& typeString) const;

    /**
     * Registers the passed type at the factory, using its class name as type name.
     *
     * @note The factory takes ownership of the passed object.
     */
    virtual void registerType(const T* type);

    /**
     * Registers the passed type for the passed typeName at the factory.
     *
     * @note The factory takes ownership of the passed object.
     */
    virtual void registerType(const std::string& typeName, const T* type);

    /**
     * Returns whether a type with the passed name has been registered at the factory.
     */
    bool isTypeKnown(const std::string& typeName) const;

    /**
     * Returns the names of all types that have been registered at the factory.
     */
    const std::vector<std::string>& getKnownTypeNames() const;

protected:
    /// List of the names of all types that have been registered at the factory.
    std::vector<std::string> typeNames_;

    /// Maps from class name to an instance of the corresponding type.
    std::map<std::string, const T*> typeMap_;
};

//---------------------------------------------------------------------------
// template definitions

template<class T>
ResourceFactory<T>::~ResourceFactory() {
    // delete all mapped type instances
    typename std::map<std::string, const T*>::const_iterator it;
    for (it = typeMap_.begin(); it != typeMap_.end(); ++it) {
        delete it->second;
    }

    typeMap_.clear();
    typeNames_.clear();
}

template<class T>
std::string ResourceFactory<T>::getSerializableTypeString(const std::type_info& type) const {
    typename std::map<std::string, const T*>::const_iterator it;
    for (it = typeMap_.begin(); it != typeMap_.end(); ++it) {
        if (type == typeid(*(it->second)))
            return it->first;
    }
    return "";
}

template<class T>
Serializable* ResourceFactory<T>::createSerializableType(const std::string& typeString) const {
    typename std::map<std::string, const T*>::const_iterator it = typeMap_.find(typeString);
    if (it == typeMap_.end())
        return 0;
    else
        return dynamic_cast<T*>(it->second->create());
}

template<class T>
void ResourceFactory<T>::registerType(const T* type) {
    tgtAssert(type, "null pointer passed");
    registerType(type->getClassName(), type);
}

template<class T>
void voreen::ResourceFactory<T>::registerType(const std::string& typeName, const T* type) {
    if (isTypeKnown(typeName)) {
        LWARNINGC("voreen.ResourceFactory", "type name '" << typeName << "' already registered");
    }
    else {
        typeMap_.insert(std::make_pair(typeName, type));
        typeNames_.push_back(typeName);
    }
}

template<class T>
bool voreen::ResourceFactory<T>::isTypeKnown(const std::string& typeName) const {
    return (typeMap_.find(typeName) != typeMap_.end());
}

template<class T>
const std::vector<std::string>& ResourceFactory<T>::getKnownTypeNames() const {
    return typeNames_;
}

} // namespace

#endif // VRN_RESOURCEFACTORY_H
