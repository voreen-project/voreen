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

#ifndef VRN_VOREENSERIALIZABLEOBJECTFACTORY_H
#define VRN_VOREENSERIALIZABLEOBJECTFACTORY_H

#include "voreen/core/io/serialization/serializablefactory.h"
#include "voreen/core/voreenobject.h"

#include "tgt/logmanager.h"

#include <typeinfo>
#include <string>
#include <vector>
#include <map>

namespace voreen {

class Serializable;

/**
 * Implementation of SerializeFactory that handles VoreenObject classes based on their classname.
 * The types to be handled are populated to the factory via the registerTypes() methods.
 */
class VRN_CORE_API VoreenSerializableObjectFactory : public SerializableFactory {
public:
    /// Deletes the registered instances.
    virtual ~VoreenSerializableObjectFactory();

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
     * Returns a new instance of the serializable class corresponding to the given @c typeString.
     *
     * @see SerializableFactory::createType
     *
     * @returns either the new instance or @c 0, if no class with this type string
     *  has been registered at this factory.
     */
    virtual VoreenSerializableObject* createSerializableType(const std::string& typeString) const;

    /**
     * Returns the registered instance of the serializable class corresponding to the given @c typeString.
     *
     * @returns either the registered instance or @c 0, if no class with this type string
     *  has been registered at this factory.
     */
    virtual const VoreenSerializableObject* getSerializableType(const std::string& typeString) const;

    /**
     * Returns whether a serializable type with the passed name has been registered at the factory.
     */
    bool isSerializableTypeKnown(const std::string& typeName) const;

    /**
     * Returns all serializable types that have been registered at the factory.
     */
    const std::vector<const VoreenSerializableObject*> getSerializableTypes() const;

    /**
     * Returns all registered class instances of the specified type.
     */
    template<typename T>
    std::vector<const T*> getSerializableTypes() const;

    /**
     * Returns the names of all serializable types that have been registered at the factory.
     */
    const std::vector<std::string>& getSerializableTypeNames() const;

protected:
    /**
     * Registers the passed type at the factory, using its class name as type name.
     *
     * @note The factory takes ownership of the passed object.
     */
    virtual void registerSerializableType(const VoreenSerializableObject* type);

    /**
     * Registers the passed type at the factory.
     *
     * @note The factory takes ownership of the passed object.
     */
    virtual void registerSerializableType(const std::string& typeName, const VoreenSerializableObject* type);

    /// List of type instances that have been registered at the factory.
    std::vector<const VoreenSerializableObject*> types_;

    /// List of the names of all types that have been registered at the factory.
    std::vector<std::string> typeNames_;

    /// Maps from class name to an instance of the corresponding type.
    std::map<std::string, const VoreenSerializableObject*> typeMap_;

    static const std::string loggerCat_;
};

//-----------------------------------------------------------------------------
// template instantiations

template<typename T>
std::vector<const T*> VoreenSerializableObjectFactory::getSerializableTypes() const {
    std::vector<const T*> result;
    for (size_t i=0; i<types_.size(); i++) {
        if (dynamic_cast<const T*>(types_.at((i))))
            result.push_back(static_cast<const T*>(types_.at(i)));
    }
    return result;
}

} // namespace

#endif // VRN_SERIALIZABLEOBJECTFACTORY_H
