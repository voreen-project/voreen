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

#ifndef VRN_SERIALIZABLEFACTORY_H
#define VRN_SERIALIZABLEFACTORY_H

#include "voreen/core/voreencoreapi.h"
#include <typeinfo>
#include <string>
#include <vector>

namespace voreen {

class VRN_CORE_API Serializable;

/**
 * The @c SerializableFactory is an interface for realizing classes which are supposed to
 * create @c Serializable objects according to a given type string.
 *
 * @par
 * Here is a short example of realizing the @c SerializableFactory interface:
 * @code
 * class Parent : public Serializable {
 * public:
 *     virtual void serialize(Serializer& s) const {}
 *     virtual void deserialize(Deserializer& s) {}
 * };
 *
 * class Child : public Parent {};
 *
 * class Factory : public SerializableFactory {
 * public:
 *     virtual const std::string getSerializableTypeString(const type_info& type) const {
 *         if (type == typeid(Parent))
 *             return "Parent";
 *         else if (type == typeid(Child))
 *             return "Child";
 *         else
 *             return "";
 *     }
 *
 *     virtual Serializable* createSerializableType(const std::string& typeString) {
 *         if (typeString == "Parent")
 *             return new Parent();
 *         else if (typeString == "Child")
 *             return new Child();
 *         else
 *             return "";
 *     }
 * };
 * @endcode
 * For more complex examples and interaction with other classes of the serialization framework
 * see the &quot;serializertest&quot; application in &quot;apps/&quot; directory.
 *
 * @note The @c SerializableFactory is necessary for implementation of (de-)serialization
 *       of polymorphic user defined data classes.
 *
 * @see Serializable
 * @see XmlSerializer
 * @see XmlDeserializer
 */
class VRN_CORE_API SerializableFactory {
public:
    /**
     * Virtual default destructor.
     */
    virtual ~SerializableFactory() {}

    /**
     * Returns a type string corresponding with the given @c type_info object,
     * which is created using @c typeid, in case this type is supported by the factory.
     *
     * @param type the @c type_info object
     *
     * @returns either the string corresponding with the given @c type_info object,
     *          or an empty string if the type is not supported by this factory.
     */
    virtual std::string getSerializableTypeString(const std::type_info& type) const = 0;

    /**
     * Returns a new instance of an serializable object corresponding with the given @c typeString.
     *
     * @param typeString the type string
     *
     * @returns either the new instance or @c 0 if the type string is not supported by this factory.
     */
    virtual Serializable* createSerializableType(const std::string& typeString) const = 0;
};

} // namespace

#endif // VRN_SERIALIZABLEFACTORY_H
