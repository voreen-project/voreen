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

#ifndef VRN_XMLSERIALIZERBASE_H
#define VRN_XMLSERIALIZERBASE_H

#include <string>
#include <vector>

#include "tinyxml/tinyxml.h"

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/serialization/serializationexceptions.h"
#include "voreen/core/io/serialization/xmlserializationconstants.h"
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/serializerbase.h"
#include "voreen/core/io/serialization/serializablefactory.h"

namespace voreen {

/**
 * The @c XmlSerializerBase class encapsulates functions that are common to @c XmlSerializer
 * and @c XmlDeserializer.
 *
 * @see XmlSerializer
 * @see XmlDeserializer
 */
class VRN_CORE_API XmlSerializerBase : public SerializerBase {
public:
    /**
     * Default constructor.
     */
    XmlSerializerBase();

    /**
     * Returns whether primitive data is serialized as XML attributes or XML nodes.
     *
     * @return @c true if primitive data is serialized as XML attributes and @c false otherwise
     */
    bool getUseAttributes() const;

    /**
     * Sets whether to serialize primitive data as XML attributes or XML nodes.
     *
     * @attention All pointer references to primitive data which are serialized
     *            using XML attributes cannot be resolved. Furthermore using the
     *            same key for different values leads to a
     *            @c SerializationAttributeNamingException.
     *
     * @param useAttributes if @c true serialize primitive data as XML attributes,
     *                      otherwise as XML nodes
     */
    void setUseAttributes(const bool& useAttributes);

    /**
     * Registers an @c SerializableFactory to support serialization of polymorphic classes.
     *
     * @attention Ensure that the factory exists during the whole serialization process.
     *
     * @param factory the factory
     */
    void registerFactory(SerializableFactory* factory);

    /**
     * Convenience method for serialization factory registration.
     *
     * @see registerFactory
     */
    void registerFactories(const std::vector<SerializableFactory*>& factories);

protected:
    /**
     * Category for logging.
     */
    static const std::string loggerCat_;

protected:
    /**
     * Converts the given data to @c std::string.
     *
     * @tparam T type of data to convert
     *
     * @param data data to convert
     *
     * @return the string
     */
    template<class T>
    inline std::string convertDataToString(const T& data);

    /**
     * Converts the given data to @c std::string.
     *
     * @tparam T type of data to convert
     *
     * @param data data to convert
     *
     * @return the string
     */
    std::string convertDataToString(const float& data);

    /**
     * Converts the given data to @c std::string.
     *
     * @tparam T type of data to convert
     *
     * @param data data to convert
     *
     * @return the string
     */
    std::string convertDataToString(const double& data);

    /**
     * Converts the given @c std::string to data.
     *
     * @tparam T type of data to convert
     *
     * @param str the string
     * @param data variable to store converted data
     */
    template<class T>
    inline void convertStringToData(const std::string& str, T& data);

    /**
     * Returns type string corresponding with given @c type,
     * when it is supported by a registered factory.
     *
     * @param type the type
     *
     * @returns either the type string or an empty string
     *     when the type is not supported by any registered factory.
     */
    std::string getTypeString(const std::type_info& type);

    /**
     * Returns if given type is a primitive data type.
     *
     * @return @c true if it is a primitive data type and @c false otherwise
     */
    bool isPrimitiveType(const std::type_info& type) const;

    /**
     * Returns if given type is a primitive data pointer type.
     *
     * @return @c true if it is a primitive data pointer type and @c false otherwise.
     */
    bool isPrimitivePointerType(const std::type_info& type) const;

    /**
     * This is a helper class to ensure correct state of the XML node
     * for inserting or reading data.
     *
     * @note As C++ does not support a finally block statement, we need this
     *       class to ensure that cleanup code concerning the XML node
     *       for inserting or reading data is executed.
     */
    class VRN_CORE_API TemporaryNodeChanger {
    public:
        /**
         * Creates a @c TemporaryNodeChange, which changes the actual XML node for
         * inserting or reading data.
         *
         * @param serializer serializer whose XML node should be changed
         * @param node the new node
         */
        TemporaryNodeChanger(XmlSerializerBase& serializer, TiXmlNode* node);

        /**
         * Destructor ensures restoring the XML node for inserting or reading data
         * which was set before this instance was created.
         */
        ~TemporaryNodeChanger();

    private:
        /**
         * Serializer whose XML node is changed.
         */
        XmlSerializerBase& serializer_;

        /**
         * XML node which was set before this @c TemporaryNodeChanger was created.
         */
        TiXmlNode* storedNode_;
    };

    /**
     * XML document that contains the already serialized or deserialized data.
     *
     * @attention Keep in mind that may not all pointer references are already resolved.
     */
    TiXmlDocument document_;

    /**
     * XML node for inserting or reading data.
     */
    friend class TemporaryNodeChanger;
    TiXmlNode* node_;

    /**
     * If @c true all primitive data is serialized as XML attributes, otherwise as XML nodes.
     */
    bool useAttributes_;

    /**
     * Type definition for an @c SerializableFactory list.
     *
     * @note The @c XmlSerializerBase does not own the @c SerializableFactory objects, so remember
     *       to delete the @c SerializableFactory objects where you have created them.
     */
    typedef std::vector<SerializableFactory*> FactoryListType;

    /**
     * List of registered @c SerializableFactory objects.
     *
     * @note @c SerializableFactory objects are necessary
     *       to support polymorphic class serialization.
     */
    FactoryListType factories_;

};

template<class T>
inline std::string XmlSerializerBase::convertDataToString(const T& data) {
    std::stringstream stream;
    stream << data;
    return stream.str();
}

template<class T>
inline void XmlSerializerBase::convertStringToData(const std::string& str, T& data) {
    std::stringstream stream;
    stream << str;
    stream >> data;
}

} // namespace

#endif // VRN_XMLSERIALIZER_H
