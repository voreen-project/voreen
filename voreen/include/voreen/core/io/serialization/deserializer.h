/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_DESERIALIZER_H
#define VRN_DESERIALIZER_H

#include <string>
#include <utility>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <iostream>
#include <typeinfo>

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"
#include "tgt/logmanager.h"

#include "voreen/core/io/serialization/xmlserializationconstants.h"
#include "voreen/core/io/serialization/serializationexceptions.h"
#include "voreen/core/io/serialization/serializable.h"

namespace voreen {

class XmlDeserializer;
class JsonDeserializer;

/**
 * This class can be understood as a superclass for all deserializers, without being
 * an actual superclass. All methods dynamically dispatch to the concrete deserializer.
 *
 * This is required because the (historic) deserializer interface makes heavy use of templates,
 * but there is no way in c++ to have virtual template methods.
 *
 * See JsonDeserializer or XmlDeserializer for the actual serializer implementation.
 */
class VRN_CORE_API Deserializer {
public:

    Deserializer(XmlDeserializer& deserializer);
    Deserializer(JsonDeserializer& deserializer);

    virtual ~Deserializer();

    /**
     * Returns the absolute working directory of the document, which is typically the path
     * to the file the document has been read from.
     */
    std::string getDocumentPath() const;

    /**
     * Sets whether to serialize pointer always as content instead of references.
     *
     * @attention This is not a cascading setting, which means that contained pointers
     *            are not serialized as content.
     *
     * @attention Serialization of all pointers as content can lead to redundant data.
     *
     * @param usePointerContentSerialization if @c true pointers are always serialized as content,
     *                                       otherwise as references when possible.
     */
    void setUsePointerContentSerialization(const bool& usePointerContentSerialization);

    /**
     * Returns whether pointers are always serialized as content instead of references.
     *
     * @return @c true if pointers are always serialized as content and @c false otherwise
     */
    bool getUsePointerContentSerialization() const;

    /**
     * Try to deserialize the given @c key/data pair.
     * If there is no such data the defaultValue is returned.
     */
    template<typename T>
    void optionalDeserialize(const std::string& key, T& data, const T& defaultValue);

    /**
     * Deserialize binary from a base64 encoded string to inputBuffer, memory has to be reserved in advance.
     */
    void deserializeBinaryBlob(const std::string& key, unsigned char* inputBuffer, size_t reservedMemory);

    /**
     * Deserialize binary @c data from a base64 encoded string.
     *
     * @note: The Objects HAVE to be primitive objects (i.e., "plain old data") because they will
     *        be deserialized using a byte wise operation! So no fancy objects with fancy constructors!
     */
    template<class T>
    void deserializeBinaryBlob(const std::string& key, std::vector<T>& data);

    /**
     * Deserialize binary @c data from a base64 encoded string.
     */
    void deserializeBinaryBlob(const std::string& key, std::vector<unsigned char>& buffer);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data variable to store deserialized data
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     */
    template<typename T>
    void deserialize(const std::string& key, T& data);

     /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data variable to store deserialized data
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     */
    template<typename T>
    void deserialize(const std::string& key, tgt::TemplateBounds<T>& data);

    /**
     * Deserializes the given pointer reference.
     *
     * @tparam type of referenced data
     *
     * @param key the Json node key
     * @param data variable to store deserialized pointer reference
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    void deserialize(const std::string& key, T*& data);

    /**
     * Deserializes a std::pair.
     *
     * @tparam S data type of first pair item
     * @tparam T data type of second pair item
     *
     * @param key the Json node key
     * @param data variable to store deserialized pair
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class S, class T>
    void deserialize(const std::string& key, std::pair<S, T>& data);

    /**
     * Deserializes a data vector.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of vector items
     *
     * @param key the Json node key
     * @param data variable to store deserialized data vector
     * @param itemKey Json node key for each Json child node
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    void deserialize(const std::string& key, std::vector<T>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Deserializes a data deque.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of vector items
     *
     * @param key the Json node key
     * @param data variable to store deserialized data vector
     * @param itemKey Json node key for each Json child node
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    void deserialize(const std::string& key, std::deque<T>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Deserializes a data list.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of vector items
     *
     * @param key the Json node key
     * @param data variable to store deserialized data vector
     * @param itemKey Json node key for each Json child node
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    void deserialize(const std::string& key, std::list<T>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Deserializes a data set.
     *
     * @note Element order of set items are not guaranteed to remains constant
     *       during serialization and deserialization due to limits of
     *       some STL containers like @c std::set.
     *
     * @tparam T data type of set items
     * @tparam C comparison class @see std::set
     *
     * @param key the Json node key
     * @param data variable to store deserialized data set
     * @param itemKey Json node key for each Json child node
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     * @throws XmlSerializationInvalidOperationException
     *         if pointer content serialization is enabled,
     *         because of possible hash problems on deserialization
     */
    template<class T, class C>
    void deserialize(const std::string& key, std::set<T, C>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Deserializes a data map.
     *
     * @note Element order of map items are not guaranteed to remains constant
     *       during serialization and deserialization due to limits of
     *       some STL containers like @c std::map.
     *
     * @tparam T data type of map keys
     * @tparam U data type of map values
     * @tparam C comparison class @see std::map
     *
     * @param key the Json node key
     * @param data variable to store deserialized data map
     * @param valueKey Json node key for each value node
     * @param keyKey key for each Json key node or attribute
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T, class U, class C>
    void deserialize(
        const std::string& key,
        std::map<T, U, C>& data,
        const std::string& valueKey = XmlSerializationConstants::VALUENODE,
        const std::string& keyKey = XmlSerializationConstants::KEYNODE);

    /**
     * Adds the given error @c message to the error list.
     *
     * @param message the error message
     */
    void addError(const std::string& message);

    /**
     * Adds the error message of the given @c exception to the error list.
     *
     * @param exception the exception
     */
    void addError(const std::exception& exception);

    /**
     * Adds the error message from the given @c exception to the error list
     * and raise the exception afterwards.
     *
     * @tparam T exception type
     *
     * @param exception the exception
     *
     * @throws SerializationException the exception is always thrown
     */
    template<class T>
    void raise(const T& exception);

    /**
     * Removes the last error message from the error list.
     */
    void removeLastError();

    /**
     * Returns the error list.
     *
     * @return the error list
     */
    const std::vector<std::string>& getErrors() const;

private:
    enum class DeserializerType {
        XML,
        JSON,
    };

    union DeserializerUnion {
        struct xmltype { DeserializerType type; XmlDeserializer* serializer; } xml;
        struct jsontype { DeserializerType type; JsonDeserializer* serializer; } json;

        DeserializerUnion(XmlDeserializer& s)
            //: xml(xmltype{ DeserializerType::XML, &s }) // Not supported by MSVC2012
        {
            xml.type = DeserializerType::XML;
            xml.serializer = &s;
        }

        DeserializerUnion(JsonDeserializer& s)
            //: json(jsontype{ DeserializerType::JSON, &s }) // Not supported by MSVC2012
        {
            json.type = DeserializerType::JSON;
            json.serializer = &s;
        }

        ~DeserializerUnion() {
            /*
            switch(xml.type) {
                case SerializerType::XML:
                    xml.~xmltype();
                    break;
                case SerializerType::JSON:
                    json.~jsontype();
                    break;
                // default omitted to emit compiler warning on addition of serializers
            }
            */
        }
    } concreteDeserializer;

    /**
     * Helper function for deserializing data collections like STL container.
     *
     * Serialization using pointer content does not make sense for constant iterator items,
     * that is why this function ignores the setting concerning usage of pointer content
     * serialization.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of collection
     *
     * @param key the XML node key
     * @param collection variable to store deserialized data collection
     * @param itemKey XML node key for each XML child node
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a XML node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    void deserializeCollectionConstIteratorItems(
        const std::string& key,
        T& collection,
        const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Helper function for deserializing data collections like STL container.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of collection
     *
     * @param key the XML node key
     * @param collection variable to store deserialized data collection
     * @param itemKey XML node key for each XML child node
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a XML node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    void deserializeCollection(
        const std::string& key,
        T& collection,
        const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Helper function for deserializing data maps like STL maps.
     *
     * @note Element order of map items are not guaranteed to remains constant
     *       during serialization and deserialization due to limits of
     *       some STL containers like @c std::map.
     *
     * @tparam T data type of map
     *
     * @param key the XML node key
     * @param map varaible to store deserialized data map
     * @param valueKey XML node key for each value node
     * @param keyKey key for each XML key node or attribute
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a XML node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws XmlSerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    void deserializeMap(
        const std::string& key,
        T& map,
        const std::string& valueKey = XmlSerializationConstants::VALUENODE,
        const std::string& keyKey = XmlSerializationConstants::KEYNODE);
};
} // namespace voreen

/// ----------------------------------------------------------------------------
/// Implementation -------------------------------------------------------------
/// ----------------------------------------------------------------------------

#include "voreen/core/io/serialization/jsondeserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

namespace voreen {

#define DISPATCH_TO_DESERIALIZER(fn) \
switch(concreteDeserializer.xml.type) {\
    case DeserializerType::XML:\
        concreteDeserializer.xml.serializer->fn ;\
        break;\
    case DeserializerType::JSON:\
        concreteDeserializer.json.serializer->fn ;\
        break;\
    /* default omitted to emit compiler warning on addition of serializers */\
}\

#define DISPATCH_TO_DESERIALIZER_RET(fn) \
switch(concreteDeserializer.xml.type) {\
    case DeserializerType::XML:\
        return concreteDeserializer.xml.serializer->fn ;\
    case DeserializerType::JSON:\
        return concreteDeserializer.json.serializer->fn ;\
    /* default omitted to emit compiler warning on addition of serializers */\
}\

template<typename T>
void Deserializer::optionalDeserialize(const std::string& key, T& data, const T& defaultValue) {
    try {
        deserialize(key, data);
    }
    catch (XmlSerializationNoSuchDataException e) {
        data = defaultValue;
        removeLastError();
    }
}

template<class T>
void Deserializer::deserializeBinaryBlob(const std::string& key, std::vector<T>& data) {
    size_t numItems;
    deserialize(key+".numItems", numItems);
    if (numItems > 0) {
        data.assign(numItems, T());
        deserializeBinaryBlob(key+".data", reinterpret_cast<unsigned char*>(&data[0]), sizeof(T) * data.size());
    }
}
template<typename T>
void Deserializer::deserialize(const std::string& key, T& data) {
    DISPATCH_TO_DESERIALIZER(deserialize(key, data));
}

template<class T>
void Deserializer::deserialize(const std::string& key, T*& data) {
    DISPATCH_TO_DESERIALIZER(deserialize(key, data));
}

template<class S, class T>
void Deserializer::deserialize(const std::string& key, std::pair<S, T>& data) {
    DISPATCH_TO_DESERIALIZER(deserialize(key, data));
}

template<class T>
void Deserializer::deserialize(const std::string& key, std::vector<T>& data, const std::string& itemKey) {
    deserializeCollection(key, data, itemKey);
}

template<class T>
void Deserializer::deserialize(const std::string& key, std::deque<T>& data, const std::string& itemKey) {
    deserializeCollection(key, data, itemKey);
}

template<class T>
void Deserializer::deserialize(const std::string& key, std::list<T>& data, const std::string& itemKey) {
    deserializeCollection(key, data, itemKey);
}

template<class T, class C>
void Deserializer::deserialize(const std::string& key, std::set<T, C>& data, const std::string& itemKey) {
    if(getUsePointerContentSerialization()) {
        raise(XmlSerializationInvalidOperationException("Set deserialization using pointer content is not permitted."));
    }
    deserializeCollectionConstIteratorItems(key, data, itemKey);
}

template<class T, class U, class C>
void Deserializer::deserialize(const std::string& key,
                                  std::map<T, U, C>& data,
                                  const std::string& valueKey,
                                  const std::string& keyKey) {
    deserializeMap(key, data, valueKey, keyKey);
}

template<class T>
void Deserializer::deserializeCollectionConstIteratorItems(const std::string& key, T& collection, const std::string& itemKey) {
    DISPATCH_TO_DESERIALIZER(deserializeCollectionConstIteratorItems(key, collection, itemKey));
}

template<class T>
void Deserializer::deserializeCollection(const std::string& key, T& collection, const std::string& itemKey) {
    DISPATCH_TO_DESERIALIZER(deserializeCollection(key, collection, itemKey));
}

template<class T>
void Deserializer::deserializeMap(const std::string& mapItemKey,
                                            T& map,
                                            const std::string& valueKey,
                                            const std::string& keyKey) {
    DISPATCH_TO_DESERIALIZER(deserializeMap(mapItemKey, map, valueKey, keyKey));
}

template<typename T>
void Deserializer::deserialize(const std::string& key, tgt::TemplateBounds<T>& data) {
    DISPATCH_TO_DESERIALIZER(deserialize(key, data));
}

template<class T>
void Deserializer::raise(const T& exception) {
    DISPATCH_TO_DESERIALIZER(raise(exception));
}


} // namespace voreen

#endif //VRN_DESERIALIZER_H
