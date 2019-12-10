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

#ifndef VRN_SERIALIZER_H
#define VRN_SERIALIZER_H

#include <string>
#include <utility>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <iostream>
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmlserializationconstants.h"
#include "voreen/core/io/serialization/serializationexceptions.h"

namespace voreen {

class XmlSerializer;
class JsonSerializer;

/**
 * This class can be understood as a superclass for all serializers, without being
 * an actual superclass. All methods dynamically dispatch to the concrete serializer.
 *
 * This is required because the (historic) serializer interface makes heavy use of templates,
 * but there is no way in c++ to have virtual template methods.
 *
 * See JsonSerializer or XmlSerializer for the actual serializer implementation.
 */
class VRN_CORE_API Serializer
{
public:
    Serializer(JsonSerializer& jsonserializer);
    Serializer(XmlSerializer& jsonserializer);

    virtual ~Serializer();

    /**
     * Returns the absolute working directory of the document, which is typically the path
     * to the file the document will be written to.
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
     * Serialize the given @c key/data pair if data != defaultValue.
     */
    template<typename T>
    void optionalSerialize(const std::string& key, const T& data, const T& defaultValue);

    /**
     * Serialize the given binary @c data as a base64 encoded string.
     */
    void serializeBinaryBlob(const std::string& key, const unsigned char* data, size_t length);

    /**
     * Binary serialize a std::vector of objects.
     *
     * @note: The Objects HAVE to be primitive objects (i.e., "plain old data") because they will
     *        be deserialized using a byte wise operation! So no fancy objects with fancy constructors!
     */
    template <class T>
    void serializeBinaryBlob(const std::string& key, const std::vector<T>& data);

    /**
     * Serialize the given binary @c data vector as a base64 encoded string.
     */
    void serializeBinaryBlob(const std::string& key, const std::vector<unsigned char>& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data, of a simple/tgt/Serializable type
     *
     */
    template<typename T>
    void serialize(const std::string& key, const T& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    template<typename T>
    void serialize(const std::string& key, const tgt::TemplateBounds<T>& data);

    /**
     * Serializes the given pointer reference.
     *
     * @tparam type of referenced data
     *
     * @param key the Json node key
     * @param data the pointer reference
     */
    template<class T>
    void serialize(const std::string& key, const T* const& data);

    /**
     * Serializes the given std::pair.
     *
     * @tparam S data type of first pair item
     * @tparam T data type of second pair item
     *
     * @param key the Json node key
     * @param data the pair to serialize
     *
     */
    template<class S, class T>
    void serialize(const std::string& key, const std::pair<S, T>& data);

    /**
     * Serializes the given data vector.
     *
     * @par
     * Iterates over the given data collection and serialize each collection item.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of vector items
     *
     * @param key the Json node key
     * @param data the data vector
     * @param itemKey Json node key for each Json child node
     *
     */
    template<class T>
    void serialize(const std::string& key, const std::vector<T>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Serializes the given data deque.
     *
     * @par
     * Iterates over the given data collection and serialize each collection item.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of vector items
     *
     * @param key the Json node key
     * @param data the data vector
     * @param itemKey Json node key for each Json child node
     *
     */
    template<class T>
    void serialize(const std::string& key, const std::deque<T>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Serializes the given data list.
     *
     * @par
     * Iterates over the given data collection and serialize each collection item.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of vector items
     *
     * @param key the Json node key
     * @param data the data vector
     * @param itemKey Json node key for each Json child node
     *
     */
    template<class T>
    void serialize(const std::string& key, const std::list<T>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Serializes the given data set.
     *
     * @par
     * Iterates over the given data collection and serialize each collection item.
     *
     * @note Element order of set items are not guaranteed to remains constant
     *       during serialization and deserialization due to limits of
     *       some STL containers like @c std::set.
     *
     * @tparam T data type of set items
     * @tparam C comparison class @see std::set
     *
     * @param key the Json node key
     * @param data the data set
     * @param itemKey Json node key for each Json child node
     *
     */
    template<class T, class C>
    void serialize(const std::string& key, const std::set<T, C>& data, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Serializes the given data map.
     *
     * @par
     * Iterates over the given map and serialize each @c key/value pair item of the map.
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
     * @param data the data map
     * @param valueKey Json node key for each value node
     * @param keyKey key for each Json key node or attribute
     *
     */
    template<class T, class U, class C>
    void serialize(
        const std::string& key,
        const std::map<T, U, C>& data,
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
    enum class SerializerType {
        XML,
        JSON,
    };

    union SerializerUnion {
        struct xmltype { SerializerType type; XmlSerializer* serializer; } xml;
        struct jsontype { SerializerType type; JsonSerializer* serializer; } json;

        SerializerUnion(XmlSerializer& s)
            : xml({ SerializerType::XML, &s })
        {
        }

        SerializerUnion(JsonSerializer& s)
            : json({ SerializerType::JSON, &s })
        {
        }

        ~SerializerUnion() {
            /*
            Not required, as we take references to serializers
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
    } concreteSerializer;

    /**
     * A list containing all error messages that were raised during serialization process.
     */
    std::vector<std::string> errors_;

    /**
     * Helper function for serializing data collections like STL container.
     *
     * @note Element order of collection items remains constant during
     *       serialization and deserialization.
     *
     * @tparam T data type of collection
     *
     * @param key the Json node key
     * @param collection the data collection
     * @param itemKey Json node key for each Json child node
     *
     */
    template<class T>
    void serializeCollection(const std::string& key, const T& collection, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Helper function for serializing data maps like STL maps.
     *
     * @note Element order of map items are not guaranteed to remains constant
     *       during serialization and deserialization due to limits of
     *       some STL containers like @c std::map.
     *
     * @tparam T data type of map
     *
     * @param key the Json node key
     * @param map the map
     * @param valueKey Json node key for each value node
     * @param keyKey key for each Json key node or attribute
     *
     */
    template<class T>
    void serializeMap(
        const std::string& key,
        const T& map,
        const std::string& valueKey = XmlSerializationConstants::VALUENODE,
        const std::string& keyKey = XmlSerializationConstants::KEYNODE);
};


/// ----------------------------------------------------------------------------
/// Implementation -------------------------------------------------------------
/// ----------------------------------------------------------------------------
} // namespace voreen

#include "voreen/core/io/serialization/jsonserializer.h"
#include "voreen/core/io/serialization/xmlserializer.h"

namespace voreen {

#define DISPATCH_TO_SERIALIZER(fn) \
switch(concreteSerializer.xml.type) {\
    case SerializerType::XML:\
        concreteSerializer.xml.serializer->fn ;\
        break;\
    case SerializerType::JSON:\
        concreteSerializer.json.serializer->fn ;\
        break;\
    /* default omitted to emit compiler warning on addition of serializers */\
}\

#define DISPATCH_TO_SERIALIZER_RET(fn) \
switch(concreteSerializer.xml.type) {\
    case SerializerType::XML:\
        return concreteSerializer.xml.serializer->fn ;\
    case SerializerType::JSON:\
        return concreteSerializer.json.serializer->fn ;\
    /* default omitted to emit compiler warning on addition of serializers */\
}\

template<typename T>
void Serializer::optionalSerialize(const std::string& key, const T& data, const T& defaultValue) {
    if(data != defaultValue)
        serialize(key, data);
}

template <class T>
void Serializer::serializeBinaryBlob(const std::string& key, const std::vector<T>& data) {
    serialize(key+".numItems", data.size());
    if (!data.empty())
        serializeBinaryBlob(key+".data", reinterpret_cast<const unsigned char*>(&data[0]), sizeof(T) * data.size());
    else
        serializeBinaryBlob(key+".data", 0, 0);
}


template<typename T>
void Serializer::serialize(const std::string& key, const T& data) {
    DISPATCH_TO_SERIALIZER(serialize(key, data));
}

template<typename T>
void Serializer::serialize(const std::string& key, const tgt::TemplateBounds<T>& data) {
    DISPATCH_TO_SERIALIZER(serialize(key, data));
}

template<class T>
void Serializer::serialize(const std::string& key, const T* const& data) {
    DISPATCH_TO_SERIALIZER(serialize(key, data));
}

template<class S, class T>
void Serializer::serialize(const std::string& key, const std::pair<S, T>& data) {
    DISPATCH_TO_SERIALIZER(serialize(key, data));
}

template<class T>
void Serializer::serialize(const std::string& key, const std::vector<T>& data, const std::string& itemKey) {
    serializeCollection(key, data, itemKey);
}

template<class T>
void Serializer::serialize(const std::string& key, const std::deque<T>& data, const std::string& itemKey) {
    serializeCollection(key, data, itemKey);
}

template<class T>
void Serializer::serialize(const std::string& key, const std::list<T>& data, const std::string& itemKey) {
    serializeCollection(key, data, itemKey);
}

template<class T, class C>
void Serializer::serialize(const std::string& key, const std::set<T, C>& data, const std::string& itemKey) {
    serializeCollection(key, data, itemKey);
}

template<class T, class U, class C>
void Serializer::serialize(const std::string& key,
                              const std::map<T, U, C>& data,
                              const std::string& valueKey,
                              const std::string& keyKey) {
    serializeMap(key, data, valueKey, keyKey);
}

template<class T>
void Serializer::serializeCollection(const std::string& key, const T& collection, const std::string& itemKey) {
    DISPATCH_TO_SERIALIZER(serializeCollection(key, collection, itemKey));
}

template<class T>
void Serializer::serializeMap(const std::string& key,
                                        const T& map,
                                        const std::string& valueKey,
                                        const std::string& keyKey) {
    DISPATCH_TO_SERIALIZER(serializeMap(key, map, valueKey, keyKey));
}


template<class T>
void Serializer::raise(const T& exception) {
    DISPATCH_TO_SERIALIZER(raise(exception));
}


} // namespace voreen
#endif //VRN_SERIALIZER_H
