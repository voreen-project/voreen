/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_JSONSERIALIZER_H
#define VRN_JSONSERIALIZER_H

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

#include "rapidjson/document.h"
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/serializerbase.h"
#include "voreen/core/io/serialization/serializationexceptions.h"
#include "voreen/core/io/serialization/xmlserializationconstants.h"

namespace voreen {

/**
 * @c JsonSerializer is responsible for serializing memory data to json documents.
 *
 * @par
 * The implementation is capable of serializing simple C++ data types, user defined classes
 * which are derived from @c Serializable, and most STL containers containing just the previously
 * mentioned data.
 *
 * @par
 * Furthermore, cycles, joints and polymorphic @c Serializable derivatives are supported which
 * makes the @c JsonSerializer quite flexible. The serialization process is independent
 * of serialization order due to the use of @c key/data pairs, provided that different
 * keys are used at each hierarchy level.
 *
 * @par
 * Finally, the Json document can be written to an arbitrary @c std::ostream. For instance, this
 * allows you to write the Json document into a @c std::stringstream for further processing.
 *
 * @par
 * You have to use @c JsonDeserializer, which is the counterpart to @c JsonSerializer, for
 * deserializing Json documents to memory data.
 *
 * @par
 * Here is a short example of using the @c JsonSerializer:
 * @code
 * int i = 1;
 *
 * std::fstream f;
 * f.open("file.xml", std::ios::out);
 *
 * JsonSerializer s;
 * Serializer serializer(s);
 * serializer.serialize("i", i);
 * s.write(f);
 *
 * f.close();
 * @endcode
 * For more complex examples and interaction with other classes of the serialization framework
 * see the &quot;serializertest&quot; application in &quot;apps/&quot; directory.
 *
 * @attention All memory data must be serialized before the Json document can be written.
 *            Otherwise, it cannot be ensured that all pointer references are correctly resolved.
 *
 * @note The implementation of reference deserialization is still incomplete!
 *
 * @see JsonDeserializer
 * @see Serializable
 */
class VRN_CORE_API JsonSerializer : public SerializerBase {
public:
    /**
     * Constructor which initializes the Json document.
     *
     * @par
     * Initialization of the Json document means creating Json declaration and root node
     * as well as adding them to the Json document.
     *
     * @param documentPath Absolute working directory of the document, which is typically the path
     *      to the Json file the document will be written to. This information is not used by the
     *      serializer itself and is therefore not required, but is intended to be accessible
     *      by serializing objects for absolute-to-relative path conversions.
     */
    JsonSerializer(std::string documentPath = "");

    /**
     * Default destructor.
     */
    ~JsonSerializer();

    /**
     * Writes the Json document that contains all already serialized data to the given stream.
     *
     * @attention Keep in mind that all memory data must be serialized before the Json document
     *            can be written. Otherwise, it cannot be ensured that all pointer references
     *            are resolved.
     *
     * @param stream the output stream
     * @param pretty formats stream with newlines, if set to true
     * @param compressed whether or not the output should be be gzip compressed
     */
    void write(std::ostream& stream, bool pretty=false, bool compressed=false);


    /// ----------------------------------------------------------------------------
    /// Generic *Serializer Interface (used by Serializer for dynamic dispatch) ----
    /// ----------------------------------------------------------------------------

    /**
     * Returns the absolute working directory of the document, which is typically the path
     * to the Json file the document will be written to.
     */
    std::string getDocumentPath() const;

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const bool& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const char& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const signed char& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const unsigned char& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const uint16_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const int16_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const uint32_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const int32_t& data);

// There seems to be no uint*_t typedef for long unsigned ints on mac, so we need to provide an implementation for this type.
#ifdef __APPLE__
    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const long unsigned int& data);
#endif

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const uint64_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const int64_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const float& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const double& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const long double& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    void serialize(const std::string& key, const std::string& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::vec2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::vec3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::vec4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::dvec2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::dvec3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::dvec4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::ivec2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::ivec3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::ivec4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::svec2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::svec3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::svec4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::col3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::col4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::mat2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::mat3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::mat4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::Matrix2d& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::Matrix3d& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::Matrix4d& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the Json node key
     * @param data the data
     */
    template<typename T>
    void serialize(const std::string& key, const tgt::TemplateBounds<T>& data);

    /**
     * Serializes the given @c Serializable interface realization.
     *
     * @note All user defined classes must realize the @c Serializable
     *       interface to be serializable.
     *
     * @param key the Json node key
     * @param data the @c Serializable realization
     *
     */
    void serialize(const std::string& key, const Serializable& data);

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

protected:
    /**
     * Category for logging.
     */
    static const std::string loggerCat_;

private:

    /**
     * Adds given unresolved reference to the unresolved reference vector.
     *
     * @tparam T type of referenced data
     *
     * @param unresolvedReferenceNode Json node of the unresolved reference
     * @param reference the unresolved reference
     */
    //template<class T>
    //void addUnresolvedReference(TiJsonElement* unresolvedReferenceNode, const T* const& reference);

    /**
     * Resolves unresolved references during @c write method.
     */
    void resolveUnresolvedReferences();

    /**
     * Helper function for serializing @c key/data pairs with simple data types.
     *
     * @tparam T data type
     *
     * @param key the Json node key
     * @param data the data
     *
     */
    template<class T>
    void serializeSimpleTypes(const std::string& key, const T& data);

    /**
     * Helper function for serializing @c key/data pairs with tgt vectors.
     *
     * @tparam T data type of vector
     *
     * @param key the Json node key
     * @param vector the tgt vector
     * @param isColor flag determine if given vector is a color
     */
    template<class T>
    void serializeTgtVector(const std::string& key, const T& vector);


    /// Path to the target Json document
    std::string documentPath_;

    rapidjson::Document document_;
    rapidjson::Document::AllocatorType& allocator_; //Allocator for document_, conveniency accessor
    rapidjson::Value* node_; // Invariant: Always an object/document node!

    class VRN_CORE_API TemporaryNodeChanger {
    public:
        TemporaryNodeChanger(JsonSerializer& serializer, rapidjson::Value& node)
            : serializer_(serializer)
            , storedNode_(serializer.node_)
        {
            serializer.node_ = &node;
        }

        ~TemporaryNodeChanger() {
            serializer_.node_ = storedNode_;
        }

    private:
        JsonSerializer& serializer_;

        rapidjson::Value* storedNode_;
    };

    void serializeMember(const std::string& key, rapidjson::Value&& val) {
        rapidjson::Value keyV;
        keyV.SetString(key.c_str(), static_cast<rapidjson::SizeType>(key.length()), allocator_);
        node_->AddMember(keyV, val, allocator_);
    }

    rapidjson::Value serializeToValue(const Serializable& data);

    rapidjson::Value serializeToValue(const std::string& data) {
        rapidjson::Value val;
        val.SetString(data.c_str(), static_cast<rapidjson::SizeType>(data.length()), allocator_);

        return val;
    }

    rapidjson::Value serializeToValue(const bool& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const uint8_t& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const uint16_t& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const uint32_t& data) {
        return serializeJsonScalarToValue(data);
    }
#ifdef __APPLE__
    rapidjson::Value serializeToValue(const unsigned long& data) {
        return serializeJsonScalarToValue(static_cast<uint64_t>(data));
    }
#endif
    rapidjson::Value serializeToValue(const uint64_t& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const int8_t& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const int16_t& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const int32_t& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const int64_t& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const float& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const double& data) {
        return serializeJsonScalarToValue(data);
    }
    rapidjson::Value serializeToValue(const tgt::vec2& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::vec3& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::vec4& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::ivec2& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::ivec3& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::ivec4& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::svec2& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::svec3& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::svec4& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::dvec2& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::dvec3& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::dvec4& data) {
        return serializeTgtVecToValue(data, false);
    }
    rapidjson::Value serializeToValue(const tgt::col3& data) {
        return serializeTgtVecToValue(data, true);
    }
    rapidjson::Value serializeToValue(const tgt::col4& data) {
        return serializeTgtVecToValue(data, true);
    }
    rapidjson::Value serializeToValue(const tgt::mat2& data) {
        return serializeTgtMatToValue(data);
    }
    rapidjson::Value serializeToValue(const tgt::mat3& data) {
        return serializeTgtMatToValue(data);
    }
    rapidjson::Value serializeToValue(const tgt::mat4& data) {
        return serializeTgtMatToValue(data);
    }
    rapidjson::Value serializeToValue(const tgt::dmat2& data) {
        return serializeTgtMatToValue(data);
    }
    rapidjson::Value serializeToValue(const tgt::dmat3& data) {
        return serializeTgtMatToValue(data);
    }
    rapidjson::Value serializeToValue(const tgt::dmat4& data) {
        return serializeTgtMatToValue(data);
    }

    template<typename T>
    rapidjson::Value serializeToValue(const tgt::TemplateBounds<T>& data);

    template<class T>
    rapidjson::Value serializeToValue(const T* const& data);

    template<class S, class T>
    rapidjson::Value serializeToValue(const std::pair<S, T>& data);

    template<class S>
    rapidjson::Value serializeToValue(const std::vector<S>& data);

    template<class T>
    rapidjson::Value serializeToValue(const std::deque<T>& data);

    template<class T>
    rapidjson::Value serializeToValue(const std::list<T>& data);

    template<class T, class C>
    rapidjson::Value serializeToValue(const std::set<T, C>& data);

    template<class T, class U, class C>
    rapidjson::Value serializeToValue(const std::map<T, U, C>& data);

    template<class S>
    rapidjson::Value serializeIterableToValue(const S& iterable);

    template<class S>
    rapidjson::Value serializeMapToValue(const S& map, const std::string& valueKey = XmlSerializationConstants::VALUENODE, const std::string& keyKey = XmlSerializationConstants::KEYNODE);

    template<class T>
    rapidjson::Value serializeJsonScalarToValue(const T& data) {
        rapidjson::Value val(data);
        return val;
    }
    template<class T>
    rapidjson::Value serializeTgtVecToValue(T&, bool isColor);

    template<class T>
    rapidjson::Value serializeTgtMatToValue(T&);
};



template<class T>
rapidjson::Value JsonSerializer::serializeTgtVecToValue(T& vector, bool isColor) {
    rapidjson::Value vec;
    vec.SetArray();

    //TODO: handle isColor? why?
    for(int i=0; i < vector.size; ++i) {
        vec.PushBack(serializeToValue(vector[i]), allocator_);
    }

    return vec;
}

template<class T>
rapidjson::Value JsonSerializer::serializeTgtMatToValue(T& matrix) {
    rapidjson::Value mat;
    mat.SetArray();

    for(int i=0; i < matrix.rows; ++i) {
        mat.PushBack(serializeToValue(matrix.elemRows[i]), allocator_);
    }

    return mat;
}

template<class T>
void JsonSerializer::serialize(const std::string& key, const T* const& data) {
    serializeMember(key, serializeToValue(data));
}

template<class S, class T>
void JsonSerializer::serialize(const std::string& key, const std::pair<S, T>& data) {
    serializeMember(key, serializeToValue(data));
}

template<class T>
void JsonSerializer::serializeSimpleTypes(const std::string& key, const T& data) {
    serializeMember(key, serializeToValue(data));
}

template<class T>
void JsonSerializer::serializeTgtVector(const std::string& key, const T& vector) {
    serializeMember(key, serializeToValue(vector));
}

template<class T>
void JsonSerializer::serializeCollection(const std::string& key, const T& collection, const std::string&) {
    serializeMember(key, serializeToValue(collection));
}

template<class T>
void JsonSerializer::serializeMap(const std::string& key,
                                        const T& map,
                                        const std::string& valueKey,
                                        const std::string& keyKey) {
    serializeMember(key, serializeMapToValue(map, valueKey, keyKey));
}

template<typename T>
void JsonSerializer::serialize(const std::string& key, const tgt::TemplateBounds<T>& data) {
    serializeMember(key, serializeToValue(data));
}

template<typename T>
rapidjson::Value JsonSerializer::serializeToValue(const tgt::TemplateBounds<T>& data) {
    std::pair<tgt::Vector3<T>, tgt::Vector3<T>> pair(data.getLLF(), data.getURB());
    return serializeToValue(pair);
}

template<class T>
rapidjson::Value JsonSerializer::serializeToValue(const T* const& data) {
    /*
    //TODO in general:
    rapidjson::Value val;
    //TODO set appropriate type

    // Is data not a null pointer?
    if (data) {
        //addUnresolvedReference(newNode, data);
    }
    return val;
    */
    throw SerializationException("Pointer/reference serialization is not supported (yet)");
}

template<class S, class T>
rapidjson::Value JsonSerializer::serializeToValue(const std::pair<S, T>& data) {
    rapidjson::Value pair;
    pair.SetObject();

    pair.AddMember("First", serializeToValue(data.first), allocator_);

    pair.AddMember("Second", serializeToValue(data.second), allocator_);

    return pair;
}

template<class S>
rapidjson::Value JsonSerializer::serializeToValue(const std::vector<S>& data) {
    return serializeIterableToValue(data);
}

template<class T>
rapidjson::Value JsonSerializer::serializeToValue(const std::deque<T>& data) {
    return serializeIterableToValue(data);
}

template<class T>
rapidjson::Value JsonSerializer::serializeToValue(const std::list<T>& data) {
    return serializeIterableToValue(data);
}

template<class T, class C>
rapidjson::Value JsonSerializer::serializeToValue(const std::set<T, C>& data) {
    return serializeIterableToValue(data);
}

template<class T, class U, class C>
rapidjson::Value JsonSerializer::serializeToValue(const std::map<T, U, C>& data) {
    return serializeMapToValue(data);
}

template<class S>
rapidjson::Value JsonSerializer::serializeIterableToValue(const S& iterable)
{
    rapidjson::Value array;
    array.SetArray();

    for(const auto& element : iterable) {
        array.PushBack(serializeToValue(element), allocator_);
    }

    return array;
}

template<class S>
rapidjson::Value JsonSerializer::serializeMapToValue(const S& map, const std::string& valueKey, const std::string& keyKey)
{
    //aka, serialize iterable of pairs
    rapidjson::Value array;
    array.SetArray();

    for(const auto& element : map) {
        rapidjson::Value pair;
        pair.SetObject();

        pair.AddMember(serializeToValue(keyKey), serializeToValue(element.first), allocator_);
        pair.AddMember(serializeToValue(valueKey), serializeToValue(element.second), allocator_);

        array.PushBack(pair, allocator_);
    }

    return array;
}

} // namespace

#endif // VRN_JSONSERIALIZER_H
