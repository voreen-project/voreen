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

#ifndef VRN_JSONDESERIALIZER_H
#define VRN_JSONDESERIALIZER_H

#include <string>
#include <utility>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <iostream>
#include <typeinfo>

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"
#include "tgt/logmanager.h"

#include "voreen/core/io/serialization/xmlserializationconstants.h"
#include "voreen/core/io/serialization/serializerbase.h"
#include "voreen/core/io/serialization/serializationexceptions.h"
#include "voreen/core/io/serialization/serializable.h"
#include "rapidjson/document.h"

namespace voreen {

/**
 * @c JsonDeserializer is responsible for deserializing Json documents to memory data.
 *
 * @par
 * The implementation is capable of deserializing simple C++ data types, user defined classes
 * which are derived from @c Serializable, and most STL containers containing just the previously
 * mentioned data.
 *
 * @par
 * Furthermore, cycles, joints and polymorphic @c Serializable derivatives are supported which
 * makes the @c JsonDeserializer quite flexible. The deserialization process is independent
 * of deserialization order due to the use of @c key/value data, provided that different
 * keys are used at each hierarchy level.
 *
 * @par
 * Initially, the Json document can be read from an arbitrary @c std::istream. For instance, this
 * allows you to read a in memory created Json document from a @c std::stringstream.
 *
 * @par
 * You have to use @c JsonSerializer, which is the counterpart to @c JsonDeserializer, for
 * serializing memory data to Json documents.
 *
 * @par
 * Here is a short example of using the @c JsonDeserializer:
 * @code
 * int i;
 *
 * std::fstream f;
 * f.open("file.xml", std::ios::in);
 *
 * JsonDeserializer d;
 * d.read(f);
 * Deserializer deserializer(d);
 * deserializer.deserialize("i", i);
 *
 * f.close();
 * @endcode
 * For more complex examples and interaction with other classes of the serialization framework
 * see the &quot;serializertest&quot; application in &quot;apps/&quot; directory.
 *
 * @attention The whole Json document must be deserialized before you can use deserialized data.
 *            Otherwise, it cannot be ensured that all pointer references are correctly resolved.
 *
 * @note The implementation of reference deserialization is still incomplete!
 *
 * @see JsonSerializer
 * @see Serializer
 * @see SerializerBase
 * @see Serializable
 */
class VRN_CORE_API JsonDeserializer : public SerializerBase {
public:
    /**
     * Constructor.
     *
     * @param documentPath Absolute working directory of the document, which is typically the path
     *      to the Json file the document has been read from. This information is not used by the
     *      serializer itself and is therefore not required, but is intended to be accessible
     *      by serializing objects for relative-to-absolute path conversions.
     */
    JsonDeserializer(std::string documentPath = "");

    /**
     * Default destructor.
     */
    ~JsonDeserializer();

    /**
     * Reads the Json document from the given input stream after an optional Json preprocessor is applied.
     *
     * @param stream the input stream
     *
     * @throws XmlSerializationFormatException if the Json document is incorrect formatted.
     * @throws XmlSerializationVersionMismatchException if the Json document serialization version
     *     does not match with the serialization version of this class.
     * @throws XmlSerializationReferenceResolvingException if there are references in
     *     the Json document which cannot be resolved.
     */
    void read(std::istream& stream);


    /// ----------------------------------------------------------------------------
    /// Generic *Deserializer Interface (used by Deserializer for dynamic dispatch)
    /// ----------------------------------------------------------------------------

    /**
     * Returns the absolute working directory of the document, which is typically the path
     * to the Json file the document has been read from.
     */
    std::string getDocumentPath() const;

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
    void deserialize(const std::string& key, bool& data);

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
    void deserialize(const std::string& key, char& data);

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
    void deserialize(const std::string& key, signed char& data);

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
    void deserialize(const std::string& key, unsigned char& data);

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
    void deserialize(const std::string& key, uint16_t& data);

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
    void deserialize(const std::string& key, int16_t& data);

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
    void deserialize(const std::string& key, uint32_t& data);

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
    void deserialize(const std::string& key, int32_t& data);

// There seems to be no uint*_t typedef for long unsigned ints on mac, so we need to provide an implementation for this type.
#ifdef __APPLE__
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
    void deserialize(const std::string& key, long unsigned int& data);
#endif

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
    void deserialize(const std::string& key, uint64_t& data);

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
    void deserialize(const std::string& key, int64_t& data);

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
    void deserialize(const std::string& key, float& data);

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
    void deserialize(const std::string& key, double& data);

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
    void deserialize(const std::string& key, long double& data);

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
    void deserialize(const std::string& key, std::string& data);

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
    void deserialize(const std::string& key, tgt::vec2& data);

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
    void deserialize(const std::string& key, tgt::vec3& data);

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
    void deserialize(const std::string& key, tgt::vec4& data);

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
    void deserialize(const std::string& key, tgt::dvec2& data);

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
    void deserialize(const std::string& key, tgt::dvec3& data);

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
    void deserialize(const std::string& key, tgt::dvec4& data);

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
    void deserialize(const std::string& key, tgt::ivec2& data);

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
    void deserialize(const std::string& key, tgt::ivec3& data);

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
    void deserialize(const std::string& key, tgt::ivec4& data);

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
    void deserialize(const std::string& key, tgt::col3& data);

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
    void deserialize(const std::string& key, tgt::col4& data);

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
    void deserialize(const std::string& key, tgt::mat2& data);

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
    void deserialize(const std::string& key, tgt::mat3& data);

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
    void deserialize(const std::string& key, tgt::mat4& data);

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
    void deserialize(const std::string& key, tgt::Matrix2d& data);

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
    void deserialize(const std::string& key, tgt::Matrix3d& data);

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
    void deserialize(const std::string& key, tgt::Matrix4d& data);

    /**
     * Deserializes the given @c Serializable interface realization.
     *
     * @note All user defined classes must realize the @c Serializable
     *       interface to be deserializable.
     *
     * @param key the Json node key
     * @param data variable to store deserialized @c Serializable realization
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a Json node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple Json nodes share same id attribute
     */
    void deserialize(const std::string& key, Serializable& data);

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
     * @param key the Json node key
     * @param collection variable to store deserialized data collection
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
     * @param key the Json node key
     * @param collection variable to store deserialized data collection
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
     * @param key the Json node key
     * @param map varaible to store deserialized data map
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
    template<class T>
    void deserializeMap(
        const std::string& key,
        T& map,
        const std::string& valueKey = XmlSerializationConstants::VALUENODE,
        const std::string& keyKey = XmlSerializationConstants::KEYNODE);

protected:
    /**
     * Category for logging.
     */
    static const std::string loggerCat_;

private:


    class VRN_CORE_API TemporaryNodeChanger {
    public:
        TemporaryNodeChanger(JsonDeserializer& deserializer, const rapidjson::Value& node)
            : deserializer_(deserializer)
            , storedNode_(deserializer.node_)
        {
            deserializer.node_ = &node;
        }

        ~TemporaryNodeChanger() {
            deserializer_.node_ = storedNode_;
        }

    private:
        JsonDeserializer& deserializer_;

        const rapidjson::Value* storedNode_;
    };

    /// Path to Json file the document was read from.
    std::string documentPath_;

    rapidjson::Document document_;
    std::vector<char> documentBuffer_;
    const rapidjson::Value* node_;

    /// ----------------------------------------------------------------------------
    /// Some Implementation --------------------------------------------------------
    /// ----------------------------------------------------------------------------
    template<class T>
    void readInt(const rapidjson::Value& val, T& data);
    template<class T>
    void readUint(const rapidjson::Value& val, T& data);
    template<class T>
    void readFloat(const rapidjson::Value& val, T& data);
    template<class T>
    void readTgtVector(const rapidjson::Value& val, T& data);
    template<class T>
    void readTgtMatrix(const rapidjson::Value& val, T& data);



    void readFromValue(const rapidjson::Value& val, bool& data) {
        if(!val.IsBool()) {
            throw SerializationException("Not a bool");
        } else {
            data = val.GetBool();
        }
    }
    void readFromValue(const rapidjson::Value& val, char& data) { //TODO: maybe remove char? we don't know if it's signed!
        readInt(val, data);
    }
    void readFromValue(const rapidjson::Value& val, uint8_t& data) {
        readUint(val, data);
    }
    void readFromValue(const rapidjson::Value& val, uint16_t& data) {
        readUint(val, data);
    }
    void readFromValue(const rapidjson::Value& val, uint32_t& data) {
        readUint(val, data);
    }
    void readFromValue(const rapidjson::Value& val, uint64_t& data) {
        readUint(val, data);
    }
    void readFromValue(const rapidjson::Value& val, int8_t& data) {
        readInt(val, data);
    }
    void readFromValue(const rapidjson::Value& val, int16_t& data) {
        readInt(val, data);
    }
    void readFromValue(const rapidjson::Value& val, int32_t& data) {
        readInt(val, data);
    }
    void readFromValue(const rapidjson::Value& val, int64_t& data) {
        readInt(val, data);
    }
    void readFromValue(const rapidjson::Value& val, float& data) {
        readFloat(val, data);
    }
    void readFromValue(const rapidjson::Value& val, double& data) {
        readFloat(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::vec2& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::vec3& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::vec4& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::dvec2& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::dvec3& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::dvec4& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::ivec2& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::ivec3& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::ivec4& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::col3& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::col4& data) {
        readTgtVector(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::mat2& data) {
        readTgtMatrix(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::mat3& data) {
        readTgtMatrix(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::mat4& data) {
        readTgtMatrix(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::dmat2& data) {
        readTgtMatrix(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::dmat3& data) {
        readTgtMatrix(val, data);
    }
    void readFromValue(const rapidjson::Value& val, tgt::dmat4& data) {
        readTgtMatrix(val, data);
    }

    void readFromValue(const rapidjson::Value& val, std::string& data) {
        if(!val.IsString()) {
            throw SerializationException("Not a string");
        } else {
            data = val.GetString();
        }
    }
    void readFromValue(const rapidjson::Value& val, Serializable& data);

    template<class T>
    void readFromValue(const rapidjson::Value& val, T*& data) {
        throw SerializationException("Pointer Deserialization is not (yet) implemented");
    }

    template<typename T>
    void readFromValue(const rapidjson::Value& val, tgt::TemplateBounds<T>& data);

    template<class S, class T>
    void readFromValue(const rapidjson::Value& val, std::pair<S, T>& data);

    template<class T>
    void readFromValue(const rapidjson::Value& val, std::vector<T>& data);

    template<class T>
    void readFromValue(const rapidjson::Value& val, std::deque<T>& data);

    template<class T>
    void readFromValue(const rapidjson::Value& val, std::list<T>& data);

    template<class T, class C>
    void readFromValue(const rapidjson::Value& val, std::set<T, C>& data);

    template<class T, class U, class C>
    void readFromValue(const rapidjson::Value& val, std::map<T, U, C>& data);

    template<class T>
    void readFromValueIntoBackInserter(const rapidjson::Value& val, T& data);

    template<class T>
    void readFromValueIntoMap(const rapidjson::Value& val, T& data, const std::string& valueKey = XmlSerializationConstants::VALUENODE, const std::string& keyKey = XmlSerializationConstants::KEYNODE);

    const rapidjson::Value& getFromValue(const rapidjson::Value& node, const std::string& key) {
        if(!node.IsObject()) {
            throw SerializationException("Node is not an object.");
        }
        if(!node.HasMember(key.c_str())) {
            throw SerializationException("Node does not have key \"" + key + "\".");
        }
        return node[key.c_str()];
    }

    const rapidjson::Value& getMember(const std::string& key) {
        return getFromValue(*node_, key);
    }

    template<typename T>
    void deserializeFromMember(const std::string& key, T& val) {
        readFromValue(getMember(key), val);
    }
};

} //namespace voreen

#include "voreen/core/io/serialization/deserializer.h"

namespace voreen {

template<class T>
void JsonDeserializer::readInt(const rapidjson::Value& val, T& data) {
    if(!val.IsInt()) {
        throw SerializationException("Not an int");
    } else {
        data = static_cast<T>(val.GetInt());
    }
}
template<class T>
void JsonDeserializer::readUint(const rapidjson::Value& val, T& data) {
    if(!val.IsUint()) {
        throw SerializationException("Not an int");
    } else {
        data = static_cast<T>(val.GetUint());
    }
}

template<class T>
void JsonDeserializer::readFloat(const rapidjson::Value& val, T& data) {
    if(!val.IsDouble()) {
        throw SerializationException("Not a float");
    } else {
        data = static_cast<T>(val.GetDouble());
    }
}

template<class T>
void JsonDeserializer::readTgtVector(const rapidjson::Value& val, T& data) {
    if(!val.IsArray()) {
        throw SerializationException("Not an array");
    }
    if(val.Size() != data.size) {
        throw SerializationException("Invalid array size");
    }
    for(int i=0; i < data.size; ++i) {
        readFromValue(val[i], data[i]);
    }
}

template<class T>
void JsonDeserializer::readTgtMatrix(const rapidjson::Value& val, T& data) {
    if(!val.IsArray()) {
        throw SerializationException("Not an array");
    } if(val.Size() != data.cols) {
        throw SerializationException("Invalid array size");
    }
    for(int i=0; i<data.cols; ++i) {
         readFromValue(val[i], data.elemRows[i]);
    }
}

template<class T>
void JsonDeserializer::deserialize(const std::string& key, T*& data) {
    //TODO
    throw SerializationException("Not implemented");
}

template<class S, class T>
void JsonDeserializer::deserialize(const std::string& key, std::pair<S, T>& data) {
    deserializeFromMember(key, data);
}

template<class T>
void JsonDeserializer::deserializeCollectionConstIteratorItems(const std::string& key, T& collection, const std::string& /*itemKey (not meaningful for json) */) {
    deserializeFromMember(key, collection);
}

template<class T>
void JsonDeserializer::deserializeCollection(const std::string& key, T& collection, const std::string& itemKey) {
    deserializeCollectionConstIteratorItems(key, collection, itemKey);
}

template<class T>
void JsonDeserializer::deserializeMap(const std::string& mapItemKey,
                                            T& map,
                                            const std::string& valueKey,
                                            const std::string& keyKey) {
    deserializeFromMember(mapItemKey, map);
}

template<typename T>
void JsonDeserializer::deserialize(const std::string& key, tgt::TemplateBounds<T>& data) {
    deserializeFromMember(key, data);
}

template<typename T>
void JsonDeserializer::readFromValue(const rapidjson::Value& val, tgt::TemplateBounds<T>& data) {
    //TODO fix in serialize
    std::pair<tgt::Vector3<T>, tgt::Vector3<T>> llfAndUrb;
    readFromValue(val, llfAndUrb);
    data = tgt::TemplateBounds<T>(llfAndUrb.first, llfAndUrb.second);
}

template<class S, class T>
void JsonDeserializer::readFromValue(const rapidjson::Value& val, std::pair<S, T>& data) {
    TemporaryNodeChanger nodeChanger(*this, val);
    Deserializer(*this).deserialize("First", data.first);
    Deserializer(*this).deserialize("Second", data.second);
}

template<class T>
void JsonDeserializer::readFromValue(const rapidjson::Value& val, std::vector<T>& data) {
    readFromValueIntoBackInserter(val, data);
}

template<class T>
void JsonDeserializer::readFromValue(const rapidjson::Value& val, std::deque<T>& data) {
    readFromValueIntoBackInserter(val, data);
}

template<class T>
void JsonDeserializer::readFromValue(const rapidjson::Value& val, std::list<T>& data) {
    readFromValueIntoBackInserter(val, data);
}

template<class T, class C>
void JsonDeserializer::readFromValue(const rapidjson::Value& val, std::set<T, C>& data) {
    readFromValueIntoBackInserter(val, data);
}

template<class T, class U, class C>
void JsonDeserializer::readFromValue(const rapidjson::Value& val, std::map<T, U, C>& data) {
    readFromValueIntoMap(val, data);
}

template<class T>
void JsonDeserializer::readFromValueIntoBackInserter(const rapidjson::Value& val, T& collection) {
    if(!val.IsArray()) {
        throw SerializationException("Value is not an Array.");
    }

    for (rapidjson::Value::ConstValueIterator itr = val.Begin(); itr != val.End(); ++itr) {
        const auto& element = *itr;

        typename T::value_type item;

        readFromValue(element, item);

        collection.insert(collection.end(), std::move(item));
    }
}
template<class T>
void JsonDeserializer::readFromValueIntoMap(const rapidjson::Value& val, T& map, const std::string& valueKey, const std::string& keyKey) {

    if(!val.IsArray()) {
        throw SerializationException("Value is not an Array.");
    }

    for (rapidjson::Value::ConstValueIterator itr = val.Begin(); itr != val.End(); ++itr) {
        const auto& element = *itr;

        typename T::key_type key;
        typename T::mapped_type value;

        TemporaryNodeChanger nodeChanger(*this, element);
        Deserializer(*this).deserialize(valueKey, key);
        Deserializer(*this).deserialize(keyKey, value);

        map[key] = value;
    }
}

} // namespace

#endif // VRN_JSONDESERIALIZER_H
