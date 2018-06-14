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

#ifndef VRN_XMLDESERIALIZER_H
#define VRN_XMLDESERIALIZER_H

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

#include "tinyxml/tinyxml.h"

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"
#include "tgt/logmanager.h"

#include "voreen/core/io/serialization/xmlserializerbase.h"
#include "voreen/core/io/serialization/serializationexceptions.h"
#include "voreen/core/io/serialization/xmlserializationconstants.h"
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/abstractserializable.h"
#include "voreen/core/io/serialization/serializablefactory.h"
#include "voreen/core/io/serialization/xmlprocessor.h"

namespace voreen {

/**
 * @c XmlDeserializer is responsible for deserializing XML documents to memory data.
 *
 * @par
 * The implementation is capable of deserializing simple C++ data types, user defined classes
 * which are derived from @c Serializable, and most STL containers containing just the previously
 * mentioned data.
 *
 * @par
 * Furthermore, cycles, joints and polymorphic @c Serializable derivatives are supported which
 * makes the @c XmlDeserializer quite flexible. The deserialization process is independent
 * of deserialization order due to the use of @c key/value data, provided that different
 * keys are used at each hierarchy level.
 *
 * @par
 * Initially, the XML document can be read from an arbitrary @c std::istream. For instance, this
 * allows you to read a in memory created XML document from a @c std::stringstream.
 *
 * @par
 * You have to use @c XmlSerializer, which is the counterpart to @c XmlDeserializer, for
 * serializing memory data to XML documents.
 *
 * @par
 * Here is a short example of using the @c XmlDeserializer:
 * @code
 * int i;
 *
 * std::fstream f;
 * f.open("file.xml", std::ios::in);
 *
 * XmlDeserializer d;
 * d.read(f);
 * Deserializer deserializer(d);
 * deserializer.deserialize("i", i);
 *
 * f.close();
 * @endcode
 * For more complex examples and interaction with other classes of the serialization framework
 * see the &quot;serializertest&quot; application in &quot;apps/&quot; directory.
 *
 * @attention The whole XML document must be deserialized before you can use deserialized data.
 *            Otherwise, it cannot be ensured that all pointer references are correctly resolved.
 *
 * @note For further information on handling cycles, joints and polymorphism, see:
 *       http://www.parashift.com/c++-faq-lite/serialization.html
 *
 * @see Deserializer
 * @see XmlSerializer
 * @see XmlSerializerBase
 * @see Serializable
 */
class VRN_CORE_API XmlDeserializer : public XmlSerializerBase
{
public:
    /**
     * Constructor.
     *
     * @param documentPath Absolute working directory of the document, which is typically the path
     *      to the XML file the document has been read from. This information is not used by the
     *      serializer itself and is therefore not required, but is intended to be accessible
     *      by serializing objects for relative-to-absolute path conversions.
     */
    XmlDeserializer(std::string documentPath = "");

    /**
     * Default destructor.
     */
    ~XmlDeserializer();

    /**
     * Reads the XML document from the given input stream after an optional XML preprocessor is applied.
     *
     * @param stream the input stream
     * @param xmlProcessor XML preprocessor
     *
     * @throws SerializationFormatException if the XML document is incorrect formatted.
     * @throws SerializationVersionMismatchException if the XML document serialization version
     *     does not match with the serialization version of this class.
     * @throws SerializationReferenceResolvingException if there are references in
     *     the XML document which cannot be resolved.
     */
    void read(std::istream& stream, XmlProcessor* xmlProcessor = 0);


    /// ----------------------------------------------------------------------------
    /// Generic *Deserializer Interface (used by Deserializer for dynamic dispatch)
    /// ----------------------------------------------------------------------------

    /**
     * Returns the absolute working directory of the document, which is typically the path
     * to the XML file the document has been read from.
     */
    std::string getDocumentPath() const;

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, bool& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, char& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, signed char& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, unsigned char& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, uint16_t& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, int16_t& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, uint32_t& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, int32_t& data);

// There seems to be no uint*_t typedef for long unsigned ints on mac, so we need to provide an implementation for this type.
#ifdef __APPLE__
    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, long unsigned int& data);
#endif

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, uint64_t& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, int64_t& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, float& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, double& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, long double& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, std::string& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::vec2& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::vec3& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::vec4& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::dvec2& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::dvec3& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::dvec4& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::ivec2& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::ivec3& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::ivec4& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a XML node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::svec2& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a XML node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::svec3& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws XmlSerializationNoSuchDataException if no data with the given key can be found.
     * @throws XmlSerializationFormatException if a XML node is incorrect formatted.
     * @throws XmlSerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::svec4& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::col3& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::col4& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::mat2& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::mat3& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::mat4& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::Matrix2d& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::Matrix3d& data);

    /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, tgt::Matrix4d& data);

     /**
     * Deserializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    template<typename T>
    void deserialize(const std::string& key, tgt::TemplateBounds<T>& data);

    /**
     * Deserializes the given @c Serializable interface realization.
     *
     * @note All user defined classes must realize the @c Serializable
     *       interface to be deserializable.
     *
     * @param key the XML node key
     * @param data variable to store deserialized @c Serializable realization
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    void deserialize(const std::string& key, Serializable& data);

    /**
     * Deserializes the given pointer reference.
     *
     * @tparam type of referenced data
     *
     * @param key the XML node key
     * @param data variable to store deserialized pointer reference
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationMemoryAllocationException
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
     * @param key the XML node key
     * @param data variable to store deserialized pair
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws SerializationMemoryAllocationException
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
     * @param key the XML node key
     * @param collection variable to store deserialized data collection
     * @param itemKey XML node key for each XML child node
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    inline void deserializeCollectionConstIteratorItems(
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
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    inline void deserializeCollection(
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
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    inline void deserializeMap(
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
    /**
     * Adds given @c id/address pair to the pointer address look up map.
     *
     * @param element XML element containing the id attribute
     * @param address the address
     *
     * @throws SerializationDuplicateIdException if id already exists
     */
    inline void addReferenceAddress(TiXmlElement* element, void* address);

    /**
     * Adds given @c id/reference pair to the unresolved reference map.
     *
     * @param element XML element containing the reference attribute
     * @param data the unresolved reference
     */
    template<class T>
    inline void addUnresolvedReference(TiXmlElement* element, T*& data);

    /**
     * Helper function allocating new memory of given type on the heap using a registered factory.
     *
     * @tparam data type of memory to allocate
     *
     * @param type parameter is needed to allocate memory using a factory
     *
     * @return pointer to allocated memory or @c null in case of memory was not allocatable
     */
    template<class T>
    inline T* allocateMemory(const std::string* type);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(bool* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @param data parameter is needed to determine how to allocate memory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(char* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(signed char* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(unsigned char* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(short* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(int* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(long* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(float* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(double* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(long double* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(std::string* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::vec2* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::vec3* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::vec4* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::ivec2* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::ivec3* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::ivec4* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::col3* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(tgt::col4* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory or @c null in case of memory was not allocatable
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(Serializable* data, const std::string* type = 0);

    /**
     * Allocates new memory of given type on the heap and returns a pointer to it.
     *
     * @tparam data type of memory to allocate
     *
     * @param data parameter is needed to determine how to allocate memory
     * @param type parameter is needed to allocate memory for an @c AbstractSerializable by factory
     *
     * @return pointer to allocated memory or @c null in case of memory was not allocatable
     *
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline T* allocateMemory(AbstractSerializable* data, const std::string* type = 0);

    typedef std::set<std::string> ReferenceIdSetType;
    typedef std::vector<std::string> ReferenceIdListType;

    /**
     * Helper function for recursive search for unresolvable references.
     *
     * To determine all unresolvable references you have to look up for
     * all @c references which are not in the @c resolvableReferences set.
     *
     * @param node the node to search for unresolvable references
     * @param references a set containing all references
     * @param resolvableReferences a set containing all resolvable references
     */
    void findUnresolvableReferences(TiXmlElement* node, ReferenceIdSetType& references, ReferenceIdSetType& resolvableReferences) const;

    /**
     * Finds and returns unresolvable references.
     *
     * @return a vector containing the reference ids of each unresolvable reference.
     */
    ReferenceIdListType findUnresolvableReferences() const;

    /**
     * Helper function for deserializing @c key/data pairs with simple data types.
     *
     * @tparam T data type
     *
     * @param key the XML node key
     * @param data variable to store deserialized data
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     */
    template<class T>
    inline void deserializeSimpleTypes(const std::string& key, T& data);

    /**
     * Helper function for deserializing tgt vectors.
     *
     * @tparam T data type of vector
     *
     * @param key the XML node key
     * @param vector variable to store deserialized vector
     * @param isColor flag determine if given vector is a color
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     */
    template<class T>
    inline void deserializeTgtVector(const std::string& key, T& vector, const bool& isColor = false);

    /**
     * Helper function deserializes the given @c key/data pair from XML attribute named @c key
     * contained in the XML node @c nodeKey.
     *
     * @tparam T data type
     *
     * @param nodeKey the XML node key
     * @param key the XML attribute key
     * @param data variable to store deserialized data
     * @param markVisited if @c true the XML node with given @c nodeKey is marked as visited
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     * @throws SerializationFormatException if a XML node is incorrect formatted.
     * @throws SerializationDuplicateIdException if multiple XML nodes share same id attribute
     * @throws SerializationMemoryAllocationException
     *         in case of trying to allocate memory for an @c AbstractSerializable
     *         or if there are not enough allocate items if pointer content serialization is enabled
     */
    template<class T>
    inline void deserializeAttributeFromNode(const std::string& nodeKey, const std::string& key, T& data, const bool& markVisited = true);

    /**
     * Resolves the reference corresponding with the reference attribute of the given @c element.
     *
     * @param element XML element containing reference attribute
     * @param data variable to store resolved reference
     */
    template<class T>
    void resolveReference(TiXmlElement* element, T*& data);

    /**
     * Returns the next XML element with given @c key as name, which is not deserialized yet.
     *
     * @param key XML node name
     *
     * @returns the next XML element with given @c key as name
     *
     * @throws SerializationNoSuchDataException if no data with the given key can be found.
     */
    TiXmlElement* getNextXmlElement(const std::string& key);

    /**
     * Type definition for set of visited XML nodes.
     */
    typedef std::set<TiXmlNode*> VisitedNodesSetType;

    /**
     * Set of already visited XML nodes.
     */
    VisitedNodesSetType visitedNodes_;

    /**
     * Type definition for pointer address look up map.
     */
    typedef std::map<std::string, void*> IdAddressMapType;

    /**
     * Map for looking up pointer addresses by reference id.
     */
    IdAddressMapType idAddressMap_;

    /**
     * The @c IReferenceResolver interface allow to treat every @c ReferenceResolver
     * template specialization in the same way.
     *
     * @see ReferenceResolver
     */
    class IReferenceResolver {
    public:

        virtual ~IReferenceResolver() {}

        /**
         * Resolves an unresolved pointer reference.
         *
         * @param resolvedAddress the resolved address
         */
        virtual void resolve(void* resolvedAddress) = 0;
    };

    /**
     * @c ReferenceResolver is a concrete implementation for resolving
     * unresolved pointer references.
     *
     * @attention Resolving of unresolved references cannot be done in
     *            case of incomplete XML documents, which does not contain
     *            the referenced data, or in case of deserialization method-calls
     *            do not correspond to the serialization method-calls.
     *
     * @tparam type of referenced data
     *
     * @see IReferenceResolver
     */
    template<class T>
    class ReferenceResolver : public IReferenceResolver {
    public:
        /**
         * Constructs a @c ReferenceResolver in accordance to the given parameters.
         *
         * @param data variable to store resolved address
         */
        ReferenceResolver(T*& data);

        /**
         * @see IReferenceResolver::resolve
         */
        virtual void resolve(void* resolvedAddress);

    private:
        /**
         * Reference to pointer which needs address resolving.
         */
        T*& data_;
    };

    /**
     * Type definition for reference resolver list.
     */
    typedef std::vector<IReferenceResolver*> ReferenceResolverListType;

    /**
     * Type definition for unresolved reference map.
     */
    typedef std::map<std::string, ReferenceResolverListType> UnresolvedReferenceMapType;

    /**
     * Unresolved reference map which contains map/reference-resolver-list pairs.
     */
    UnresolvedReferenceMapType unresolvedReferenceMap_;

    /// Path to XML file the document was read from.
    std::string documentPath_;
};

} // namespace voreen

#include "voreen/core/io/serialization/deserializer.h"

namespace voreen {

template<class T>
XmlDeserializer::ReferenceResolver<T>::ReferenceResolver(T*& data) : data_(data) {
}

template<class T>
void XmlDeserializer::ReferenceResolver<T>::resolve(void* resolvedAddress) {
    data_ = static_cast<T*>(resolvedAddress);
}

inline void XmlDeserializer::addReferenceAddress(TiXmlElement* element, void* address) {
    const std::string* id = element->Attribute(XmlSerializationConstants::IDATTRIBUTE);
    // Has element no id attribute?
    if (!id)
        return;

    // Does id already exists?
    if (idAddressMap_.find(*id) != idAddressMap_.end())
        raise(SerializationDuplicateIdException("Multiple XML nodes share the id attribute '" + *id + "'"));

    // If there are unresolved references to this id, resolve them...
    UnresolvedReferenceMapType::iterator idIterator = unresolvedReferenceMap_.find(*id);
    if (idIterator != unresolvedReferenceMap_.end()) {
        for (ReferenceResolverListType::iterator it = idIterator->second.begin(); it != idIterator->second.end(); ++it) {
            (*it)->resolve(address);
            delete *it;
        }

        unresolvedReferenceMap_.erase(idIterator);
    }

    idAddressMap_[*id] = address;
}

template<class T>
inline void XmlDeserializer::addUnresolvedReference(TiXmlElement* element, T*& data) {
    const std::string* reference = element->Attribute(XmlSerializationConstants::REFERENCEATTRIBUTE);

    // Has element no reference attribute?
    if (!reference)
        return;

    // Are there other unresolved references with same id?
    if (unresolvedReferenceMap_.find(*reference) != unresolvedReferenceMap_.end()) {
        // Insert reference resolver into existing list...
        unresolvedReferenceMap_.find(*reference)->second.push_back(new ReferenceResolver<T>(data));
    }
    else {
        // Create new reference resolver list and append a new reference resolver...
        ReferenceResolverListType referenceResolverList;
        referenceResolverList.push_back(new ReferenceResolver<T>(data));
        unresolvedReferenceMap_[*reference] = referenceResolverList;
    }
}

template<class T>
inline T* XmlDeserializer::allocateMemory(const std::string* type) {
    T* data = 0;

    // Does type attribute exists?
    if (type)
        // Look for a matching factory to create object instance...
        for (FactoryListType::iterator it = factories_.begin(); it != factories_.end(); ++it) {
            try {
                Serializable* serializable = (*it)->createSerializableType(*type);
                data = dynamic_cast<T*>(serializable);
            }
            catch (std::exception& e) {
                raise(SerializationMemoryAllocationException("Object of type '" + *type + "' could not be created: " + std::string(e.what())));
            }
            // Was factory capable of instance creation?
            if (data)
                break;
        }

    return data;
}

template<class T> inline T* XmlDeserializer::allocateMemory(bool*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(char*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(signed char*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(unsigned char*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(short*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(int*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(long*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(float*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(double*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(long double*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(std::string*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::vec2*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::vec3*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::vec4*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::ivec2*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::ivec3*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::ivec4*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::col3*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(tgt::col4*, const std::string*) { return new T(); }
template<class T> inline T* XmlDeserializer::allocateMemory(Serializable*, const std::string* type) {
    // Is a factory registered to allocate memory for objects of this type?
    T* data = allocateMemory<T>(type);
    if (data) {
        return data;
    }
    // No factory is registered to allocate memory for this object...
    else {
        // Use default constructor to allocate memory...
        return new T();
    }
}

template<class T> inline T* XmlDeserializer::allocateMemory(AbstractSerializable*, const std::string* type) {
    // Is a factory registered to allocate memory for objects of this type?
    T* data = allocateMemory<T>(type);
    if (data) {
        return data;
    }
    // No factory is registered to allocate memory for this object...
    else {
        // Raise exception since we cannot allocate memory for abstract classes without a factory...
        if (type)
            raise(SerializationMemoryAllocationException("Object of type '" + *type + "' could not be created (no matching factory)."));
        else
            raise(SerializationMemoryAllocationException("Object of type '" + std::string(typeid(T).name()) + "' could not be created (no matching factory)."));
        return 0;
    }
}

template<class T>
void XmlDeserializer::deserialize(const std::string& key, T*& data) {
    // Always deserialize pointer from content?
    if (usePointerContentSerialization_) {
        TemporaryUsePointerContentSerializationChanger usePointerContentSerializationChanger(*this, false);

        Deserializer(*this).deserialize(key, *data);

        return;
    }

    TiXmlElement* element = getNextXmlElement(key);

    const std::string* reference = element->Attribute(XmlSerializationConstants::REFERENCEATTRIBUTE);

    // Set pointer to null to be on the safe side(delete on 0 is allowed)...
    data = 0;

    // Has element reference attribute?
    if (reference) {
        resolveReference(element, data);
    }
    else { // No reference attribute...
        // Does element contain serialized content?
        if (useAttributes_ || element->Attribute(XmlSerializationConstants::VALUEATTRIBUTE) || !element->NoChildren()) {
            // Is AbstractSerializable 0 pointer?
            if (!element->FirstAttribute() && element->NoChildren())
                return;

            // Gather type information...
            const std::string* type = element->Attribute(XmlSerializationConstants::TYPEATTRIBUTE);

            // Allocate memory for the object...
            data = allocateMemory<T>(data, type);

            // ATTENTION: We have to set element unvisited in order to deserialize it.
            visitedNodes_.erase(visitedNodes_.find(element));

            deserialize(key, *data);
        }
    }
}

template<class S, class T>
void XmlDeserializer::deserialize(const std::string& key, std::pair<S, T>& data) {
    // Always deserialize pointer from content?
    if (usePointerContentSerialization_) {
        TemporaryUsePointerContentSerializationChanger usePointerContentSerializationChanger(*this, false);

        Deserializer(*this).deserialize(key, data.first);
        Deserializer(*this).deserialize(key, data.second);

        return;
    }

    TiXmlElement* element = getNextXmlElement(key);

    TemporaryNodeChanger nodeChanger(*this, element);

    // first item
    if (useAttributes_ && isPrimitiveType(typeid(S))) {
        // deserialize primitive type from XML attribute
        deserializeAttributeFromNode("First", XmlSerializationConstants::VALUEATTRIBUTE, data.first);
    }
    else {
        Deserializer(*this).deserialize("First", data.first);
    }

    // second item
    if (useAttributes_ && isPrimitiveType(typeid(T))) {
        // deserialize primitive type from XML attribute
        deserializeAttributeFromNode("Second", XmlSerializationConstants::VALUEATTRIBUTE, data.second);
    }
    else {
        Deserializer(*this).deserialize("Second", data.second);
    }

    addReferenceAddress(element, &data);
}

template<class T>
inline void XmlDeserializer::deserializeSimpleTypes(const std::string& key, T& data) {
    TiXmlElement* element;
    if (useAttributes_)
        element = node_->ToElement();
    else
        element = getNextXmlElement(key);

    const std::string* value;

    if (useAttributes_) {
        value = element->Attribute(key);
        if (!value)
            raise(SerializationNoSuchDataException("No attribute with key '" + key + "' found."));
    }
    else {
        // Is there no value attribute in XML element?
        value = element->Attribute(XmlSerializationConstants::VALUEATTRIBUTE);
        if (!value)
            raise(SerializationFormatException("XML node '" + key + "' does not contain value attribute."));
    }

    convertStringToData(*value, data);

    if (!useAttributes_)
        addReferenceAddress(element, &data);
}

template<>
inline void XmlDeserializer::deserializeSimpleTypes(const std::string& key, std::string& data) {
    TiXmlElement* element;
    // Deserialize from XML attributes?
    if (useAttributes_)
        element = node_->ToElement();
    else
        element = getNextXmlElement(key);

    const std::string* value;
    // Deserialize from XML attributes?
    if (useAttributes_)
        value = element->Attribute(key);
    else
        value = element->Attribute(XmlSerializationConstants::VALUEATTRIBUTE);

    // Is there no matching attribute in XML element?
    if (value) {
        // ATTENTION: No special handling of the given string is needed that is why this block
        //            has to correspond with the not specialized serializeSimpleTypes method.
        data = *value;
    }
    else {
        // Deserialize from XML attributes?
        if (useAttributes_)
            element = getNextXmlElement(key);

        // Is there no text in the XML element?
        if (!element->FirstChild() || element->FirstChild()->Type() != TiXmlElement::TINYXML_TEXT) {
            // Deserialize from XML attributes?
            if (useAttributes_)
                raise(SerializationNoSuchDataException("Neither a '" + key + "' attribute nor text content data exists."));
            else
                raise(SerializationFormatException("XML node '" + key + "' neither has a value attribute nor contains text data."));
        }

        data = element->FirstChild()->ToText()->ValueStr();
    }

    if (!useAttributes_)
        addReferenceAddress(element, &data);
}

template<class T>
inline void XmlDeserializer::deserializeTgtVector(const std::string& key, T& vector, const bool& isColor) {
    TiXmlElement* element = getNextXmlElement(key);

    const std::string* x;
    const std::string* y;
    const std::string* z;
    const std::string* w;
    if (isColor) {
        x = element->Attribute(XmlSerializationConstants::COLORRATTRIBUTE);
        y = element->Attribute(XmlSerializationConstants::COLORGATTRIBUTE);
        z = element->Attribute(XmlSerializationConstants::COLORBATTRIBUTE);
        w = element->Attribute(XmlSerializationConstants::COLORAATTRIBUTE);
    }
    else {
        x = element->Attribute(XmlSerializationConstants::VECTORXATTRIBUTE);
        y = element->Attribute(XmlSerializationConstants::VECTORYATTRIBUTE);
        z = element->Attribute(XmlSerializationConstants::VECTORZATTRIBUTE);
        w = element->Attribute(XmlSerializationConstants::VECTORWATTRIBUTE);
    }

    // Is there no x attribute?
    if (!x)
        raise(SerializationFormatException("Xml node '" + key + "' has no x-attribute."));

    if (isColor) {
        // ATTENTION: The string to int to type conversion is mandatory,
        //            since tgt::col3, tgt::col4 are vector of UCHARs.
        //            Otherwise, strings like "153" are converted to '1' instead of 153.
        int value;
        convertStringToData(*x, value);
        vector[0] = static_cast<typename T::ElemType>(value);
    }
    else
        convertStringToData(*x, vector[0]);

    // Y attribute expected?
    if (vector.size >= 2) {
        // Is there no y attribute?
        if (!y)
            raise(SerializationFormatException("Xml node '" + key + "' has no y-attribute."));

        if (isColor) {
            // ATTENTION: The string to int to type conversion is mandatory,
            //            since tgt::col3, tgt::col4 are vector of UCHARs.
            //            Otherwise, strings like "153" are converted to '1' instead of 153.
            int value;
            convertStringToData(*y, value);
            vector[1] = static_cast<typename T::ElemType>(value);
        }
        else
            convertStringToData(*y, vector[1]);
    }

    // Z attribute expected?
    if (vector.size >= 3) {
        // Is there no z attribute?
        if (!z)
            raise(SerializationFormatException("Xml node '" + key + "' has no z-attribute."));

        if (isColor) {
            // ATTENTION: The string to int to type conversion is mandatory,
            //            since tgt::col3, tgt::col4 are vector of UCHARs.
            //            Otherwise, strings like "153" are converted to '1' instead of 153.
            int value;
            convertStringToData(*z, value);
            vector[2] = static_cast<typename T::ElemType>(value);
        }
        else
            convertStringToData(*z, vector[2]);
    }

    // W attribute expected?
    if ((int)vector.size >= 4) {
        // Is there no w attribute?
        if (!w)
            raise(SerializationFormatException("Xml node '" + key + "' has no w-attribute."));

        if (isColor) {
            // ATTENTION: The string to int to type conversion is mandatory,
            //            since tgt::col3, tgt::col4 are vector of UCHARs.
            //            Otherwise, strings like "153" are converted to '1' instead of 153.
            int value;
            convertStringToData(*w, value);
            vector[3] = static_cast<typename T::ElemType>(value);
        }
        else
            convertStringToData(*w, vector[3]);
    }

    addReferenceAddress(element, &vector);
}

template<class T>
inline void XmlDeserializer::deserializeAttributeFromNode(const std::string& nodeKey, const std::string& key, T& data, const bool& markVisited) {
    TemporaryNodeChanger nodeChanger(*this, getNextXmlElement(nodeKey));

    Deserializer(*this).deserialize(key, data);

    if (!markVisited)
        visitedNodes_.erase(node_);
}

template<class T>
inline void XmlDeserializer::deserializeCollectionConstIteratorItems(const std::string& key, T& collection, const std::string& itemKey) {
    TiXmlElement* element = getNextXmlElement(key);

    TemporaryNodeChanger nodeChanger(*this, element);

    try {
        collection.clear();

        while (true) {
            typename T::value_type item;

            try {
                // Deserialize primitive type from XML attribute?
                if (useAttributes_ && isPrimitiveType(typeid(typename T::value_type))) {
                    deserializeAttributeFromNode(itemKey, XmlSerializationConstants::VALUEATTRIBUTE, item);
                }
                else
                    Deserializer(*this).deserialize(itemKey, item);

                collection.insert(collection.end(), std::move(item));
            }
            catch (SerializationMemoryAllocationException&) {
                // Just ignore items for which no memory was allocatable...
            }
        }
    }
    catch (SerializationNoSuchDataException&) {
        // There are just no more collection items, just ignore the exception...
        removeLastError();
    }

    addReferenceAddress(element, &collection);
}

template<class T>
inline void XmlDeserializer::deserializeCollection(const std::string& key, T& collection, const std::string& itemKey) {
    // Never deserialize pointer from content?
    if (!usePointerContentSerialization_) {
        deserializeCollectionConstIteratorItems(key, collection, itemKey);
        return;
    }

    TiXmlElement* element = getNextXmlElement(key);

    TemporaryNodeChanger nodeChanger(*this, element);

    for (typename T::iterator it = collection.begin(); it != collection.end(); ++it) {
        // Deserialize primitive type from XML attribute?
        if (useAttributes_ && isPrimitiveType(typeid(typename T::value_type))) {
            deserializeAttributeFromNode(itemKey, XmlSerializationConstants::VALUEATTRIBUTE, *it);
        }
        else
            Deserializer(*this).deserialize(itemKey, *it);
    }

    try {
        getNextXmlElement(itemKey);
        // Collection is too small for storing all serialized values. Log it, but do not treat as error.
        //addError(SerializationMemoryAllocationException("Not enough allocated items in collection."));
        LWARNING("Not enough allocated items in collection with key: " << key);
    }
    catch (SerializationNoSuchDataException&) {
        removeLastError();
    }

    addReferenceAddress(element, &collection);
}

template<class T>
inline void XmlDeserializer::deserializeMap(const std::string& key,
                                            T& map,
                                            const std::string& valueKey,
                                            const std::string& keyKey) {
    TiXmlElement* element = getNextXmlElement(key);

    TemporaryNodeChanger nodeChanger(*this, element);

    try {
        if (!usePointerContentSerialization_)
            map.clear();

        while (true) {
            typename T::key_type key;
            typename T::mapped_type value;

            bool memoryAllocationExceptionCatched = false;

            try {
                // Deserialize primitive key type as XML attribute?
                if (useAttributes_ && isPrimitiveType(typeid(typename T::key_type))) {
                    // ATTENTION: We have to set element unvisited in order to deserialize its content later.
                    deserializeAttributeFromNode(valueKey, keyKey, key, false);
                }
                else
                    Deserializer(*this).deserialize(keyKey, key);
            }
            catch (SerializationMemoryAllocationException&) {
                memoryAllocationExceptionCatched = true;
            }

            // Always deserialize pointer from content?
            if (usePointerContentSerialization_) {
                if (map.find(key) == map.end()) {
                    // There is no map item with deserialized key. Log it, but do not treat it as error.
                    //addError(SerializationMemoryAllocationException("No allocated map item with key '" + convertDataToString(key) + "'."));
                    LWARNING("No allocated map item with key '" << convertDataToString(key) << "'.");
                    getNextXmlElement(valueKey);
                    continue;
                }

                // Deserialize primitive mapped type from XML attribute?
                if (useAttributes_
                    && (isPrimitiveType(typeid(typename T::mapped_type)) || isPrimitivePointerType(typeid(typename T::mapped_type))))
                {
                    deserializeAttributeFromNode(valueKey, XmlSerializationConstants::VALUEATTRIBUTE, map.find(key)->second);
                } // Deserialize primitive mapped type from XML node
                else
                    Deserializer(*this).deserialize(valueKey, map.find(key)->second);;
            }
            else {
                try {
                    // Deserialize primitive mapped type as XML attribute?
                    if (useAttributes_ && isPrimitiveType(typeid(typename T::mapped_type)))
                        deserializeAttributeFromNode(valueKey, XmlSerializationConstants::VALUEATTRIBUTE, value);
                    // Deserialize primitive mapped type from XML node
                    else
                        Deserializer(*this).deserialize(valueKey, value);
                }
                catch (SerializationMemoryAllocationException&) {
                    memoryAllocationExceptionCatched = true;
                }
            }

            if (!memoryAllocationExceptionCatched && !usePointerContentSerialization_)
                map[key] = value;
        }
    }
    catch (SerializationNoSuchDataException&) {
        // There are just no more collection items, just ignore the exception...
        removeLastError();
    }

    addReferenceAddress(element, &map);
}

template<class T>
void XmlDeserializer::resolveReference(TiXmlElement* element, T*& data) {
    const std::string* reference = element->Attribute(XmlSerializationConstants::REFERENCEATTRIBUTE);

    // Has XML element no reference attribute?
    if (!reference)
        return;

    // Is referenced data already deserialized?
    if (idAddressMap_.find(*reference) != idAddressMap_.end()) {
        data = static_cast<T*>(idAddressMap_.find(*reference)->second);
    }
    else { // Have to wait for deserialization of referenced data...
        const std::string* type = element->Attribute(XmlSerializationConstants::TYPEATTRIBUTE);

        // Has XML element an type attribute?
        if (type) {
            // Check if later resolving is possible at all?
            delete allocateMemory<T>(data, type);

            // As no exception was raise later resolving is possible, so add unresolved reference...
            addUnresolvedReference(element, data);
        }
        else
            addUnresolvedReference(element, data);
    }
}


template<typename T>
void XmlDeserializer::deserialize(const std::string& key, tgt::TemplateBounds<T>& data) {
    tgt::Vector3<T> llf, urb;
    deserializeTgtVector(key+".llf", llf);
    deserializeTgtVector(key+".urb", urb);
    data = tgt::TemplateBounds<T>(llf, urb);
}

} // namespace

#endif // VRN_XMLDESERIALIZER_H
