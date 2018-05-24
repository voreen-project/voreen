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

#ifndef VRN_XMLSERIALIZER_H
#define VRN_XMLSERIALIZER_H

#include <string>
#include <utility>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <iostream>

#include "tinyxml/tinyxml.h"

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"

#include "voreen/core/io/serialization/xmlserializerbase.h"
#include "voreen/core/io/serialization/serializationexceptions.h"
#include "voreen/core/io/serialization/xmlserializationconstants.h"
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/serializablefactory.h"

namespace voreen {

/**
 * @c XmlSerializer is responsible for serializing memory data to XML documents.
 *
 * @par
 * The implementation is capable of serializing simple C++ data types, user defined classes
 * which are derived from @c Serializable, and most STL containers containing just the previously
 * mentioned data.
 *
 * @par
 * Furthermore, cycles, joints and polymorphic @c Serializable derivatives are supported which
 * makes the @c XmlSerializer quite flexible. The serialization process is independent
 * of serialization order due to the use of @c key/data pairs, provided that different
 * keys are used at each hierarchy level.
 *
 * @par
 * Finally, the XML document can be written to an arbitrary @c std::ostream. For instance, this
 * allows you to write the XML document into a @c std::stringstream for further processing.
 *
 * @par
 * You have to use @c XmlDeserializer, which is the counterpart to @c XmlSerializer, for
 * deserializing XML documents to memory data.
 *
 * @par
 * Here is a short example of using the @c XmlSerializer:
 * @code
 * int i = 1;
 *
 * std::fstream f;
 * f.open("file.xml", std::ios::out);
 *
 * XmlSerializer s;
 * Serializer(s) serializer;
 * serializer.serialize("i", i);
 * s.write(f);
 *
 * f.close();
 * @endcode
 * For more complex examples and interaction with other classes of the serialization framework
 * see the &quot;serializertest&quot; application in &quot;apps/&quot; directory.
 *
 * @attention All memory data must be serialized before the XML document can be written.
 *            Otherwise, it cannot be ensured that all pointer references are correctly resolved.
 *
 * @note For further information on handling cycles, joints and polymorphism, see:
 *       http://www.parashift.com/c++-faq-lite/serialization.html
 *
 * @see Serializer
 * @see XmlDeserializer
 * @see XmlSerializerBase
 * @see Serializable
 */
class VRN_CORE_API XmlSerializer : public XmlSerializerBase
{
public:
    /**
     * Constructor which initializes the XML document.
     *
     * @par
     * Initialization of the XML document means creating XML declaration and root node
     * as well as adding them to the XML document.
     *
     * @param documentPath Absolute working directory of the document, which is typically the path
     *      to the XML file the document will be written to. This information is not used by the
     *      serializer itself and is therefore not required, but is intended to be accessible
     *      by serializing objects for absolute-to-relative path conversions.
     */
    XmlSerializer(std::string documentPath = "");

    /**
     * Default destructor.
     */
    ~XmlSerializer();

    /**
     * Writes the XML document that contains all already serialized data to the given stream.
     *
     * @attention Keep in mind that all memory data must be serialized before the XML document
     *            can be written. Otherwise, it cannot be ensured that all pointer references
     *            are resolved.
     *
     * @param stream the output stream
     */
    void write(std::ostream& stream);

    /// ----------------------------------------------------------------------------
    /// Generic *Serializer Interface (used by Serializer for dynamic dispatch) ----
    /// ----------------------------------------------------------------------------

    /**
     * Returns the absolute working directory of the document, which is typically the path
     * to the XML file the document will be written to.
     */
    std::string getDocumentPath() const;

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const bool& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const char& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const signed char& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const unsigned char& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const uint16_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const int16_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const uint32_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const int32_t& data);

// There seems to be no uint*_t typedef for long unsigned ints on mac, so we need to provide an implementation for this type.
#ifdef __APPLE__
    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const long unsigned int& data);
#endif

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const uint64_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const int64_t& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const float& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const double& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const long double& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void serialize(const std::string& key, const std::string& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::vec2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::vec3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::vec4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::dvec2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::dvec3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::dvec4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::ivec2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::ivec3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::ivec4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::col3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::col4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::mat2& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::mat3& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::mat4& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::Matrix2d& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::Matrix3d& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    void serialize(const std::string& key, const tgt::Matrix4d& data);

    /**
     * Serializes the given @c Serializable interface realization.
     *
     * @note All user defined classes must realize the @c Serializable
     *       interface to be serializable.
     *
     * @param key the XML node key
     * @param data the @c Serializable realization
     *
     * @throws SerializationAttributeNamingException
     *     if @c serialize method of given @c Serializable raises this exception
     */
    void serialize(const std::string& key, const Serializable& data);

    /**
     * Serializes the given @c key/data pair.
     *
     * @param key the XML node key
     * @param data the data
     */
    template<typename T>
    void serialize(const std::string& key, const tgt::TemplateBounds<T>& data);

    /**
     * Serializes the given pointer reference.
     *
     * @tparam type of referenced data
     *
     * @param key the XML node key
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
     * @param key the XML node key
     * @param data the pair to serialize
     *
     * @throws SerializationAttributeNamingException
     *     if serialization of pair items raises this exception
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
     * @param key the XML node key
     * @param collection the data collection
     * @param itemKey XML node key for each XML child node
     *
     * @throws SerializationAttributeNamingException
     *     if serialization of collection items raises this exception
     */
    template<class T>
    inline void serializeCollection(const std::string& key, const T& collection, const std::string& itemKey = XmlSerializationConstants::ITEMNODE);

    /**
     * Helper function for serializing data maps like STL maps.
     *
     * @note Element order of map items are not guaranteed to remains constant
     *       during serialization and deserialization due to limits of
     *       some STL containers like @c std::map.
     *
     * @tparam T data type of map
     *
     * @param key the XML node key
     * @param map the map
     * @param valueKey XML node key for each value node
     * @param keyKey key for each XML key node or attribute
     *
     * @throws SerializationAttributeNamingException
     *     if serialization of @c key/value pair items raises the exception
     */
    template<class T>
    inline void serializeMap(
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
     * Adds given @c data/node pair to the XML node look up map.
     *
     * @param data data address
     * @param node the XML node
     */
    void addDataNode(const void* data, TiXmlElement* node);

    /**
     * Adds given unresolved reference to the unresolved reference vector.
     *
     * @tparam T type of referenced data
     *
     * @param unresolvedReferenceNode XML node of the unresolved reference
     * @param reference the unresolved reference
     */
    template<class T>
    inline void addUnresolvedReference(TiXmlElement* unresolvedReferenceNode, const T* const& reference);

    /**
     * Checks if the given XML attribute key is not a reserved key and unique.
     *
     * @param key the XML attrubte key
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    void checkAttributeKey(const std::string& key);

    /**
     * Makes the XML attribute key of the given @c element the first XML attribute.
     *
     * @param key key of XML attribute to move
     * @param element the @c element
     */
    void moveAttributeToFront(const std::string& key, TiXmlElement* element);

    /**
     * Resolves unresolved references during @c write method.
     */
    void resolveUnresolvedReferences();

    /**
     * Helper function for serializing @c key/data pairs with simple data types.
     *
     * @tparam T data type
     *
     * @param key the XML node key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    template<class T>
    inline void serializeSimpleTypes(const std::string& key, const T& data);

    /**
     * Helper function for serializing @c key/data pairs with tgt vectors.
     *
     * @tparam T data type of vector
     *
     * @param key the XML node key
     * @param vector the tgt vector
     * @param isColor flag determine if given vector is a color
     */
    template<class T>
    inline void serializeTgtVector(const std::string& key, const T& vector, const bool& isColor = false);

    /**
     * Helper function creates a XML node with the given @c nodeKey and
     * serializes the given @c key/data pair as an attribute.
     *
     * @tparam T data type
     *
     * @param nodeKey the XML node key
     * @param key the XML attribute key
     * @param data the data
     *
     * @throws SerializationAttributeNamingException
     *     if primitive data is serialized as XML attributes and key is reserved or not unique
     */
    template<class T>
    inline void serializeAttributeInNode(const std::string& nodeKey, const std::string& key, const T& data);


    /**
     * Type definition for XML node look up map.
     */
    typedef std::map<void*, TiXmlElement*> DataNodeMapType;

    /**
     * Helper function for base64 encoding.
     */
    std::string base64Encode(const std::vector<unsigned char>& inputBuffer);

    /**
     * Map for looking up XML nodes by pointer address.
     */
    DataNodeMapType dataNodeMap_;

    /**
     * The @c IReferenceContentSerializer interface allows to treat every
     * @c ReferenceContentSerializer template specialization in the same way.
     *
     * @see ReferenceContentSerializer
     * @see UnresolvedReference
     */
    class IReferenceContentSerializer {
    public:

        virtual ~IReferenceContentSerializer() {}

        /**
         * Serializes the content of a reference.
         *
         * @param key the XML node key
         *
         * @returns serialized content
         */
        virtual TiXmlElement* serialize(const std::string& key) const = 0;
    };

    template<class T> class ReferenceContentSerializer;
    template<class T> friend class ReferenceContentSerializer;
    /**
     * @c ReferenceContentSerializer is a concrete implementation for serializing
     * the content of pointer references.
     *
     * @note Serialization of pointer reference content is just necessary when
     *       the pointer data is not serialized, i.e. it is created using @c new.
     *
     * @tparam type of referenced data
     *
     * @see IReferenceContentSerializer
     * @see UnresolvedReference
     */
    template<class T>
    class ReferenceContentSerializer : public IReferenceContentSerializer {
    public:
        /**
         * Constructs a @c ReferenceContentSerializer in accordance to the given parameters.
         *
         * @param serializer Serializer to use for content serialization
         * @param reference the unresolved reference
         */
        ReferenceContentSerializer(XmlSerializer* serializer, const T* const& reference);

        /**
         * @see IReferenceContentSerializer::serialize
         */
        virtual TiXmlElement* serialize(const std::string& key) const;

    private:
        /**
         * Serializer to use for content serialization.
         */
        XmlSerializer* serializer_;

        /**
         * The unresolved reference.
         */
        const T* const reference_;
    };

    /**
     * The @c UnresolvedReference class encapsulates data for resolving unresolved references.
     */
    struct UnresolvedReference {
        /**
         * The XML node of the unresolved pointer reference.
         */
        TiXmlElement* unresolvedReferenceNode;

        /**
         * Type string corresponding with the reference type.
         *
         * @note The string is empty if there is no corresponding type string.
         */
        std::string typeString;

        /**
         * The pointer address of the reference for XML node look up.
         */
        void* reference;

        /**
         * A @c IReferenceContentSerializer instance for serializing reference content
         * when pointed data is not serialized during serialization process.
         *
         * @attention Remember to free assigned memory for @c referenceContentSerializer.
         */
        IReferenceContentSerializer* referenceContentSerializer;
    };

    /**
     * Type definition for unresolved references vector.
     */
    typedef std::vector<UnresolvedReference> UnresolvedReferencesType;

    /**
     * Vector of unresolved references which are resolved during the @c write method.
     */
    UnresolvedReferencesType unresolvedReferences_;

    /**
     * Current id for resolving unresolved references.
     */
    long id_;

    /// Path to the target XML document
    std::string documentPath_;

};
} // namespace voreen

/// ----------------------------------------------------------------------------
/// Implementation -------------------------------------------------------------
/// ----------------------------------------------------------------------------

#include "voreen/core/io/serialization/serializer.h"

namespace voreen {

template<class T>
XmlSerializer::ReferenceContentSerializer<T>::ReferenceContentSerializer(
    XmlSerializer* serializer, const T* const& reference)
    : serializer_(serializer)
    , reference_(reference)
{
}

template<class T>
TiXmlElement* XmlSerializer::ReferenceContentSerializer<T>::serialize(const std::string& key) const {
    TiXmlElement* result = 0;

    // Create temporary XML node...
    TemporaryNodeChanger nodeChange(*serializer_, new TiXmlElement(XmlSerializationConstants::TEMPNODE));

    // Serialize reference content to temporary XML node...
    Serializer(*serializer_).serialize(key, *reference_);
    // ATTENTION: You must not clone or delete temporary node respectively its inner nodes,
    //            otherwise created pointers will point to incorrect address or become invalid.
    //            Remember to delete temporary XML node later.
    result = serializer_->node_->ToElement();

    return result;
}

template<class T>
inline void XmlSerializer::addUnresolvedReference(TiXmlElement* unresolvedReferenceNode,
                                                  const T* const& reference)
{
    // Is null pointer?
    if (!reference)
        // It is refering to nowhere, so we cannot do anything.
        return;

    // Assemble unresolved reference for resolving later on...
    XmlSerializer::UnresolvedReference unresolvedReference;
    unresolvedReference.unresolvedReferenceNode = unresolvedReferenceNode;
    unresolvedReference.typeString = getTypeString(typeid(*reference));
    unresolvedReference.reference = const_cast<T*>(reference);
    unresolvedReference.referenceContentSerializer = new ReferenceContentSerializer<T>(this, reference);

    // As we want ascending ids so far as possible,
    // we have to insert the unresolved references at the beginning.
    unresolvedReferences_.insert(unresolvedReferences_.begin(), unresolvedReference);
}

template<class T>
void XmlSerializer::serialize(const std::string& key, const T* const& data) {
    // Always serialize pointer content?
    if (usePointerContentSerialization_) {
        TemporaryUsePointerContentSerializationChanger usePointerContentSerializationChanger(*this, false);

        Serializer(*this).serialize(key, *data);

        return;
    }

    TiXmlElement* newNode = new TiXmlElement(key);

    // Is data not a null pointer?
    if (data) {
        // Set type attribute if a registered factory supports the data type T...
        std::string typeString = getTypeString(typeid(*data));
        if (!typeString.empty())
            newNode->SetAttribute(XmlSerializationConstants::TYPEATTRIBUTE, typeString);
    }

    node_->LinkEndChild(newNode);

    addUnresolvedReference(newNode, data);
}

template<class S, class T>
void XmlSerializer::serialize(const std::string& key, const std::pair<S, T>& data) {
    TiXmlElement* newNode = new TiXmlElement(key);
    node_->LinkEndChild(newNode);

    TemporaryNodeChanger nodeChanger(*this, newNode);

    // first item
    if (useAttributes_ && isPrimitiveType(typeid(S))) {
        // serialize primitive type as XML attribute
        serializeAttributeInNode("First", XmlSerializationConstants::VALUEATTRIBUTE, data.first);
    }
    else {
        // serialize as XML node
        Serializer(*this).serialize("First", data.first);
    }

    // second item
    if (useAttributes_ && isPrimitiveType(typeid(T))) {
        // serialize primitive type as XML attribute
        serializeAttributeInNode("Second", XmlSerializationConstants::VALUEATTRIBUTE, data.second);
    }
    else {
        // serialize as XML node
        Serializer(*this).serialize("Second", data.second);
    }

    addDataNode(&data, newNode);
}

template<class T>
inline void XmlSerializer::serializeSimpleTypes(const std::string& key, const T& data) {
    // Serialize as XML attribute?
    if (useAttributes_) {
        checkAttributeKey(key);
        node_->ToElement()->SetAttribute(key, convertDataToString(data));
    }
    // Serialize as XML node?
    else {
        TiXmlElement* newNode = new TiXmlElement(key);
        node_->LinkEndChild(newNode);
        newNode->SetAttribute(XmlSerializationConstants::VALUEATTRIBUTE, convertDataToString(data));

        addDataNode(&data, newNode);
    }
}

template<>
inline void XmlSerializer::serializeSimpleTypes(const std::string& key, const std::string& data) {
    // check, if we have to serialize the string as CDATA
    bool requireCDATA = false;
    requireCDATA |= data.find("\n") != std::string::npos;
    requireCDATA |= data.find("\r") != std::string::npos;

    // Serialize as XML attribute wanted and possible?
    if (useAttributes_ && !requireCDATA) {
        checkAttributeKey(key);
        node_->ToElement()->SetAttribute(key, data);
        return;
    }

    // Serialize as XML node...
    TiXmlElement* newNode = new TiXmlElement(key);
    node_->LinkEndChild(newNode);
    if (!requireCDATA) {
        // ATTENTION: No special handling of the given string is needed that is why this block
        //            has to correspond with the not specialized serializeSimpleTypes method.
        newNode->SetAttribute(XmlSerializationConstants::VALUEATTRIBUTE, data);
    }
    else {
        // Serialize string as CDATA to prevent possible conversion or other XML errors...
        TiXmlText* text = new TiXmlText(data);
        text->SetCDATA(true);
        newNode->LinkEndChild(text);
    }

    addDataNode(&data, newNode);
}

template<class T>
inline void XmlSerializer::serializeTgtVector(const std::string& key, const T& vector, const bool& isColor) {
    TiXmlElement* newNode = new TiXmlElement(key);
    node_->LinkEndChild(newNode);

    newNode->SetAttribute(
        isColor ? XmlSerializationConstants::COLORRATTRIBUTE : XmlSerializationConstants::VECTORXATTRIBUTE,
        convertDataToString(isColor ? static_cast<int>(vector[0]) : vector[0]));

    if (vector.size >= 2) {
        newNode->SetAttribute(
            isColor ? XmlSerializationConstants::COLORGATTRIBUTE : XmlSerializationConstants::VECTORYATTRIBUTE,
            convertDataToString(isColor ? static_cast<int>(vector[1]) : vector[1]));
    }

    if (vector.size >= 3) {
        newNode->SetAttribute(
            isColor ? XmlSerializationConstants::COLORBATTRIBUTE : XmlSerializationConstants::VECTORZATTRIBUTE,
            convertDataToString(isColor ? static_cast<int>(vector[2]) : vector[2]));
    }

    if ((int)vector.size >= 4) {
        newNode->SetAttribute(
            isColor ? XmlSerializationConstants::COLORAATTRIBUTE : XmlSerializationConstants::VECTORWATTRIBUTE,
            convertDataToString(isColor ? static_cast<int>(vector[3]) : vector[3]));
    }

    addDataNode(&vector, newNode);
}

template<class T>
inline void XmlSerializer::serializeAttributeInNode(const std::string& nodeKey, const std::string& key, const T& data) {
    TiXmlNode* newNode = new TiXmlElement(nodeKey);
    node_->LinkEndChild(newNode);

    TemporaryNodeChanger nodeChanger(*this, newNode);

    Serializer(*this).serialize(key, data);
}

template<class T>
inline void XmlSerializer::serializeCollection(const std::string& key, const T& collection, const std::string& itemKey) {
    TiXmlElement* newNode = new TiXmlElement(key);
    node_->LinkEndChild(newNode);

    TemporaryNodeChanger nodeChanger(*this, newNode);

    for (typename T::const_iterator it = collection.begin(); it != collection.end(); ++it) {
        // Serialize primitive type as XML attribute?
        if (useAttributes_ && isPrimitiveType(typeid(typename T::value_type))) {
            serializeAttributeInNode(itemKey, XmlSerializationConstants::VALUEATTRIBUTE, *it);
        }
        // Serialize as XML node...
        else {
            Serializer(*this).serialize(itemKey, *it);
        }
    }

    addDataNode(&collection, newNode);
}

template<class T>
inline void XmlSerializer::serializeMap(const std::string& key,
                                        const T& map,
                                        const std::string& valueKey,
                                        const std::string& keyKey) {
    TiXmlElement* newNode = new TiXmlElement(key);
    node_->LinkEndChild(newNode);

    TemporaryNodeChanger nodeChanger(*this, newNode);

    for (typename T::const_iterator it = map.begin(); it != map.end(); ++it)
    {
        // Serialize primitive key type as XML node?
        if (!useAttributes_ || !isPrimitiveType(typeid(typename T::key_type)))
            Serializer(*this).serialize(keyKey, it->first);

        // Serialize primitive mapped type as XML attribute?
        if (useAttributes_ && isPrimitiveType(typeid(typename T::mapped_type))) {
            serializeAttributeInNode(valueKey, XmlSerializationConstants::VALUEATTRIBUTE, it->second);
        }
        // Serialize as XML node...
        else {
            if (usePointerContentSerialization_ && useAttributes_ && isPrimitivePointerType(typeid(typename T::mapped_type))) {
                serializeAttributeInNode(valueKey, XmlSerializationConstants::VALUEATTRIBUTE, it->second);
            }
            else
                Serializer(*this).serialize(valueKey, it->second);
        }

        // Serialize primitive key type as XML attribute?
        if (useAttributes_ && isPrimitiveType(typeid(typename T::key_type))) {
            TemporaryNodeChanger innerNodeChanger(*this, node_->LastChild());

            Serializer(*this).serialize(keyKey, it->first);

            moveAttributeToFront(keyKey, node_->ToElement());
        }
    }

    addDataNode(&map, newNode);
}

template<typename T>
void XmlSerializer::serialize(const std::string& key, const tgt::TemplateBounds<T>& data) {
    tgt::Vector3<T> llf = data.getLLF();
    tgt::Vector3<T> urb = data.getURB();
    serializeTgtVector(key+".llf", llf);
    serializeTgtVector(key+".urb", urb);
}
} // namespace

#endif // VRN_XMLSERIALIZER_H
