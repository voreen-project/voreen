/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "voreen/core/io/serialization/xmldeserializer.h"
#include "voreen/core/io/serialization/deserializer.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/animation/animation.h"

namespace voreen {

const std::string XmlDeserializer::loggerCat_ = "voreen.XmlDeserializer";

XmlDeserializer::XmlDeserializer(std::string documentPath)
    : XmlSerializerBase()
    , documentPath_(documentPath)
{
    // register application (as proxy for modules)
    if (VoreenApplication::app()) {
        registerFactory(VoreenApplication::app());
        registerFactories(VoreenApplication::app()->getSerializerFactories());
    }
    else {
        LWARNING("VoreenApplication not instantiated");
    }
}

XmlDeserializer::~XmlDeserializer() {
    for (UnresolvedReferenceMapType::iterator mapIt = unresolvedReferenceMap_.begin();
        mapIt != unresolvedReferenceMap_.end(); ++mapIt)
    {
        for (ReferenceResolverListType::iterator listIt = mapIt->second.begin();
            listIt != mapIt->second.end(); ++listIt)
        {
            delete *listIt;
        }
    }
}

std::string XmlDeserializer::getDocumentPath() const {
    return documentPath_;
}

void XmlDeserializer::deserialize(const std::string& key, bool& data) {
    std::string boolValue;
    deserializeSimpleTypes(key, boolValue);

    std::transform(boolValue.begin(), boolValue.end(), boolValue.begin(), tolower);

    if (boolValue == "true" || boolValue == "1")
        data = true;
    else if (boolValue == "false" || boolValue == "0")
        data = false;
    else
        raise(SerializationFormatException("XML node with key '" + key + "' contains unknown bool value."));
}

void XmlDeserializer::deserialize(const std::string& key, char& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, signed char& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, unsigned char& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, uint16_t& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, int16_t& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, uint32_t& data) {
    deserializeSimpleTypes(key, data);
}

#ifdef __APPLE__
void XmlDeserializer::deserialize(const std::string& key, long unsigned int& data) {
    deserializeSimpleTypes(key, data);
}
#endif

void XmlDeserializer::deserialize(const std::string& key, int32_t& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, uint64_t& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, int64_t& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, float& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, double& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, long double& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, std::string& data) {
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::vec2& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::vec3& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::vec4& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::dvec2& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::dvec3& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::dvec4& data) {
    deserializeTgtVector(key, data);
}


void XmlDeserializer::deserialize(const std::string& key, tgt::ivec2& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::ivec3& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::ivec4& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::svec2& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::svec3& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::svec4& data) {
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::col3& data) {
    deserializeTgtVector(key, data, true);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::col4& data) {
    deserializeTgtVector(key, data, true);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::mat2& data) {
    tgt::vec2 row0, row1;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    data = tgt::mat2(row0, row1);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::mat3& data) {
    tgt::vec3 row0, row1, row2;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    data = tgt::mat3(row0, row1, row2);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::mat4& data) {
    tgt::vec4 row0, row1, row2, row3;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    deserializeTgtVector(key+".row3", row3);
    data = tgt::mat4(row0, row1, row2, row3);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::Matrix2d& data) {
    tgt::dvec2 row0, row1;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    data = tgt::Matrix2d(row0, row1);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::Matrix3d& data) {
    tgt::dvec3 row0, row1, row2;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    data = tgt::Matrix3d(row0, row1, row2);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::Matrix4d& data) {
    tgt::dvec4 row0, row1, row2, row3;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    deserializeTgtVector(key+".row3", row3);
    data = tgt::Matrix4d(row0, row1, row2, row3);
}

void XmlDeserializer::deserialize(const std::string& key, Serializable& data) {
    TiXmlElement* element = getNextXmlElement(key);

    TemporaryNodeChanger nodeChanger(*this, element);

    Deserializer deserializer(*this);
    data.deserialize(deserializer);

    addReferenceAddress(element, &data);
}

TiXmlElement* XmlDeserializer::getNextXmlElement(const std::string& key) {
    TiXmlElement* element = node_->FirstChildElement(key);
    while (element) {
        // Was node not visited before?
        if (visitedNodes_.find(element) == visitedNodes_.end())
        {
            visitedNodes_.insert(element);
            return element;
        }

        element = element->NextSiblingElement(key);
    }

    raise(SerializationNoSuchDataException("No further XML node with key '" + key + "' found."));
    return 0;
}

void XmlDeserializer::findUnresolvableReferences(
    TiXmlElement* node,
    ReferenceIdSetType& references,
    ReferenceIdSetType& resolvableReferences) const
{
    const std::string* id = node->Attribute(XmlSerializationConstants::IDATTRIBUTE);
    const std::string* ref = node->Attribute(XmlSerializationConstants::REFERENCEATTRIBUTE);

    if (id)
        resolvableReferences.insert(*id);

    if (ref)
        references.insert(*ref);

    for (TiXmlElement* child = node->FirstChildElement(); child != 0; child = child->NextSiblingElement())
        findUnresolvableReferences(child, references, resolvableReferences);
}

XmlDeserializer::ReferenceIdListType XmlDeserializer::findUnresolvableReferences() const
{
    ReferenceIdListType unresolvableReferences;

    ReferenceIdSetType references;
    ReferenceIdSetType resolvableReferences;
    findUnresolvableReferences(node_->ToElement(), references, resolvableReferences);

    for (ReferenceIdSetType::const_iterator it = references.begin(); it != references.end(); ++it)
        if (resolvableReferences.find(*it) == resolvableReferences.end())
            unresolvableReferences.push_back(*it);

    return unresolvableReferences;
}

void XmlDeserializer::read(std::istream& stream, XmlProcessor* xmlProcessor) {
    // Read input stream...
    std::stringbuf buffer;
    do
    {
        // Use 0 character instead of '\n' to minimize the number of get-calls...
        stream.get(buffer, 0);
    } while (stream.good() && !stream.eof()
        && (buffer.sputc(stream.get()) != std::stringbuf::traits_type::eof()));

    // Parse input...
    document_.Parse(buffer.str().c_str());

    TiXmlElement* root = document_.RootElement();

    // Is there no root element?
    if (!root)
        raise(SerializationFormatException(std::string("No root node found.")));

    // Has root node incorrect name?
    if (root->ValueStr() != XmlSerializationConstants::ROOTNODE) {
        raise(SerializationFormatException("XML root node name is '" + root->ValueStr()
            + "' instead of '" + XmlSerializationConstants::ROOTNODE + "'."));
    }

    const std::string* version = root->Attribute(XmlSerializationConstants::VERSIONATTRIBUTE);
    // Is serialization version not set?
    if (!version)
        raise(SerializationFormatException("XML root node has no version attribute."));
    // Does XmlSerializer and XmlDeserializer version not match the XML document version?
    if (*version != XmlSerializationConstants::VERSION) {
        raise(SerializationVersionMismatchException("XML document has version " + *version
            + " instead of " + XmlSerializationConstants::VERSION + "."));
    }

    node_ = root;

    // Apply preprocessor, if a XML preprocessor is given...
    if (xmlProcessor)
        xmlProcessor->process(document_);

    std::vector<std::string> unresolvableReferences = findUnresolvableReferences();
    if (!unresolvableReferences.empty()) {
        std::stringstream idStream;
        for (ReferenceIdListType::iterator it = unresolvableReferences.begin(); it != unresolvableReferences.end(); ++it)
            idStream << (idStream.str().empty() ? "" : ", ") << "'" << *it << "'";

        raise(SerializationReferenceResolvingException(
            "XML document contains the following unresolvable references: " + idStream.str() + "."));
    }
}

} // namespace
