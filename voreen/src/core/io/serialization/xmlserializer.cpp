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

#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/animation/animation.h"

namespace voreen {

const std::string XmlSerializer::loggerCat_ = "voreen.XmlSerializer";

XmlSerializer::XmlSerializer(std::string documentPath)
    : XmlSerializerBase()
    , id_(0)
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

    // Insert XML declaration...
    document_.LinkEndChild(new TiXmlDeclaration(
        XmlSerializationConstants::XMLVERSION,
        XmlSerializationConstants::XMLENCODING,
        XmlSerializationConstants::XMLSTANDALONE));

    // Create XML root node...
    TiXmlElement* root = new TiXmlElement(XmlSerializationConstants::ROOTNODE);
    root->SetAttribute(
        XmlSerializationConstants::VERSIONATTRIBUTE,
        XmlSerializationConstants::VERSION);
    document_.LinkEndChild(root);

    node_ = root;
}

XmlSerializer::~XmlSerializer() {
    for (UnresolvedReferencesType::iterator it = unresolvedReferences_.begin();
        it != unresolvedReferences_.end(); ++it) {
        delete it->referenceContentSerializer;
    }
}

std::string XmlSerializer::getDocumentPath() const {
    return documentPath_;
}

void XmlSerializer::addDataNode(const void* data, TiXmlElement* node) {
    dataNodeMap_[const_cast<void*>(data)] = node;
}

void XmlSerializer::checkAttributeKey(const std::string& key) {
    // Is key a reserved XML attribute name?
    if (key == XmlSerializationConstants::IDATTRIBUTE
        || key == XmlSerializationConstants::REFERENCEATTRIBUTE
        || key == XmlSerializationConstants::TYPEATTRIBUTE)
    {
        raise(SerializationAttributeNamingException("'" + key + "' is a reserved XML attribute name."));
    }

    // Does key already exists?
    if (node_->ToElement()->Attribute(key))
        raise(SerializationAttributeNamingException("XML attribute with name '" + key + "' already exists."));
}

void XmlSerializer::moveAttributeToFront(const std::string& key, TiXmlElement* element) {
    // Is there no XML attribute with given key?
    if (!element->Attribute(key))
        return;

    // Store value of XML attribute with given key for insertion later on...
    std::string value = *element->Attribute(key);

    // Remove attribute...
    element->RemoveAttribute(key);

    typedef std::vector<std::pair<std::string, std::string> > AttributeListType;
    AttributeListType attributes;

    // Store all existing XML attributes...
    for (TiXmlAttribute* attribute = element->FirstAttribute(); attribute != 0; attribute = attribute->Next())
        attributes.push_back(std::make_pair(attribute->NameTStr(), attribute->ValueStr()));

    // Remove all existing XML attributes...
    for (AttributeListType::iterator it = attributes.begin(); it != attributes.end(); ++it)
        element->RemoveAttribute(it->first);

    // First, insert attribute which should be at front...
    element->SetAttribute(key, value);

    // Restore other attributes...
    for (AttributeListType::iterator it = attributes.begin(); it != attributes.end(); ++it)
        element->SetAttribute(it->first, it->second);
}

void XmlSerializer::resolveUnresolvedReferences() {
    typedef std::pair<TiXmlNode*, TiXmlNode*> ReplaceListEntryType;
    typedef std::vector<ReplaceListEntryType> ReplaceListType;
    // ATTENTION: You have to serialize all reference contents before replacing the XML nodes.
    //            Otherwise, you will invalid some pointers.
    ReplaceListType replaceList;

    // ATTENTION: Simple iteration over all unresolved references is not possible due to
    //            the possibility of adding further unresolved references during resolving process.
    while (!unresolvedReferences_.empty()) {
        UnresolvedReference unresolvedReference = unresolvedReferences_.back();
        unresolvedReferences_.pop_back();

        // Is referenced data already serialized?
        if (dataNodeMap_.find(unresolvedReference.reference) != dataNodeMap_.end()) {
            TiXmlElement* referencedDataNode = dataNodeMap_.find(unresolvedReference.reference)->second;

            // Is id attribute not set on referenced data node?
            if (!referencedDataNode->Attribute(XmlSerializationConstants::IDATTRIBUTE)) {
                // Add id attribute to the referenced data node...
                referencedDataNode->SetAttribute(
                    XmlSerializationConstants::IDATTRIBUTE,
                    XmlSerializationConstants::IDPREFIX + convertDataToString(id_++));
            }

            // Add reference attribute to the unresolved reference node...
            unresolvedReference.unresolvedReferenceNode->SetAttribute(
                XmlSerializationConstants::REFERENCEATTRIBUTE,
                *referencedDataNode->Attribute(XmlSerializationConstants::IDATTRIBUTE));
        }
        else {
            const std::string key = unresolvedReference.unresolvedReferenceNode->ValueStr();

            TiXmlNode* contentNode = unresolvedReference.referenceContentSerializer->serialize(key);

            typedef std::vector<std::pair<std::string, std::string> > AttributeListType;
            AttributeListType attributes;

            // Store all newly created attributes in contentNode...
            for (TiXmlAttribute* attribute = contentNode->FirstChildElement()->FirstAttribute()
                ; attribute != 0
                ; attribute = attribute->Next())
            {
                attributes.push_back(std::make_pair(attribute->NameTStr(), attribute->ValueStr()));
            }

            // Remove all attributes from contentNode...
            for (AttributeListType::iterator it = attributes.begin(); it != attributes.end(); ++it)
                contentNode->FirstChildElement()->RemoveAttribute(it->first);


            // Transfer all existing attributes to the newly created contentNode...
            for (TiXmlAttribute* attribute = unresolvedReference.unresolvedReferenceNode->FirstAttribute()
                ; attribute != 0
                ; attribute = attribute->Next())
            {
                contentNode->FirstChildElement()->SetAttribute(
                    attribute->NameTStr(),
                    attribute->ValueStr());
            }

            // Restore all saved attributes from contentNode...
            for (AttributeListType::iterator it = attributes.begin(); it != attributes.end(); ++it)
                contentNode->FirstChildElement()->SetAttribute(it->first, it->second);

            // When type string corresponding with reference type exists,
            // set type attribute appropriatly...
            if (!unresolvedReference.typeString.empty()) {
                contentNode->FirstChildElement()->SetAttribute(
                    XmlSerializationConstants::TYPEATTRIBUTE,
                    unresolvedReference.typeString);
            }

            replaceList.push_back(std::make_pair(unresolvedReference.unresolvedReferenceNode, contentNode));
        }

        // ATTENTION: Do not forget to free memory of referenceContentSerializer
        delete unresolvedReference.referenceContentSerializer;
    }

    // Replace all unresolved reference XML nodes by its resolved counterparts...
    // ATTENTION: You have to replace XML nodes bottom up, otherwise serialized
    //            data will be lost and you may try to replace an already deleted XML node.
    // ATTENTION: Remember to delete temporary XML node after replacement!
    while (!replaceList.empty()) {
        ReplaceListEntryType replaceEntry = replaceList.back();
        replaceList.pop_back();

        TiXmlNode* unresolvedReferenceParentNode = replaceEntry.first->Parent();
        // ATTENTION: The resolved reference XML node is the first child of a temporary XML node!
        TiXmlNode* resolvedReferenceNode = replaceEntry.second->FirstChildElement();

        // ATTENTION: ReplaceChild-call copies the resolved reference node, that is why
        //            we are responsible for deleting the resolved reference node,
        //            respectively the temporary XML node.
        unresolvedReferenceParentNode->ReplaceChild(replaceEntry.first, *resolvedReferenceNode);

        // Delete temporary XML node, especially this deletes the inner resolved reference node...
        delete replaceEntry.second;
    }
}

void XmlSerializer::serialize(const std::string& key, const bool& data) {
    serializeSimpleTypes(key, (data ? "true" : "false"));
}

void XmlSerializer::serialize(const std::string& key, const char& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const signed char& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const unsigned char& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const uint16_t& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const int16_t& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const uint32_t& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const int32_t& data) {
    serializeSimpleTypes(key, data);
}

#ifdef __APPLE__
void XmlSerializer::serialize(const std::string& key, const long unsigned int& data) {
    serializeSimpleTypes(key, data);
}
#endif

void XmlSerializer::serialize(const std::string& key, const uint64_t& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const int64_t& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const float& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const double& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const long double& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const std::string& data) {
    serializeSimpleTypes(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::vec2& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::vec3& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::vec4& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::dvec2& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::dvec3& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::dvec4& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::ivec2& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::ivec3& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::ivec4& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::svec2& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::svec3& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::svec4& data) {
    serializeTgtVector(key, data);
}

void XmlSerializer::serialize(const std::string& key, const tgt::col3& data) {
    serializeTgtVector(key, data, true);
}

void XmlSerializer::serialize(const std::string& key, const tgt::col4& data) {
    serializeTgtVector(key, data, true);
}

void XmlSerializer::serialize(const std::string& key, const tgt::mat2& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
}

void XmlSerializer::serialize(const std::string& key, const tgt::mat3& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
}

void XmlSerializer::serialize(const std::string& key, const tgt::mat4& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
    serializeTgtVector(key+".row3", data[3]);
}

void XmlSerializer::serialize(const std::string& key, const tgt::Matrix2d& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
}

void XmlSerializer::serialize(const std::string& key, const tgt::Matrix3d& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
}

void XmlSerializer::serialize(const std::string& key, const tgt::Matrix4d& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
    serializeTgtVector(key+".row3", data[3]);
}

void XmlSerializer::serialize(const std::string& key, const Serializable& data) {
    TiXmlElement* newNode = new TiXmlElement(key);
    node_->LinkEndChild(newNode);

    TemporaryNodeChanger nodeChanger(*this, newNode);

    Serializer serializer(*this);
    data.serialize(serializer);

    addDataNode(&data, newNode);
}

void XmlSerializer::write(std::ostream& stream) {
    resolveUnresolvedReferences();

    TiXmlPrinter printer;
    document_.Accept(&printer);
    stream << printer.Str();
}


} // namespace
