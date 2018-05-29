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

#include "voreen/core/properties/propertyowner.h"
#include "voreen/core/properties/property.h"

#include "tgt/logmanager.h"

#include <fstream>

namespace voreen {

void PropertyOwnerObserver::preparePropertyRemoval(Property*) {}
void PropertyOwnerObserver::propertiesChanged(const PropertyOwner* owner) {}
void PropertyOwnerObserver::propertyValueHasBeenModified(Property* property) {}

const std::string PropertyOwner::loggerCat_("voreen.PropertyOwner");

PropertyOwner::PropertyOwner()
    : VoreenSerializableObject()
    , invalidationLevel_(50)
{}

PropertyOwner::PropertyOwner(const std::string& id, const std::string& guiName)
    : VoreenSerializableObject(id, guiName)
    , invalidationLevel_(50)
{}

PropertyOwner::~PropertyOwner() {
}

void PropertyOwner::addProperty(Property* prop) {
    tgtAssert(prop, "Null pointer passed");
    if (prop->getID().empty()) {
        if (tgt::Singleton<tgt::LogManager>::isInited())
            LERROR(getID() << ": Property with empty id (type: " << prop->getClassName() << ")");
        return;
    }
    if (getProperty(prop->getID())) {
        if (tgt::Singleton<tgt::LogManager>::isInited())
            LERROR(getID() << ": Duplicate property id '" << prop->getID() << "'");
        return;
    }
    properties_.push_back(prop);
    prop->setOwner(this);
}

void PropertyOwner::addProperty(Property& prop) {
    addProperty(&prop);
}

void PropertyOwner::removeProperty(Property* prop) {
    tgtAssert(prop, "Null pointer passed");
    if (!getProperty(prop->getID())) {
        LERROR(getID() << ": Property '" << prop->getID() << "' cannot be removed, it does not exist");
    }
    // inform the observers to prepare property removal
    // thus all links can be removed in the processornetwork
    std::vector<PropertyOwnerObserver*> observers = getObservers();
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->preparePropertyRemoval(prop);

    prop->setOwner(0);
    properties_.erase(find(properties_.begin(), properties_.end(), prop));
}

void PropertyOwner::removeProperty(Property& prop) {
    removeProperty(&prop);
}

void PropertyOwner::notifyPropertiesChanged() const {
    std::vector<PropertyOwnerObserver*> propObservers = /*Observable<PropertyOwnerObserver>::*/getObservers();
    for (size_t i = 0; i < propObservers.size(); ++i)
        propObservers[i]->propertiesChanged(this);
}

void PropertyOwner::notifyPropertyValueHasBeenModified(Property* prop) const {
    std::vector<PropertyOwnerObserver*> propObservers = /*Observable<PropertyOwnerObserver>::*/getObservers();
    for (size_t i = 0; i < propObservers.size(); ++i)
        propObservers[i]->propertyValueHasBeenModified(prop);
}

void PropertyOwner::resetAllProperties(){
    for (size_t i = 0; i < properties_.size(); ++i)
        properties_[i]->reset();
}

const std::vector<Property*>& PropertyOwner::getProperties() const {
    return properties_;
}

Property* PropertyOwner::getProperty(const std::string& id) const {
    for (size_t i = 0; i < properties_.size(); ++i) {
        if (properties_[i]->getID() == id)
            return properties_[i];
    }
    return 0;
}

Property* PropertyOwner::getPropertyByName(const std::string& guiName) const {
    for (size_t i = 0; i < properties_.size(); ++i) {
        if (properties_[i]->getGuiName() == guiName)
            return properties_[i];
    }
    return 0;
}

void PropertyOwner::setPropertyGroupGuiName(const std::string& groupID, const std::string& name) {
    for(size_t i = 0; i < properties_.size() ; ++i)
        if(properties_[i]->getGroupID().compare(groupID) == 0)
            properties_[i]->setGroupName(name);
}

std::string PropertyOwner::getPropertyGroupGuiName(const std::string& id) const {
    for(size_t i = 0; i < properties_.size(); ++i)
        if (properties_[i]->getGroupID().compare(id) == 0)
            return properties_[i]->getGroupName();

    return "";
}

void PropertyOwner::setPropertyGroupVisible(const std::string& id, bool visible) {
    for(size_t i = 0; i < properties_.size() ; ++i)
        if (properties_[i]->getGroupID().compare(id) == 0)
            properties_[i]->setVisibleFlag(visible);
}

int PropertyOwner::getInvalidationLevel() const {
    return invalidationLevel_;
}

void PropertyOwner::invalidate(int inv) {
    invalidationLevel_ = std::max(invalidationLevel_, inv);
}

bool PropertyOwner::isValid() const {
    return (getInvalidationLevel() == 0);
}

void PropertyOwner::setValid() {
    invalidationLevel_ = 0;
}

void PropertyOwner::serialize(Serializer& s) const {
    // create temporary property map for serialization
    std::map<std::string, Property*> propertyMap;
    for (std::vector<Property*>::const_iterator it = properties_.begin(); it != properties_.end(); ++it)
        propertyMap[(*it)->getID()] = *it;

    // serialize properties
    const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    s.setUsePointerContentSerialization(true);
    try {
        s.serialize("Properties", propertyMap, "Property", "mapKey");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }
    s.setUsePointerContentSerialization(usePointerContentSerialization);
}

void PropertyOwner::deserialize(Deserializer& s) {
    // create temporary property map for deserialization
    std::map<std::string, Property*> propertyMap;
    for (std::vector<Property*>::const_iterator it = properties_.begin(); it != properties_.end(); ++it)
        propertyMap[(*it)->getID()] = *it;

    // deserialize properties
    const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    s.setUsePointerContentSerialization(true);
    try {
        s.deserialize("Properties", propertyMap, "Property", "mapKey");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }
    s.setUsePointerContentSerialization(usePointerContentSerialization);
}

bool serializeSettings(const PropertyOwner* po, const std::string& filename) {
    std::ofstream stream(filename.c_str(), std::ios_base::out);
    if (stream.fail()) {
        LWARNINGC("VoreenSettings", "Unable to open file " << filename << " for writing.");
        return false;
    }

    bool success = true;
    try {
        XmlSerializer xmlSerializer;
        Serializer serializer(xmlSerializer);
        po->serialize(serializer);
        xmlSerializer.write(stream);
        if (stream.bad()) {
            LWARNINGC("VoreenSettings", "Unable to write to file: " << filename);
            success = false;
        }
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNINGC("VoreenSettings", "SerializationException: " << e.what());
        stream.close();
        success = false;
    }
    return success;
}

bool deserializeSettings(PropertyOwner* po, const std::string& filename) {
    std::ifstream stream;
    stream.open(filename.c_str(), std::ios_base::in);
    if(stream.fail()) {
        stream.close();
        return false;
    }
    else {
        XmlDeserializer xmlDeserializer;
        Deserializer deserializer(xmlDeserializer);
        try {
            xmlDeserializer.read(stream);
            po->deserialize(deserializer);
            stream.close();
        }
        catch (SerializationNoSuchDataException) {
            // no data present => ignore
            xmlDeserializer.removeLastError();
            return false;
        }
        catch (SerializationException &e) {
            LWARNINGC("VoreenSettingsDialog", std::string("Deserialization failed: ") + e.what());
            return false;
        }
    }
    return true;
}

} // namespace voreen
