/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_TEMPLATEFILTERPROPERTIES_H
#define VRN_TEMPLATEFILTERPROPERTIES_H

#include "filterproperties.h"
#include "../volumefiltering/slicereader.h"
#include "../volumefiltering/volumefilter.h"

#include "voreen/core/io/serialization/xmldeserializer.h"
#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/serializer.h"
#include <sstream>

namespace {

// Helper frequently used in FilterProperties subclasses
template<typename T>
void copyPropertyValue(const T& src, T& dst) {
    voreen::XmlSerializer xser("");

    voreen::Serializer ser(xser);
    src.serializeValue(ser);

    std::stringstream s;
    xser.write(s);

    voreen::XmlDeserializer xdeser("");

    xdeser.read(s);

    voreen::Deserializer deser(xdeser);
    dst.deserializeValue(deser);

    dst.invalidate();
}

template<typename Settings>
std::string settingsId(const std::string& id) {
    std::string name = Settings::getVolumeFilterName();
    for(auto& c : name) {
        if(c == ' ') {
            c = '_';
        }
    }
    return name + "_" + id;
}

template<typename T>
void deserializeTemplatePropertyWithValueFallback(voreen::Deserializer& s, std::string id, voreen::TemplateProperty<T>& prop) {
    try {
        s.deserialize(id, prop);
    }
    catch (voreen::SerializationException&) {
        s.removeLastError();
        T val;
        s.deserialize(id, val);
        prop.set(val);
    }
}
}

namespace voreen {

template<typename Settings>
class TemplateFilterProperties : public FilterProperties {

public:
    TemplateFilterProperties();
    virtual ~TemplateFilterProperties() {}

    virtual std::string getVolumeFilterName() const;
    virtual void adjustPropertiesToInput(const SliceReaderMetaData& input, int instanceId);
    virtual VolumeFilter* getVolumeFilter(const SliceReaderMetaData& inputmetadata, int instanceId) const;
    virtual void storeInstance(int instanceId);
    virtual void restoreInstance(int instanceId);
    virtual void removeInstance(int instanceId);
    virtual std::vector<int> getStoredInstances() const;

    void serialize(Serializer& s) const;
    void deserialize(Deserializer& s);

private:
    Settings current_;
    std::map<int, Settings> instanceSettings_;
};

template<typename Settings>
TemplateFilterProperties<Settings>::TemplateFilterProperties()
    : current_()
    , instanceSettings_()
{
    // Add properties to list.
    current_.addProperties(properties_);

    // Store default settings.
    storeInstance(DEFAULT_SETTINGS);
}

template<typename Settings>
std::string TemplateFilterProperties<Settings>::getVolumeFilterName() const {
    return Settings::getVolumeFilterName();
}

template<typename Settings>
void TemplateFilterProperties<Settings>::adjustPropertiesToInput(const SliceReaderMetaData& input, int instanceId) {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }
    instanceSettings_[instanceId].adjustPropertiesToInput(input);
}

template<typename Settings>
VolumeFilter* TemplateFilterProperties<Settings>::getVolumeFilter(const SliceReaderMetaData& inputmetadata, int instanceId) const {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        return nullptr;
    }
    return instanceSettings_.at(instanceId).getVolumeFilter(inputmetadata);
}

template<typename Settings>
void TemplateFilterProperties<Settings>::storeInstance(int instanceId) {
    instanceSettings_[instanceId] = current_;
}
template<typename Settings>
void TemplateFilterProperties<Settings>::restoreInstance(int instanceId) {
    if (instanceSettings_.find(instanceId) == instanceSettings_.end()) {
        instanceSettings_[instanceId] = instanceSettings_[DEFAULT_SETTINGS];
    }

    current_ = instanceSettings_[instanceId];
}
template<typename Settings>
void TemplateFilterProperties<Settings>::removeInstance(int instanceId) {
    instanceSettings_.erase(instanceId);
}

template<typename Settings>
std::vector<int> TemplateFilterProperties<Settings>::getStoredInstances() const {
    std::vector<int> output;
    for(auto& kv : instanceSettings_) {
        if(kv.first != DEFAULT_SETTINGS) {
            output.push_back(kv.first);
        }
    }
    return output;
}

template<typename Settings>
void TemplateFilterProperties<Settings>::serialize(Serializer& s) const {
    s.serialize(getId("instanceSettings"), instanceSettings_);
}
template<typename Settings>
void TemplateFilterProperties<Settings>::deserialize(Deserializer& s) {
    try {
        s.deserialize(getId("instanceSettings"), instanceSettings_);
    }
    catch (SerializationException&) {
        s.removeLastError();
        LERROR("You need to reconfigure " << getVolumeFilterName() << " instances of " << ( properties_[0]->getOwner() ? properties_[0]->getOwner()->getGuiName() : "VolumeFilterList"));
    }
}

}
#endif
