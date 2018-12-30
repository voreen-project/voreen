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

#include "voreen/core/datastructures/meta/metadatacontainer.h"

namespace voreen {

MetaDataContainer::MetaDataContainer(const MetaDataContainer& mdc) {
    for (std::map<std::string, MetaDataBase*>::const_iterator it = mdc.metaData_.begin(); it != mdc.metaData_.end(); ++it) {
        const MetaDataBase* md = it->second;
        addMetaData(it->first, md->clone());
    }
}

MetaDataContainer::~MetaDataContainer() {
    clearMetaData();
}

void MetaDataContainer::addMetaData(const std::string& key, MetaDataBase* meta) {
    removeMetaData(key);
    metaData_[key] = meta;
}

bool MetaDataContainer::hasMetaData(const std::string& key) const {
    return metaData_.find(key) != metaData_.end();
}

MetaDataBase* MetaDataContainer::getMetaData(const std::string& key) {
    if (hasMetaData(key))
        return metaData_[key];
    else
        return 0;
}

const MetaDataBase* MetaDataContainer::getMetaData(const std::string& key) const {
    if (hasMetaData(key))
        return metaData_.find(key)->second;//metaData_.at(key); the function "at" does not exists under Visual Studio 2008
    else
        return 0;
}

void MetaDataContainer::removeMetaData(const std::string& key) {
    if (!hasMetaData(key))
        return;

    delete metaData_[key];
    metaData_.erase(metaData_.find(key));
}

void MetaDataContainer::renameMetaData(const std::string& oldKey, const std::string& newKey) {
    if (!hasMetaData(oldKey) || hasMetaData((newKey)))
        return;

    metaData_[newKey] = metaData_[oldKey];
    metaData_.erase(metaData_.find(oldKey));
}

void MetaDataContainer::clearMetaData() {
    for (std::map<std::string, MetaDataBase*>::iterator it = metaData_.begin(); it != metaData_.end(); ++it)
        delete it->second;

    metaData_.clear();
}

void MetaDataContainer::serialize(Serializer& s) const {
    if (metaData_.empty())
        return;

    s.serialize("MetaData", metaData_, "MetaItem", "name");
}

void MetaDataContainer::deserialize(Deserializer& s) {
    clearMetaData();

    try {
        s.deserialize("MetaData", metaData_, "MetaItem", "name");
    }
    catch (SerializationNoSuchDataException&) {
        // There were no meta data during serialization, so we can just ignore the exception...
        s.removeLastError();
    }
}

std::vector<std::string> MetaDataContainer::getKeys() const {
    std::vector<std::string> keys;
    for (std::map<std::string, MetaDataBase*>::const_iterator it = metaData_.begin(); it != metaData_.end(); ++it) {
        keys.push_back(it->first);
    }
    return keys;
}

MetaDataContainer& MetaDataContainer::operator=(const MetaDataContainer& mdc) {
    clearMetaData();
    for (std::map<std::string, MetaDataBase*>::const_iterator it = mdc.metaData_.begin(); it != mdc.metaData_.end(); ++it) {
        const MetaDataBase* md = it->second;
        addMetaData(it->first, md->clone());
    }

    return *this;
}

} // namespace
