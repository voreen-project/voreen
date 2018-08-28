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

#ifndef VRN_METADATACONTAINER_H
#define VRN_METADATACONTAINER_H

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/datastructures/meta/metadatabase.h"

namespace voreen {

/**
 * @c MetaDataContainer stores meta data.
 *
 * @see Serializable
 */
class VRN_CORE_API MetaDataContainer : public Serializable {

    friend class XmlSerializer;
    friend class XmlDeserializer;

public:
    MetaDataContainer() {}
    MetaDataContainer(const MetaDataContainer& mdc);
    MetaDataContainer& operator=(const MetaDataContainer& mdc);
    virtual ~MetaDataContainer();

    /**
     * Adds the given @c key/meta-data pair to the meta data associated with the @c Processor.
     *
     * @attention Already existing meta data with the same key is overwritten.
     *
     * @param key the key
     * @param meta the meta-data to add
     */
    virtual void addMetaData(const std::string& key, MetaDataBase* meta);

    /**
     * Returns if meta data to the given key exists.
     *
     * @param key the key
     * @returns @c true if meta data exists, otherwise @c false.
     */
    virtual bool hasMetaData(const std::string& key) const;

    virtual std::vector<std::string> getKeys() const;

    /**
     * Returns meta data to the given key.
     *
     * @param key the key
     * @returns the meta data to the given key if it exists, otherwise @c 0.
     */
    virtual MetaDataBase* getMetaData(const std::string& key);

    ///@overload
    virtual const MetaDataBase* getMetaData(const std::string& key) const;

    /**
     * Removes meta data with the given key. If no meta data exists, nothing happens
     *
     * @param key the key
     */
    virtual void removeMetaData(const std::string& key);

    /**
     * @brief Changes the key of a MetaData entry.
     * Does nothing if the key doesn't exist.
     */
    virtual void renameMetaData(const std::string& oldKey, const std::string& newKey);

    /**
     * Deletes all associated meta data.
     */
    virtual void clearMetaData();

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

private:
    /**
     * Stored meta data.
     */
    std::map<std::string, MetaDataBase*> metaData_;
};

} // namespace

#endif // VRN_METADATACONTAINER_H
