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

#ifndef VRN_VVDFORMAT_H
#define VRN_VVDFORMAT_H

#include "voreen/core/io/serialization/serializable.h"

#include "voreen/core/datastructures/volume/volume.h"

namespace voreen {


/**
 * This is the base volume format of volumes in Voreen.
 * Each vvd-file can point to a raw-file or an octree-folder.
 */

// See below
class VVDSliceWriter;

/** Base class to combine raw files and octree files */
class VRN_CORE_API VvdDataObjectBase : public Serializable {
    /**
     * Unfortunately VVD Serialization is not made to save Volumes
     * that are just created on the fly as is done in OctreeSave
     * via VVDSliceWriter so that we need to access and overwrite
     * both dimensions_ and format_ in the VVDObject after creating
     * it with another volume. :(
     */
    friend class VVDSliceWriter;
public:
    VvdDataObjectBase() {}
    VvdDataObjectBase(const VolumeBase* volume);

    tgt::ivec3 getDimensions() { return dimensions_; }
    std::string getFormat() { return format_; }

    virtual void serialize(Serializer& s) const = 0;
    virtual void deserialize(Deserializer& s) = 0;
protected:
    tgt::ivec3 dimensions_;
    std::string format_;
};

class VRN_CORE_API VvdRawDataObject : public VvdDataObjectBase {
public:
    VvdRawDataObject() {}
    /** If no absRawOa */
    VvdRawDataObject(const VolumeBase* volume, std::string absRawPath);

    std::string getRawPath() { return absRawPath_; }

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
public:
    std::string absRawPath_;
};

/** A RAW data object can serialize the volume data in the same file.
 *  This is used in the Flowreen module.
 */
class VRN_CORE_API VvdRawDataNoExternalFileObject : public VvdDataObjectBase {
public:
    VvdRawDataNoExternalFileObject() {}
    VvdRawDataNoExternalFileObject(const VolumeBase* volume);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
public:
    std::vector<char> rawData_;
};



class VRN_CORE_API VvdOctreeDataObject : public VvdDataObjectBase {
public:
    VvdOctreeDataObject() {}
    VvdOctreeDataObject(const VolumeBase* volume, std::string filename);

    std::string getFilename() { return filename_; }

    virtual void serialize(Serializer& s) const {};
    virtual void deserialize(Deserializer& s) {};
public:
    std::string filename_;
};


/**
 * This class is used to serialize or deserialize vvd-files.
 */
class VRN_CORE_API VvdObject : public Serializable {
public:
    /**
     * Constructor
     * @param vb VolumeBase used to copy MetaData and DerivedData
     * @param absRAWPath absolut path to the .raw file (empty in octree case)
     * @param useExternalRAWFile if false, the RAW data will be serialized in the VVD file itself
     */
    VvdObject(const VolumeBase* vb, std::string absRAWPAth = "", bool useExternalRAWFile = true);
    /** For deserialization */
    VvdObject() {}
    /** Destructor */
    ~VvdObject();
    /** Copy-constructor to  repare representationData pointer*/
    VvdObject(const VvdObject& other);

    /**
     * Main function
     * Used to return a new volume with a disk or octree representation
     */
    Volume* createVolume();

    /** @see Serializable::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Serializable::deserialize */
    virtual void deserialize(Deserializer& s);

public:
    VvdDataObjectBase* representationData_;             ///< object containing all needed information to load/save the representaion
    MetaDataContainer metaData_;                        ///< contains the meta data of the volume
    std::set<VolumeDerivedData*> derivedData_;          ///< contains the derived data of the volume

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
