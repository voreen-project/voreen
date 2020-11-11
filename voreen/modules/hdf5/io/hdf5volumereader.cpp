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

#include "hdf5volumereader.h"

#include "../utils/hdf5utils.h"
#include "volumediskhdf5.h"
#include "hdf5volumewriter.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "tgt/exception.h"
#include "tgt/assert.h"
#include "tgt/vector.h"

#include <boost/thread.hpp>

namespace { // anonymous helper functions

/// Try to identify a single volume from the given DataSet and add it to urls
static void collectVolume(const H5::DataSet& dataSet, const std::string& fileName, const std::string& pathInFile, std::vector<voreen::VolumeURL>& urls, bool separated) {
    H5::DataSpace dataSpace = dataSet.getSpace();
    hsize_t dimensions[4];
    dataSpace.getSimpleExtentDims(dimensions);
    int channels;
    switch(dataSpace.getSimpleExtentNdims()) {
        case 3: //One channel
            channels = 1;
            break;
        case 4: //Multiple channels
            channels = dimensions[3];
            break;
        default: //Not a suitable DataSet
            return;
    }
    if(separated) {
        for (int channel = 0; channel < channels; ++channel) {
            voreen::VolumeURL subURL("hdf5", fileName, "");
            subURL.addSearchParameter("channel", std::to_string(channel));
            subURL.addSearchParameter("path", pathInFile);
            subURL.getMetaDataContainer().addMetaData("Location", new voreen::StringMetaData(pathInFile));
            subURL.getMetaDataContainer().addMetaData("Channel", new voreen::IntMetaData(channel));
            subURL.getMetaDataContainer().addMetaData("Volume Dimensions", new voreen::IVec3MetaData(
                    static_cast<tgt::ivec3>(voreen::vec3HDF5ToTgt(dimensions))));
            urls.push_back(subURL);
        }
    }
    else {
        voreen::VolumeURL subURL("hdf5", fileName, "");
        subURL.addSearchParameter("path", pathInFile);
        subURL.addSearchParameter("channels", std::to_string(channels));
        subURL.getMetaDataContainer().addMetaData("Location", new voreen::StringMetaData(pathInFile));
        subURL.getMetaDataContainer().addMetaData("Channels", new voreen::IntMetaData(channels));
        subURL.getMetaDataContainer().addMetaData("Volume Dimensions", new voreen::IVec3MetaData(
                static_cast<tgt::ivec3>(voreen::vec3HDF5ToTgt(dimensions))));
        urls.push_back(subURL);
    }
}
/// Extract any possible volumes from the CommonFG and add them to urls
static void collectVolumes(
#if H5_VERSION_GE(1, 10, 1)
        const H5::Group& fg,
#else
        const H5::CommonFG& fg,
#endif
        const std::string& fileName, const std::string& pathInFile, std::vector<voreen::VolumeURL>& urls, bool separated) {
    for(hsize_t i=0; i<fg.getNumObjs(); ++i) {
        H5O_type_t type = fg.childObjType(i);
        switch(type) {
            case H5O_TYPE_GROUP:
                {
                    std::string groupName = fg.getObjnameByIdx(i);
                    collectVolumes(fg.openGroup(groupName), fileName, pathInFile + "/" + groupName, urls, separated);
                }
                break;
            case H5O_TYPE_DATASET:
                {
                    std::string dataSetName = fg.getObjnameByIdx(i);
                    collectVolume(fg.openDataSet(dataSetName), fileName, pathInFile + "/" + dataSetName, urls, separated);
                }
                break;
            default:
                continue;
        }
    }
}

} // anonymous helper functions

// ---------------------------------------------------------------------------------------------------------------------

namespace voreen {

const std::string HDF5VolumeReaderBase::loggerCat_ = "voreen.hdf5.HDF5VolumeReader";

HDF5VolumeReaderBase::HDF5VolumeReaderBase(bool separatedChannels)
    : VolumeReader()
    , separatedChannels_(separatedChannels)
{
    extensions_.push_back("hdf5");
    extensions_.push_back("h5");

    protocols_.push_back("hdf5");
}

VolumeList* HDF5VolumeReaderBase::read(const std::string& url) {
    std::vector<VolumeURL> urls = listVolumes(url);
    std::vector<std::unique_ptr<VolumeBase>> volumes;

    try {
        for(const auto& url : urls) {
            volumes.push_back(std::unique_ptr<VolumeBase>(read(url)));
        }
    } catch(tgt::IOException& e) {
        throw;
    }

    // Transfer ownership to output list.
    VolumeList* volumeList = new VolumeList();
    for(auto& volume : volumes) {
        volumeList->add(volume.release());
    }
    return volumeList;
}

VolumeBase* HDF5VolumeReaderBase::read(const VolumeURL& origin) {
    std::string fileName = origin.getPath();
    std::string inFilePath = origin.getSearchParameter("path");
    if(inFilePath.empty()) {
        inFilePath = HDF5VolumeWriter::VOLUME_DATASET_NAME;
        LWARNING("Could not find path in VolumeURL. Using default:" << inFilePath);
    }
    int firstChannel = 0;
    int numberOfChannels = 1;
    if(separatedChannels_) {
        std::string firstChannelStr = origin.getSearchParameter("channel");
        if(firstChannelStr.empty()) {
            firstChannel = 0;
            LWARNING("Could not find first channel in VolumeURL. Using default: " << firstChannel);
        } else {
            firstChannel = std::stoi(firstChannelStr);
        }
    }
    else {
        std::string numberOfChannelsStr = origin.getSearchParameter("channels");
        if(numberOfChannelsStr.empty()) {
            numberOfChannels = 1;
            LWARNING("Could not find number of channels in VolumeURL. Using default: " << numberOfChannels);
        } else {
            numberOfChannels = std::stoi(numberOfChannelsStr);
        }
    }
    LINFO("Loading " << fileName);

    // Open the in file volume.
    std::unique_ptr<HDF5FileVolume> fileVolume = HDF5FileVolume::openVolume(fileName, inFilePath, true);

    // Retrieve all needed information before handing the file volume over to volume disk
    std::vector<VolumeDerivedData*> derivedData = fileVolume->readDerivedData(firstChannel, numberOfChannels);
    std::unique_ptr<tgt::vec3> spacing(fileVolume->tryReadSpacing());
    std::unique_ptr<tgt::vec3> offset(fileVolume->tryReadOffset());
    std::unique_ptr<tgt::mat4> physicalToWorldTransformation(fileVolume->tryReadPhysicalToWorldTransformation());
    std::unique_ptr<RealWorldMapping> rwm(fileVolume->tryReadRealWorldMapping());
    std::string baseType = fileVolume->getBaseType();

    // Retrieve meta data.
    MetaDataContainer metaData = fileVolume->readMetaData();

    // Now create the VolumeDiskHDF5 by handing over fileVolume...
    VolumeRepresentation* volume = new VolumeDiskHDF5(std::move(fileVolume), firstChannel, numberOfChannels);
    // ... and construct the Volume.
    Volume* vol = new Volume(volume, tgt::vec3::one, tgt::vec3::zero, tgt::mat4::identity);

    // Set explicitly stored meta data and prioritize them.
    Vec3MetaData* spacingMetaData = dynamic_cast<Vec3MetaData*>(metaData.getMetaData(VolumeBase::META_DATA_NAME_SPACING));
    if (spacing) {
        // If both are available we take meta data value since the HDF5 spacing may suffer from precision loss.
        if (spacingMetaData) {
            tgt::vec3 relativeError = tgt::abs((spacingMetaData->getValue() - *spacing) / spacingMetaData->getValue());
            if(tgt::hor(tgt::greaterThan(relativeError, tgt::vec3(HDF5FileVolume::MM_PER_HDF5_UNIT_OF_LENGTH)))) {
                LWARNING("Spacing has been stored explicitly (" << *spacing
                                                                << ") and differs from value stored in MetaDataContainer ("
                                                                << spacingMetaData->getValue()
                                                                << "). Taking Spacing from MetaDataContainer.");
            }
            vol->setSpacing(spacingMetaData->getValue());
            metaData.removeMetaData(VolumeBase::META_DATA_NAME_SPACING);
        }
        else
            vol->setSpacing(*spacing);
    }
    else if (!spacingMetaData)
        LWARNING("HDF5 File contains no spacing information. Assuming (1, 1, 1)");

    Vec3MetaData* offsetMetaData = dynamic_cast<Vec3MetaData*>(metaData.getMetaData(VolumeBase::META_DATA_NAME_OFFSET));
    if (offset) {
        vol->setOffset(*offset);
        if (offsetMetaData && offsetMetaData->getValue() != *offset) {
            LWARNING("Offset has been stored explicitly (" << *offset
                                                           << ") and differs from value stored in MetaDataContainer ("
                                                           << offsetMetaData->getValue() << "). Taking explicit one.");
        }
        metaData.removeMetaData(VolumeBase::META_DATA_NAME_OFFSET);
    }
    else if (!offsetMetaData)
        LWARNING("HDF5 File contains no offset information. Assuming (0, 0, 0)");

    Mat4MetaData* physicalToWorldTransformationMetaData = dynamic_cast<Mat4MetaData*>(metaData.getMetaData(VolumeBase::META_DATA_NAME_TRANSFORMATION));
    if (physicalToWorldTransformation) {
        vol->setPhysicalToWorldMatrix(*physicalToWorldTransformation);
        if (physicalToWorldTransformationMetaData && physicalToWorldTransformationMetaData->getValue() != *physicalToWorldTransformation) {
            LWARNING("PhysicalToWorldTransformation has been stored explicitly (" << *physicalToWorldTransformation
                                                           << ") and differs from value stored in MetaDataContainer ("
                                                           << physicalToWorldTransformationMetaData->getValue()
                                                           << "). Taking explicit one.");
        }
        metaData.removeMetaData(VolumeBase::META_DATA_NAME_TRANSFORMATION);
    }
    else if (!physicalToWorldTransformationMetaData)
        LWARNING("HDF5 File contains no physical to world transformation. Assuming identity.");

    // Set other meta data and resolve conflicts.
    for (const std::string& key : metaData.getKeys()) {
        MetaDataBase* value = metaData.getMetaData(key);
        if (vol->getMetaDataContainer().hasMetaData(key)) {
            if (vol->getMetaDataContainer().getMetaData(key)->toString() != value->toString())
                LWARNING("Meta data '" << key << "' has been stored separately and differs from value stored in MetaDataContainer. Ignoring.");
        }
        else
            vol->getMetaDataContainer().addMetaData(key, value->clone());
    }

    // Set real world mapping...
    if(rwm) {
        // ... from file if available
        vol->setRealWorldMapping(*rwm);

        // Add all available derived data
        for(VolumeDerivedData* d : derivedData) {
            vol->addDerivedData(d);
        }
    } else {
        // ... otherwise a static denomalizing mapping.
        vol->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(baseType));

        // Note:
        // As previously the default RealWorldMapping had a scale of one, and all DerivedData created at that
        // time used the default RealWorldMapping, we explicitely don't set the DerivedData we read from the
        // HDF5 file and thus force them to be recreated.
        for(VolumeDerivedData* d : derivedData) {
            delete d;
        }
    }

    // Add origin information to the volume.
    vol->setOrigin(origin);

    // Finally return the volume.
    return vol;
}

std::vector<VolumeURL> HDF5VolumeReaderBase::listVolumes(const std::string& urlStr) const {
    VolumeURL url(urlStr);
    std::string filepath = url.getPath();
    std::vector<VolumeURL> urls;
    // HDF5: ---------------------------
    try {
        // Lock library for HDF5 operations
        boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

        H5::H5File file(filepath, H5F_ACC_RDONLY);
        collectVolumes(file, filepath, "", urls, separatedChannels_);
    } catch(H5::Exception& error) { // catch HDF5 exceptions
        LERROR(error.getFuncName() + ": " + error.getDetailMsg());
        throw tgt::IOException("An Error occured while reading file " + filepath);
    }
    return urls;
}

bool HDF5VolumeReaderBase::canSupportFileWatching() const {
    return true;
}

HDF5VolumeReader::HDF5VolumeReader()
    : HDF5VolumeReaderBase(true)
{
}

VolumeReader* HDF5VolumeReader::create(ProgressBar* progress) const {
    return new HDF5VolumeReader();
}

HDF5VolumeReaderOriginal::HDF5VolumeReaderOriginal()
    : HDF5VolumeReaderBase(false)
{
}

VolumeReader* HDF5VolumeReaderOriginal::create(ProgressBar* progress) const {
    return new HDF5VolumeReaderOriginal();
}

} // namespace voreen
