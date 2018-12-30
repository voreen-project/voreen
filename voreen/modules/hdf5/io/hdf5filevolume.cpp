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

#include "hdf5filevolume.h"

#include "../utils/hdf5utils.h" //Conversion and hdf5liblock
#include "hdf5volumewriter.h" //Attribute names and conversion constants
#include "voreen/core/datastructures/volume/volumefactory.h"

#include "tgt/filesystem.h"

#include <boost/thread.hpp>

namespace voreen {

const std::string HDF5FileVolume::loggerCat_("voreen.hdf5.HDF5FileVolume");

const std::string HDF5FileVolume::SPACING_ATTRIBUTE_NAME("element_size_um");
const std::string HDF5FileVolume::OFFSET_ATTRIBUTE_NAME("offset_um");
const std::string HDF5FileVolume::PHYSICAL_TO_WORLD_TRANSFORMATION_ATTRIBUTE_NAME("physical_to_world_transformation");
const std::string HDF5FileVolume::METADATA_STRING_NAME("metaData");
const std::string HDF5FileVolume::REALWORLDMAPPING_SCALE_ATTRIBUTE_NAME("realWorldMapping_scale");
const std::string HDF5FileVolume::REALWORLDMAPPING_OFFSET_ATTRIBUTE_NAME("realWorldMapping_offset");
const std::string HDF5FileVolume::REALWORLDMAPPING_UNIT_ATTRIBUTE_NAME("realWorldMapping_unit");
const std::string HDF5FileVolume::REPRESENTATION_MINMAX_ATTRIBUTE_NAME("representation_minmax");
const std::string HDF5FileVolume::REPRESENTATION_HISTOGRAMINTENSITY_ATTRIBUTE_NAME("representation_histogramintensity");
const std::string HDF5FileVolume::REPRESENTATION_HISTOGRAMINTENSITY_METADATA_ATTRIBUTE_NAME("representation_histogramintensity_metadata");
const std::string HDF5FileVolume::REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_ATTRIBUTE_NAME("representation_histogramintensitygradient");
const std::string HDF5FileVolume::REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_METADATA_ATTRIBUTE_NAME("representation_histogramintensitygradient_metadata");
const std::string HDF5FileVolume::REPRESENTATION_PREVIEW_ATTRIBUTE_NAME("representation_preview");
const std::string HDF5FileVolume::REPRESENTATION_PREVIEW_METADATA_ATTRIBUTE_NAME("representation_preview_metadata");
// Convention seems to use µm. => 0.001 mm per µm.
const float HDF5FileVolume::MM_PER_HDF5_UNIT_OF_LENGTH(0.001f);

std::unique_ptr<HDF5FileVolume> HDF5FileVolume::openVolume(const std::string& fileName, const std::string& volumeLocation, bool readOnly) {
    // We absolutely need to make sure that we lock the library during HDF5 operations.
    // Therefore we create the HDF5 objects in heap memory.
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Open the file ======================================================================================
    std::unique_ptr<H5::H5File> file = nullptr;

    try {
        unsigned int flags = 0;
        if(readOnly) {
            flags |= H5F_ACC_RDONLY;
        } else {
            flags |= H5F_ACC_RDWR;
        }
        file = std::unique_ptr<H5::H5File>(new H5::H5File(fileName, flags));
    } catch(H5::Exception error) {
        throw tgt::IOException("Error opening hdf5 file:" + error.getFuncName() + ": " + error.getDetailMsg());
    }

    // Open the dataSet ===================================================================================
    std::unique_ptr<H5::DataSet> dataSet = nullptr;

    try {
        dataSet = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet(volumeLocation)));
    } catch(H5::Exception error) {
        throw tgt::IOException("Error opening dataset:" + error.getFuncName() + ": " + error.getDetailMsg());
    }
    return std::unique_ptr<HDF5FileVolume>(new HDF5FileVolume(std::move(file), std::move(dataSet), fileName, volumeLocation));
}

std::unique_ptr<HDF5FileVolume> HDF5FileVolume::createVolume(const std::string& fileName, const std::string& volumeLocation, const std::string& baseType, tgt::svec3 dimensions, size_t numChannels, bool truncateFile, int deflateLevel, tgt::svec3 chunkSize, bool shuffle) {
    tgtAssert(numChannels > 0, "Number of channels must be > 0");
    tgtAssert(dimensions.x >= chunkSize.x && dimensions.y >= chunkSize.y && dimensions.z >= chunkSize.z, "Chunksize cannot be greater than volume dimensions.");
    tgtAssert(deflateLevel >= 0 && deflateLevel <= 9, "Deflate level must be in [0, 9].");

    // We absolutely need to make sure that we lock the library during HDF5 operations.
    // Therefore we create the HDF5 objects in heap memory.
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);


    // Open the file ======================================================================================
    std::unique_ptr<H5::H5File> file = nullptr;

    try {
        // Flags for opening/creating the file
        // Although the documentation states that the constructor of H5File "opens or creates" the file
        // it only decides what to to based on the flags given as second argument...
        // Because of that if we explicitly want to truncate the file, the file does not exist yet or is
        // not a hdf5 file, we set the truncate flag.
        unsigned int flags = 0;
        if(truncateFile || !tgt::FileSystem::fileExists(fileName) || !H5::H5File::isHdf5(fileName)) {
            flags |= H5F_ACC_TRUNC;
        } else {
            // If we do not truncate the file, we have to specify that we open it in Read/Write mode.
            flags |= H5F_ACC_RDWR;
        }

        // Setup file access property list for file creation.
        H5::FileAccPropList accList(H5::FileAccPropList::DEFAULT);

        // We require at least format version 1.8 to be able to write arbitrarily large attribute data.
        // Wrapper not available in 1.8.13, so we use the c version:
        //accList.setLibverBounds(H5F_LIBVER_18, H5F_LIBVER_LATEST);
        // ...
        // Unfortunately libhdf5 v.1.10 removed H5F_LIBVER_18 (which for 1.8 was an alias for H5F_LIBVER_LATEST)
        // so for now we set the minimum version to H5F_LIBVER_LATEST and hope that there are not incompatiblitities
        // => TODO: Set the min version accordingly if there is (backwarts compatible) support in a future version
        H5Pset_libver_bounds(accList.getId(), H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

        // Create the hdf5 file using the propertylists created earlier.
        file = std::unique_ptr<H5::H5File>(new H5::H5File(fileName, flags, H5::FileCreatPropList::DEFAULT, accList));
    } catch(H5::Exception error) {
        throw tgt::IOException("Error opening/creating hdf5 file:" + error.getFuncName() + ": " + error.getDetailMsg());
    }


    // Create the dataSet =================================================================================
    std::unique_ptr<H5::DataSet> dataSet = nullptr;

    try {
        // Prepare to create the volume inside the file
        H5::DataType h5type = getPredType(baseType);
        H5::DataSpace dataSpace = createDataSpace(dimensions, numChannels);

        // Set up properties for the volume data set to be created.
        H5::DSetCreatPropList propList;

        // Set Shuffle if wanted.
        if(shuffle) {
            propList.setShuffle();
        }

        // Set chunking
        // We need to set chunking even if we want to use the default chunk size
        // because otherwise the chunk size of propList will be uninitiated
        // apparently.
        // We also cannot create a copy of the default property list because the
        // HDF5 lib does not really create a copy (although it states that it does),
        // but still references the original interal c-style property list...
        // Because of that we:

        // ... Set the chunk size to be an xy-slice of the volume manually, if
        // none was specified.
        if(chunkSize == tgt::svec3::zero) {
            chunkSize = dimensions;
            chunkSize.z = 1;
        }
        // ... And set the chunk size of the previously created property list.
        tgtAssert(chunkSize.x > 0 && chunkSize.y > 0 && chunkSize.z > 0, "All dimensions of chunkSize have to be > 0");
        hsize_t chunkdimHDF5[4];
        // If we have a multidimensionale Volume, we definitely want the channels separated
        vec4TgtToHDF5(tgt::svec4(1, chunkSize), chunkdimHDF5);
        if(numChannels == 1) {
            propList.setChunk(3, chunkdimHDF5);
        } else {
            propList.setChunk(4, chunkdimHDF5);
        }

        // Set compression level
        // setDeflate(0) still adds a filter altough it does not compress.
        if(deflateLevel != 0) {
            propList.setDeflate(deflateLevel);
        }

        // Finally try to create the data set inside the file
        dataSet = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(volumeLocation, h5type, dataSpace, propList)));
    } catch(H5::Exception error) {
        throw tgt::IOException("Error constructing dataset:" + error.getFuncName() + ": " + error.getDetailMsg());
    }

    return std::unique_ptr<HDF5FileVolume>(new HDF5FileVolume(std::move(file), std::move(dataSet), fileName, volumeLocation));
}

HDF5FileVolume::~HDF5FileVolume() {
    // Lock before calling destructors on hdf5 objects
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    dataSet_ = nullptr;
    file_ = nullptr;
}

const std::string& HDF5FileVolume::getFileName() const {
    return fileName_;
}

const std::string& HDF5FileVolume::getVolumeLocation() const {
    return volumeLocation_;
}

tgt::svec3 HDF5FileVolume::getDimensions() const {
    return dimensions_;
}

tgt::svec3 HDF5FileVolume::getChunkSize() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    H5::DSetCreatPropList propList = dataSet_->getCreatePlist();
    // If we have multiple channels, the channel dimension will be fastest changing => in the back
    hsize_t chunkdimHDF5[3];
    try {
        propList.getChunk(3, chunkdimHDF5);
    } catch (H5::PropListIException e) {
        // If there is an error retrieving the chunk dim: Assume the volume is not split into chunks
        return dimensions_;
    }
    return vec3HDF5ToTgt(chunkdimHDF5);
}

tgt::vec3* HDF5FileVolume::tryReadSpacing() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Check if there is a spacing attribute:
    if(!dataSet_->attrExists(SPACING_ATTRIBUTE_NAME)) {
        //LWARNING("HDF5 DataSet does not contain spacing information. Assuming spacing (1,1,1).");
        return nullptr;
    }
    try {
        // Convert to Voreen's mm.
        return new tgt::vec3(readVec3Attribute<float>(*dataSet_, SPACING_ATTRIBUTE_NAME) * MM_PER_HDF5_UNIT_OF_LENGTH);

    } catch(H5::AttributeIException error) {
        LERROR("Error reading spacing attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }
}
void HDF5FileVolume::writeSpacing(const tgt::vec3& spacing) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // We convert voreen's mm to hdf5s unit of length.
    writeVec3Attribute(*dataSet_, SPACING_ATTRIBUTE_NAME, spacing/MM_PER_HDF5_UNIT_OF_LENGTH);
}

tgt::vec3* HDF5FileVolume::tryReadOffset() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Check if there is an offset attribute:
    if(!dataSet_->attrExists(OFFSET_ATTRIBUTE_NAME)) {
        //LWARNING("HDF5 DataSet does not contain offset information. Assuming offset (0,0,0).");
        return nullptr;
    }
    try {
        // Convert to Voreen's mm.
        return new tgt::vec3(readVec3Attribute<float>(*dataSet_, OFFSET_ATTRIBUTE_NAME) * MM_PER_HDF5_UNIT_OF_LENGTH);

    } catch(H5::AttributeIException error) {
        LERROR("Error reading offset attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        //LERROR("\tAssuming offset (0,0,0).");
        return nullptr;
    }
}
void HDF5FileVolume::writeOffset(const tgt::vec3& offset) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // We convert voreen's mm to hdf5s unit of length.
    writeVec3Attribute(*dataSet_, OFFSET_ATTRIBUTE_NAME, offset/MM_PER_HDF5_UNIT_OF_LENGTH);
}

tgt::mat4* HDF5FileVolume::tryReadPhysicalToWorldTransformation() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    if(!dataSet_->attrExists(PHYSICAL_TO_WORLD_TRANSFORMATION_ATTRIBUTE_NAME)) {
        LDEBUG("No Physical-to-world matrix information found.");
        return nullptr;
    }

    tgt::mat4* transformation = new tgt::mat4;
    readArrayAttribute(*dataSet_, PHYSICAL_TO_WORLD_TRANSFORMATION_ATTRIBUTE_NAME, transformation->elem, tgt::mat4::size);
    return transformation;
}

void HDF5FileVolume::writePhysicalToWorldTransformation(const tgt::mat4& transformation) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    writeArrayAttribute(*dataSet_, PHYSICAL_TO_WORLD_TRANSFORMATION_ATTRIBUTE_NAME, transformation.elem, tgt::mat4::size);
}


RealWorldMapping* HDF5FileVolume::tryReadRealWorldMapping() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Check if there is no RWM information at all..
    if(   !dataSet_->attrExists(REALWORLDMAPPING_SCALE_ATTRIBUTE_NAME)
       && !dataSet_->attrExists(REALWORLDMAPPING_OFFSET_ATTRIBUTE_NAME)
       && !dataSet_->attrExists(REALWORLDMAPPING_UNIT_ATTRIBUTE_NAME)) {
        LDEBUG("No RealWorldMapping information found.");
        return nullptr;
    }
    try {
        // Read values for real world mapping from the data set.
        float scale = readScalarAttribute<float>(*dataSet_, REALWORLDMAPPING_SCALE_ATTRIBUTE_NAME);
        float offset = readScalarAttribute<float>(*dataSet_, REALWORLDMAPPING_OFFSET_ATTRIBUTE_NAME);
        std::string unit = readStringAttribute(*dataSet_, REALWORLDMAPPING_UNIT_ATTRIBUTE_NAME);

        // real world mapping is complete, so write it to rwm.
        return new RealWorldMapping(scale, offset, unit);

    } catch(H5::AttributeIException error) {
        // If there are no attributes for RWM at all we would not have gotten here, so something has really gone wrong:
        // Issue an error message!
        LERROR("Error reading RealWorldMapping attributes: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }
}

void HDF5FileVolume::writeRealWorldMapping(const RealWorldMapping& rwm) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    writeScalarAttribute(*dataSet_, REALWORLDMAPPING_SCALE_ATTRIBUTE_NAME , rwm.getScale() );
    writeScalarAttribute(*dataSet_, REALWORLDMAPPING_OFFSET_ATTRIBUTE_NAME, rwm.getOffset());
    writeStringAttribute(*dataSet_, REALWORLDMAPPING_UNIT_ATTRIBUTE_NAME  , rwm.getUnit()  );
}

std::string HDF5FileVolume::getBaseType() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    return getBaseTypeFromDataType(dataSet_->getDataType());
}

VolumeMinMax* HDF5FileVolume::tryReadVolumeMinMax(size_t channel) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Check if min-max information is present and exit early if not.
    if(!dataSet_->attrExists(REPRESENTATION_MINMAX_ATTRIBUTE_NAME)) {
        LDEBUG("No VolumeMinMax information found.");
        return nullptr;
    }
    // VolumeMinMax is serialized into 4 floats:
    // min, max, minNormalized, maxNormalized.
    const size_t perChannelDataSize = 4;

    // First read VolumeMinMax data for all channels...
    std::vector<float> values(getNumberOfChannels() * perChannelDataSize);
    try {
        readArrayAttribute(*dataSet_, REPRESENTATION_MINMAX_ATTRIBUTE_NAME, values.data(), values.size());
    } catch(H5::AttributeIException error) {
        LERROR("Error reading VolumeMinMax attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }
    // ... then find the data for the channel we are looking for.
    float* channelData = &values.data()[perChannelDataSize * channel];

    // Before reading  make sure we only read from the memory previously allocated.
    tgtAssert(channelData + perChannelDataSize <= values.data() + values.size(), "VolumeMinMax unpack error.");
    // We are safe: Get the values and construct the VolumeMinMax!
    return new VolumeMinMax(channelData[0], channelData[1], channelData[2], channelData[3]);
}


void HDF5FileVolume::writeVolumeMinMax(const VolumeMinMax* mm) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    tgtAssert(mm, "No VolumeMinMax.");

    std::vector<float> values;
    // We need to store the data for each channel.
    for(size_t channel = 0; channel < mm->getNumChannels(); ++channel) {
        // For each channel there are 4 floats to store: min and max, each as real world and normalized version.
        values.push_back(mm->getMin(channel));
        values.push_back(mm->getMax(channel));
        values.push_back(mm->getMinNormalized(channel));
        values.push_back(mm->getMaxNormalized(channel));
    }
    writeArrayAttribute(*dataSet_, REPRESENTATION_MINMAX_ATTRIBUTE_NAME, values.data(), values.size());
}

VolumeHistogramIntensity* HDF5FileVolume::tryReadVolumeHistogramIntensity(size_t channel) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Check if meta data xor bucket data is available => data is incomplete so we issue a warning.
    if(!dataSet_->attrExists(REPRESENTATION_HISTOGRAMINTENSITY_ATTRIBUTE_NAME)
            ^ !dataSet_->attrExists(REPRESENTATION_HISTOGRAMINTENSITY_METADATA_ATTRIBUTE_NAME)) {
        LWARNING("Incomplete VolumeHistogramIntensity data.");
        return nullptr;
    }
    // Either both or none of the attributes are present. We check for the metadata attribute early
    // and return if it is not present.
    if(!dataSet_->attrExists(REPRESENTATION_HISTOGRAMINTENSITY_METADATA_ATTRIBUTE_NAME)) {
        LDEBUG("No VolumeHistogramIntensity information found.");
        return nullptr;
    }

    // VolumeHistogramIntensity's meta data is serialized into 3 floats:
    // min, max and number of buckets.
    const size_t perChannelMetaDataSize = 3;

    // First read all meta data...
    std::vector<float> metaData(getNumberOfChannels() * perChannelMetaDataSize);
    try {
        readArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITY_METADATA_ATTRIBUTE_NAME, metaData.data(), metaData.size());
    } catch(H5::AttributeIException error) {
        LERROR("Error reading VolumeHistogramIntensity-MetaData attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }
    // ... then iterate over all preceeding channels and collect numberOfBuckets
    // to find this channels start position in the data attribute.
    size_t startBucketPos = 0;
    for(size_t c = 0; c<channel; ++c) {
        startBucketPos += static_cast<size_t>(metaData.at(perChannelMetaDataSize*c + 2));
    }
    // Calculate the size of the attributes data
    size_t bufferSize = 0;
    for(size_t c = 0; c<getNumberOfChannels(); ++c) {
        bufferSize += static_cast<size_t>(metaData.at(perChannelMetaDataSize*c + 2)); /* numBuckets */
    }
    //Now read this channels metadata:
    float minValue    =                     metaData.at(perChannelMetaDataSize*channel + 0 ) ;
    float maxValue    =                     metaData.at(perChannelMetaDataSize*channel + 1 ) ;
    size_t numBuckets = static_cast<size_t>(metaData.at(perChannelMetaDataSize*channel + 2 ));

    Histogram1D h(minValue, maxValue, numBuckets);
    std::vector<uint64_t> values(bufferSize);
    try {
        readArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITY_ATTRIBUTE_NAME, values.data(), values.size());
    } catch(H5::AttributeIException error) {
        LERROR("Error reading VolumeHistogramIntensity attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }

    // Finally collect actual histogram data:
    // First check that we will not read from data outside of the vector:
    tgtAssert(startBucketPos + numBuckets <= values.size(), "Histogram unpack error.");
    // Now that we are sure the memory is safe, find the first bucket and iterate over the rest to fill the previously constructed Histogram1D
    uint64_t* startBucket = &values.data()[startBucketPos];
    for(size_t bucket = 0; bucket < numBuckets; ++bucket) {
        h.increaseBucket(bucket, startBucket[bucket]);
    }
    return new VolumeHistogramIntensity(h);
}

void HDF5FileVolume::writeVolumeHistogramIntensity(const VolumeHistogramIntensity* hi) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    tgtAssert(hi, "No VolumeHistogramIntensity.");

    std::vector<uint64_t> values;
    std::vector<float> metaData;
    // Store each channel of the histogram after the previous one
    for(size_t channel = 0; channel < hi->getNumChannels(); ++channel) {
        const Histogram1D h = hi->getHistogram(channel);

        // For meta data we save the min and max value and the number of buckets.
        metaData.push_back(h.getMinValue());
        metaData.push_back(h.getMaxValue());
        // Min and max are floats so we store numBuckets as float so that we dont have to create another hdf5 attribute.
        metaData.push_back(h.getNumBuckets());

        // Iterate over the whole histogram and store each bucket value in the vector.
        for(size_t bucket = 0; bucket < h.getNumBuckets(); ++bucket) {
            values.push_back(h.getBucket(bucket));
        }
    }
    writeArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITY_METADATA_ATTRIBUTE_NAME, metaData.data(), metaData.size());
    writeArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITY_ATTRIBUTE_NAME, values.data(), values.size());
}

VolumeHistogramIntensityGradient* HDF5FileVolume::tryReadVolumeHistogramIntensityGradient(size_t channel) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Check if meta data xor bucket data is available => data is incomplete so we issue a warning.
    if(!dataSet_->attrExists(REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_ATTRIBUTE_NAME)
            ^ !dataSet_->attrExists(REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_METADATA_ATTRIBUTE_NAME)) {
        LWARNING("Incomplete VolumeHistogramIntensityGradient data.");
        return nullptr;
    }
    // Either both or none of the attributes are present. We check for the metadata attribute early
    // and return if it is not present.
    if(!dataSet_->attrExists(REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_METADATA_ATTRIBUTE_NAME)) {
        LDEBUG("No VolumeHistogramIntensityGradient information found.");
        return nullptr;
    }

    // VolumeHistogramIntensity's meta data is serialized into 6 floats:
    // min, max, number of channels for intensity dimension
    // min, max, number of channels for gradient dimension
    const size_t perChannelMetaDataSize = 6;

    // First read all meta data...
    std::vector<float> metaData(getNumberOfChannels() * perChannelMetaDataSize);
    try {
        readArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_METADATA_ATTRIBUTE_NAME, metaData.data(), metaData.size());
    } catch(H5::AttributeIException error) {
        LERROR("Error reading VolumeHistogramIntensityGradient-MetaData attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }
    // ... then iterate over all preceeding channels and collect numberOfBuckets
    // to find this channels start position in the data attribute.
    size_t maxBucketPos = 0;
    for(size_t c = 0; c<channel; ++c) {
        maxBucketPos += 1; //Account for maxbucket of other channels.
        maxBucketPos += static_cast<size_t>(metaData.at(perChannelMetaDataSize*c + 2)); /* numBuckets0 */
    }
    // Calculate the size of the attributes data
    size_t bufferSize = 0;
    for(size_t c = 0; c<getNumberOfChannels(); ++c) {
        bufferSize += 1; //Account for maxbucket of other channels.
        bufferSize += static_cast<size_t>(metaData.at(perChannelMetaDataSize*c + 2)) /* numBuckets0 */
                    * static_cast<size_t>(metaData.at(perChannelMetaDataSize*c + 5)) /* numBuckets1 */;
    }
    // Now read this channels metadata: We have min, max and the number of buckets for each of the two dimensions 0 and 1.
    float minValue0    =                     metaData.at(perChannelMetaDataSize*channel + 0) ;
    float maxValue0    =                     metaData.at(perChannelMetaDataSize*channel + 1) ;
    size_t numBuckets0 = static_cast<size_t>(metaData.at(perChannelMetaDataSize*channel + 2));
    float minValue1    =                     metaData.at(perChannelMetaDataSize*channel + 3) ;
    float maxValue1    =                     metaData.at(perChannelMetaDataSize*channel + 4) ;
    size_t numBuckets1 = static_cast<size_t>(metaData.at(perChannelMetaDataSize*channel + 5));

    // Construct the Histogram2D with the values just read...
    Histogram2D h(minValue0, maxValue0, numBuckets0, minValue1, maxValue1, numBuckets1);

    // ... and read the actual values to fill the bucket.
    std::vector<uint64_t> values(bufferSize);
    try {
        readArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_ATTRIBUTE_NAME, values.data(), values.size());
    } catch(H5::AttributeIException error) {
        LERROR("Error reading VolumeHistogramIntensityGradient attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }

    // The first positions contains the maxvalue of all buckets.
    uint64_t maxBucket = values.at(maxBucketPos);
    // Histogram starts just after maxBucket
    size_t startBucketPos = maxBucketPos + 1;

    // Finally collect actual histogram data:
    // First check that we will not read from data outside of the vector:
    tgtAssert(startBucketPos + numBuckets0*numBuckets1 <= values.size(), "Histogram unpack error.");
    // Now that we are sure the memory is safe, find the first bucket and iterate over the rest to fill the previously constructed Histogram2D
    uint64_t* startBucket = &values.data()[startBucketPos];
    for(size_t bucket = 0; bucket < numBuckets0*numBuckets1; ++bucket) {
        h.increaseBucket(bucket, startBucket[bucket]);
    }

    // Construct the VolumeHistogramIntensityGradient from the Histogram2D and the maxBucket value.
    return new VolumeHistogramIntensityGradient(h, maxBucket);
}

void HDF5FileVolume::writeVolumeHistogramIntensityGradient(const VolumeHistogramIntensityGradient* hig) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    tgtAssert(hig, "No VolumeHistogramIntensityGradient.");

    std::vector<uint64_t> values;
    std::vector<float> metaData;
    // Store each channel of the histogram after the previous one
    for(size_t channel = 0; channel < hig->getNumChannels(); ++channel) {
        const Histogram2D h = hig->getHistogram();

        // For meta data we save the min and max value as well as the number of buckets for both dimensions:
        metaData.push_back(h.getMinValue(0));
        metaData.push_back(h.getMaxValue(0));
        // Min and max are floats so we store numBuckets as float so that we dont have to create another hdf5 attribute.
        metaData.push_back(h.getNumBuckets(0));

        metaData.push_back(h.getMinValue(1));
        metaData.push_back(h.getMaxValue(1));
        // Min and max are floats so we store numBuckets as float so that we dont have to create another hdf5 attribute.
        metaData.push_back(h.getNumBuckets(1));

        // Put max bucket data just in front of the histogram.
        // This should probably considered meta data, but the floats might not be able to save the content of a uint64_t properly...
        values.push_back(hig->getMaxBucket(channel));

        // Iterate over the whole histogram and store each bucket value in the vector.
        for(size_t bucket = 0; bucket < h.getNumBuckets(); ++bucket) {
            values.push_back(h.getBucket(bucket));
        }
    }
    writeArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_METADATA_ATTRIBUTE_NAME, metaData.data(), metaData.size());
    writeArrayAttribute(*dataSet_, REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_ATTRIBUTE_NAME, values.data(), values.size());
}

VolumePreview* HDF5FileVolume::tryReadVolumePreview() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    // Check if meta data xor bucket data is available => data is incomplete so we issue a warning.
    if(!dataSet_->attrExists(REPRESENTATION_PREVIEW_ATTRIBUTE_NAME)
            ^ !dataSet_->attrExists(REPRESENTATION_PREVIEW_METADATA_ATTRIBUTE_NAME)) {
        LWARNING("Incomplete VolumePreview data.");
        return nullptr;
    }
    // Either both or none of the attributes are present. We check for the metadata attribute early
    // and return if it is not present.
    if(!dataSet_->attrExists(REPRESENTATION_PREVIEW_ATTRIBUTE_NAME)) {
        LDEBUG("No VolumePreview information found.");
        return nullptr;
    }
    // meta data is serialized into 2 uint64_t!
    const size_t metaDataSize = 2;
    std::vector<uint64_t> metaData(metaDataSize);
    try {
        readArrayAttribute(*dataSet_, REPRESENTATION_PREVIEW_METADATA_ATTRIBUTE_NAME, metaData.data(), metaData.size());
    } catch(H5::AttributeIException error) {
        LERROR("Error reading VolumePreview-MetaData attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }
    // Meta data consists of width and height, so read them...
    const uint64_t width = metaData.at(0);
    const uint64_t height = metaData.at(1);

    // ... and prepare the storage for the actual preview values.
    std::vector<unsigned char> values(width*height);
    try {
        readArrayAttribute(*dataSet_, REPRESENTATION_PREVIEW_ATTRIBUTE_NAME, values.data(), values.size());
    } catch(H5::AttributeIException error) {
        LERROR("Error reading VolumePreview attribute: " + error.getFuncName() + ": " + error.getDetailMsg());
        return nullptr;
    }
    return new VolumePreview(height, values);
}

void HDF5FileVolume::writeVolumePreview(const VolumePreview* preview) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    tgtAssert(preview, "No VolumePreview.");

    // Check if the volumepreview's data has and unfinished column. This is not allowed.
    tgtAssert((preview->getData().size() % preview->getHeight()) == 0, "Volumepreview is not rectangular.");

    // Negative height does not make sense:
    tgtAssert(preview->getHeight() > 0, "Volumepreview's height is not > 0.");

    // Because we know that the preview is rectangular we can get the actual width.
    uint64_t width = preview->getData().size()/preview->getHeight();

    // Store width and height as meta data:
    std::vector<uint64_t> metaData;
    metaData.push_back(width);
    metaData.push_back(preview->getHeight());

    writeArrayAttribute(*dataSet_, REPRESENTATION_PREVIEW_METADATA_ATTRIBUTE_NAME, metaData.data(), metaData.size());
    writeArrayAttribute(*dataSet_, REPRESENTATION_PREVIEW_ATTRIBUTE_NAME, preview->getData().data(), preview->getData().size());
}

std::vector<VolumeDerivedData*> HDF5FileVolume::readDerivedData(size_t channel) const {
    // Just try to read all currently known DerivedData types (except VolumeHash)
    // and if found, add them to the vector which will be returned.
    std::vector<VolumeDerivedData*> data;

    VolumeDerivedData* volumeMinMax = tryReadVolumeMinMax(channel);
    if(volumeMinMax) {
        data.push_back(volumeMinMax);
    }

    VolumeHistogramIntensity* volumeHistogramIntensity = tryReadVolumeHistogramIntensity(channel);
    if(volumeHistogramIntensity) {
        data.push_back(volumeHistogramIntensity);
    }

    VolumeHistogramIntensityGradient* volumeHistogramIntensityGradient = tryReadVolumeHistogramIntensityGradient(channel);
    if(volumeHistogramIntensityGradient) {
        data.push_back(volumeHistogramIntensityGradient);
    }

    // VolumePreview only shows a preview of the first channel.
    // Thus when we extract DerivedData for a single channel other than the first that volume preview does not make sense.
    if(channel == 0) {
        VolumePreview* volumePreview = tryReadVolumePreview();
        if(volumePreview) {
            data.push_back(volumePreview);
        }
    }

    return data;
}

void HDF5FileVolume::writeDerivedData(const VolumeBase* vol) const {
    // TODO: Iterate over all derived data of vol and thus notice if there is a
    // new type. In that case issue a warning.
    if(vol->hasDerivedData<VolumeMinMax>()) {
        writeVolumeMinMax(vol->getDerivedData<VolumeMinMax>());
    }
    if(vol->hasDerivedData<VolumeHistogramIntensity>()) {
        writeVolumeHistogramIntensity(vol->getDerivedData<VolumeHistogramIntensity>());
    }
    if(vol->hasDerivedData<VolumeHistogramIntensityGradient>()) {
        writeVolumeHistogramIntensityGradient(vol->getDerivedData<VolumeHistogramIntensityGradient>());
    }
    if(vol->hasDerivedData<VolumePreview>()) {
        writeVolumePreview(vol->getDerivedData<VolumePreview>());
    }
}

MetaDataContainer HDF5FileVolume::readMetaData() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    if (!dataSet_->attrExists(METADATA_STRING_NAME))
        return MetaDataContainer();

    std::stringstream stream;
    stream << readStringAttribute(*dataSet_, METADATA_STRING_NAME);

    XmlDeserializer deserialier;
    deserialier.read(stream);

    MetaDataContainer metaData;
    deserialier.deserialize("metaData", metaData);
    return metaData;
}

void HDF5FileVolume::writeMetaData(const VolumeBase* vol) const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    const Volume* volume = dynamic_cast<const Volume*>(vol);
    if(!volume)
        return;

    std::stringstream stream;
    XmlSerializer serializer;
    serializer.serialize("metaData", volume->getMetaDataContainer());
    serializer.write(stream);

    writeStringAttribute(*dataSet_, METADATA_STRING_NAME, stream.str());
}

VolumeRAM* HDF5FileVolume::loadVolume(size_t channel) const {
    return loadBrick(tgt::vec3(0,0,0), dimensions_, channel);
}

void HDF5FileVolume::writeVolume(const VolumeRAM* vol) const  {
    writeVolume(vol, 0, getNumberOfChannels());
}

void HDF5FileVolume::writeVolume(const VolumeRAM* vol, const size_t firstChannel, const size_t numberOfChannels) const  {
    tgtAssert(vol->getDimensions() == dimensions_, "Volume dimensions differ");
    writeBrick(vol, tgt::svec3::zero, firstChannel, numberOfChannels);
}


VolumeRAM* HDF5FileVolume::loadSlices(const size_t firstZSlice, const size_t lastZSlice, size_t channel) const {
    tgt::vec3 offset = tgt::vec3(0, 0, firstZSlice);
    tgt::vec3 dimensions = getDimensions();
    // lastZSlice is included => +1
    dimensions.z = lastZSlice - firstZSlice + 1;
    return loadBrick(offset, dimensions, channel);
}

void HDF5FileVolume::writeSlices(const VolumeRAM* vol, const size_t firstSlice) const {
    writeSlices(vol, firstSlice, 0, getNumberOfChannels());
}

void HDF5FileVolume::writeSlices(const VolumeRAM* vol, const size_t firstSlice, const size_t firstChannel, const size_t numberOfChannels) const {
    tgtAssert(vol->getDimensions().x == dimensions_.x, "Volume dimensions differ in x.");
    tgtAssert(vol->getDimensions().y == dimensions_.y, "Volume dimensions differ in y.");
    tgtAssert(vol->getDimensions().z + firstSlice <= dimensions_.z, "Slices do not fit into file volume with the given offset.");
    writeBrick(vol, tgt::svec3(0, 0, firstSlice), firstChannel, numberOfChannels);
}

VolumeRAM* HDF5FileVolume::loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions, size_t channel) const {
    // HDF5-cxx-libs are NOT threadsafe!
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    tgtAssert(dimensions.x + offset.x <= dimensions_.x, "Brick does not fit into file volume with the given offset in x dimension.");
    tgtAssert(dimensions.y + offset.y <= dimensions_.y, "Brick does not fit into file volume with the given offset in y dimension.");
    tgtAssert(dimensions.z + offset.z <= dimensions_.z, "Brick does not fit into file volume with the given offset in z dimension.");

    const size_t numChannels = getNumberOfChannels();
    tgtAssert(channel < numChannels, "Invalid channel");

    try {
        //Select the hyperslab (Brick + fourth dimension)

        //Just in case there are more than one channel: Allocate 4 dimensions.
        //If there is only one channel the 4. dimension will not do anything.
        hsize_t count[4];
        vec4TgtToHDF5(tgt::svec4(1, dimensions), count);
        hsize_t start[4];
        vec4TgtToHDF5(tgt::svec4(channel, offset), start);

        H5::DataSpace fileSpace(dataSet_->getSpace());
        fileSpace.selectHyperslab(H5S_SELECT_SET, count, start);

        //Memory space does not have an offset
        H5::DataSpace memSpace(3, count);

        //Create the VolumeRAM and write the dataset selection to it.
        H5::DataType h5type = dataSet_->getDataType();
        VolumeRAM* data = VolumeFactory().create(getBaseTypeFromDataType(h5type), dimensions);
        if(!data) {
            throw tgt::IOException("Could not create VolumeRAM of type " + getBaseTypeFromDataType(h5type));
        }
        dataSet_->read(data->getData(), h5type, memSpace, fileSpace);
        return data;
    } catch(H5::Exception error) { // catch HDF5 exceptions
        LERROR(error.getFuncName() + ": " + error.getDetailMsg());
        throw tgt::IOException("An Error occured while reading volume from file " + getFileName());
    }
}

void HDF5FileVolume::writeBrick(const VolumeRAM* vol, const tgt::svec3& offset) const {
    writeBrick(vol, offset, 0, getNumberOfChannels());
}

void HDF5FileVolume::writeBrick(const VolumeRAM* vol, const tgt::svec3& offset, const size_t firstChannel, const size_t numberOfChannels) const {
    // HDF5-cxx-libs are NOT threadsafe!
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    tgtAssert(vol, "No VolumeRAM");
    tgtAssert(vol->getDimensions().x + offset.x <= dimensions_.x, "Brick does not fit into file volume with the given offset in x dimension.");
    tgtAssert(vol->getDimensions().y + offset.y <= dimensions_.y, "Brick does not fit into file volume with the given offset in y dimension.");
    tgtAssert(vol->getDimensions().z + offset.z <= dimensions_.z, "Brick does not fit into file volume with the given offset in z dimension.");

    const size_t numChannels = getNumberOfChannels();
    tgtAssert(firstChannel + numberOfChannels <= numChannels, "Invalid channel range.");
    tgtAssert(vol->getNumChannels() == numberOfChannels, "vol's number of channels must be exactly numberOfChannels.");
    tgtAssert(vol->getBaseType() == getBaseType(), "vol's base type does not match the file volume's");

    try {
        //Select the hyperslab (Brick + fourth dimension)

        //Just in case there are more than one channel: Allocate 4 dimensions.
        //If there is only one channel the 4. dimension will not do anything.
        hsize_t count[4];
        vec4TgtToHDF5(tgt::svec4(numberOfChannels, vol->getDimensions()), count);
        hsize_t start[4];
        vec4TgtToHDF5(tgt::svec4(firstChannel, offset), start);

        H5::DataSpace fileSpace(dataSet_->getSpace());
        fileSpace.selectHyperslab(H5S_SELECT_SET, count, start);

        //Memory space does not have an offset
        H5::DataSpace memSpace(numChannels == 1 ? 3 : 4, count);

        //Write the volume to disk.
        dataSet_->write(vol->getData(), dataSet_->getDataType(), memSpace, fileSpace);

    } catch(H5::Exception error) { // catch HDF5 exceptions
        LERROR(error.getFuncName() + ": " + error.getDetailMsg());
        throw tgt::IOException("An Error occured while writing volume to file " + getFileName());
    }
}

HDF5FileVolume::HDF5FileVolume(std::unique_ptr<H5::H5File> file, std::unique_ptr<H5::DataSet> dataSet, const std::string& fileName, const std::string& volumeLocation)
    : file_(std::move(file))
    , dataSet_(std::move(dataSet))
    , fileName_(fileName)
    , volumeLocation_(volumeLocation)
{
    dimensions_ = this->readDimensions();
}

tgt::svec3 HDF5FileVolume::readDimensions() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    H5::DataSpace dataSpace = dataSet_->getSpace();
    hsize_t nDim = dataSpace.getSimpleExtentNdims();
    // check if the data set is at least 3 dimensional
    if(nDim != 3 && nDim != 4) {
        throw tgt::IOException("HDF5 dataset is not 3 or 4 dimensional");
    }
    // Get the exact dimensions of the data set
    hsize_t dim[4];
    dataSpace.getSimpleExtentDims(dim);
    // Another channel (in the 4-dim case) would be the fastest changing channel => dim[3]).
    return vec3HDF5ToTgt(dim);
}

size_t HDF5FileVolume::getNumberOfChannels() const {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);

    H5::DataSpace dataSpace = dataSet_->getSpace();
    hsize_t nDim = dataSpace.getSimpleExtentNdims();
    if(nDim == 3) {
        // 3 dimensional => only one channel
        return 1;
    } else if(nDim == 4) {
        // 4 dimensional => We have to take a closer look at the dataSpace
        hsize_t dim[4];
        dataSpace.getSimpleExtentDims(dim);
        // Fastest changing dimension (channel) is the first.
        return vec4HDF5ToTgt(dim).x;
    } else {
        throw tgt::IOException("HDF5 dataset is not 3 or 4 dimensional");
    }
}
} // namespace voreen
