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

#ifndef VRN_HDF5FILEVOLUME_H
#define VRN_HDF5FILEVOLUME_H

#include "modules/hdf5/utils/hdf5utils.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumepreview.h"

namespace voreen {

/**
 * Represents a volume inside a HDF5 file.
 *
 * Every 3 or 4 dimensional dataset in the file can be considered a volume
 * For a 4 dimensional dataset the additional fastest changing dimension
 * is considered the channel dimension so the layout in hdf5 style indexing
 * is: zyxc
 */
class HDF5FileVolume {
public:

    /**
     * Constant used for converting HDF5-spacing into Voreen-spacing anc vice-versa.
     */
    static const float MM_PER_HDF5_UNIT_OF_LENGTH;

    /**
     * Open a data set inside a hdf5 file as a HDF5FileVolume.
     * @throws tgt::IOException if the file or the dataset could not be opened.
     */
    static std::unique_ptr<HDF5FileVolume> openVolume(const std::string& fileName, const std::string& volumeLocation, bool readOnly = true);

    /**
     * Create a data set inside a hdf5 file as a HDF5FileVolume.
     * @param baseType The voreen style identifier for the volume base type.
     * @param truncateFile deletes an existing file if true. Otherwise the dataSet will be created within the existing file.
     *          This way a file can contain more than one volume.
     * @param chunkSize The size of the chunks the volume will be saved as in memory.
     *        The default is zero, which is symbolic for a xy-slice of the volume.
     *
     * @throws tgt::IOException if the file could not be created/opened or the dataset could not be created.
     *
     * @note Note that HDF5 does not provide means to delete or replace a dataset inside the file. Thus, if the caller does
     *          truncate the file and a dataset at volumeLocation already exists, an IOException will be thrown caused by
     *          a HDF5 error.
     */
    static std::unique_ptr<HDF5FileVolume> createVolume(const std::string& fileName, const std::string& volumeLocation, const std::string& baseType, tgt::svec3 dimensions, size_t numChannels = 1, bool truncateFile = true, int deflateLevel = 0, tgt::svec3 chunkSize=tgt::svec3::zero, bool shuffle = false);


    /**
     * Destructor.
     * Frees allocated HDF5 ressources.
     * @note Locks hdf5libMutex.
     */
    ~HDF5FileVolume();

    /**
     * Get the dimensions of this volume.
     */
    tgt::svec3 getDimensions() const;

    /**
     * Read the chunk size (i.e., a unit of organization (a block within the volume) in hdf5 files) of this volume.
     * @note Locks hdf5libMutex.
     */
    tgt::svec3 getChunkSize() const;

    /**
     * Tries to read the spacing of this volume.
     * @return the spacing of this volume or 0 if this information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    tgt::vec3* tryReadSpacing() const;

    /**
     * Write a spacing for this volume.
     * @param spacing The new spacing.
     * @note Locks hdf5libMutex.
     */
    void writeSpacing(const tgt::vec3& spacing) const;

    /**
     * Tries to read the offset of this volume.
     * @return the offset of this volume or 0 if this information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    tgt::vec3* tryReadOffset() const;

    /**
     * Write an offset for this volume.
     * @param offset The new offset.
     * @note Locks hdf5libMutex.
     */
    void writeOffset(const tgt::vec3& offset) const;

    /**
     * Tries to read the physical-to-world transformation of this volume.
     * @return The physical-to-world transformation of this volume or 0 if this information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    tgt::mat4* tryReadPhysicalToWorldTransformation() const;

    /**
     * Write a physical-to-world transformation for this volume.
     * @param transformation The new transformation.
     * @note Locks hdf5libMutex.
     */
    void writePhysicalToWorldTransformation(const tgt::mat4& transformation) const;

    /**
     * Try to read RealWorldMapping information.
     * @return A RealWorldMapping for this volume or 0 if that information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    RealWorldMapping* tryReadRealWorldMapping() const;

    /**
     * Write a real world mapping to this file volume.
     * @param rwm The new real world mapping.
     * @note Locks hdf5libMutex.
     */
    void writeRealWorldMapping(const RealWorldMapping& rwm) const;

    /**
     * Get the base type (float, uint8, ...) of this volume
     * @note Locks hdf5libMutex.
     */
    std::string getBaseType() const;

    /**
     * Try to read VolumeMinMax information.
     * @return A VolumeMinMax for this volume or 0 if that information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    VolumeMinMax* tryReadVolumeMinMax(size_t channel = 0) const;

    /**
     * Write the given VolumeMinMax to this file volume.
     * @param mm The VolumeMinMax.
     * @note mm must have the same number of channels as this file volume.
     * @note Locks hdf5libMutex.
     */
    void writeVolumeMinMax(const VolumeMinMax* mm) const;

    /**
     * Try to read VolumeHistogramIntensity information.
     * @return A VolumeHistogramIntensity for this volume or 0 if that information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    VolumeHistogramIntensity* tryReadVolumeHistogramIntensity(size_t channel = 0) const;

    /**
     * Write the given VolumeHistogramIntensity to this file volume.
     * @param hi The VolumeHistogramIntensity.
     * @note hi must have the same number of channels as this file volume.
     * @note Locks hdf5libMutex.
     */
    void writeVolumeHistogramIntensity(const VolumeHistogramIntensity* hi) const;

    /**
     * Try to read VolumeHistogramIntensityGradient information.
     * @return A VolumeHistogramIntensityGradient for this volume or 0 if that information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    VolumeHistogramIntensityGradient* tryReadVolumeHistogramIntensityGradient(size_t channel = 0) const;

    /**
     * Write the given VolumeHistogramIntensityGradient to this file volume.
     * @param hig The VolumeHistogramIntensityGradient.
     * @note hig must have the same number of channels as this file volume.
     * @note Locks hdf5libMutex.
     */
    void writeVolumeHistogramIntensityGradient(const VolumeHistogramIntensityGradient* hig) const;

    /**
     * Try to read VolumePreview information.
     * @return A VolumePreview for this volume or 0 if that information is not available.
     * @note The caller is responsible for deleting the returned object.
     * @note Locks hdf5libMutex.
     */
    VolumePreview* tryReadVolumePreview() const;

    /**
     * Write the given VolumePreview to this file volume.
     * @param preview The VolumePreview.
     * @note Locks hdf5libMutex.
     */
    void writeVolumePreview(const VolumePreview* preview) const;

    /**
     * Read all available derived data from this volume.
     * @return A collection of derived data from this volume.
     * @note The caller is responsible for deleting the VolumeDerivedData objects.
     * @note Locks hdf5libMutex.
     */
    std::vector<VolumeDerivedData*> readDerivedData(size_t channel = 0) const;

    /**
     * Write each DerivedData of vol to this file volume if applicable.
     * @param vol The volume this DerivedData will be read from.
     * @note vol must have the same number of channels as DerivedData.
     * @note Locks hdf5libMutex.
     */
    void writeDerivedData(const VolumeBase* vol) const;

    /**
     * Reads meta data from this volume, if existing.
     * @return A collection of meta data from this volume.
     * @note The caller is responsible for deleting the MetaDataBase objects.
     * @note Locks hdf5libMutex.
     */
    MetaDataContainer readMetaData() const;

    /**
     * Write each MetaData of vol to this file volume if applicable.
     * @param vol The volume the meta data will be read from.
     * @note Locks hdf5libMutex.
     */
    void writeMetaData(const VolumeBase* vol) const;

    /**
     * Get the (HDF5) file's name.
     */
    const std::string& getFileName() const;

    /**
     * Get the DataSet location that contains this volume in the HDF5 file.
     */
    const std::string& getVolumeLocation() const;

    /**
     * Get the number of channels the dataSet contains.
     * @note Locks hdf5libMutex.
     */
    size_t getNumberOfChannels() const;

    /**
     * Loads the channel/timestep from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @note Locks hdf5libMutex.
     */
    VolumeRAM* loadVolume(size_t channel = 0) const;

    /**
     * Write the data of vol to this file volume.
     * @param vol The volume to read from.
     *
     * @note vol and this file volume must have the same dimensions and number of channels.
     * @note vol's base type must match the one reported by getBaseType()
     * @note Locks hdf5libMutex.
     */
    void writeVolume(const VolumeRAM* vol) const;

    /**
     * Write the data of vol to a range of channels of this file volume.
     * @param vol The volume to read from.
     * @param firstChannel The first channel the volume will be written to.
     * @param numberOfChannels number of channels to write to.
     *
     * @note vol and this file volume must have the same dimensions.
     * @note vol must have exactly numberOfChannels channels.
     * @note firstChannel must be in [0, getNumberOfChannels()-numberOfChannels].
     * @note numberOfChannels must be in [1, getNumberOfChannels()].
     * @note vol's base type must match the one reported by getBaseType()
     * @note Locks hdf5libMutex.
     */
    void writeVolume(const VolumeRAM* vol, const size_t firstChannel, const size_t numberOfChannels = 1) const;

    /**
     * Loads a set of consecutive z slices of the HDF5 channel/timestep from disk
     * and returns them as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param firstZSlice first slice of the slice range to load (inclusive)
     * @param lastZSlice last slice of the slice range (inclusive)
     * @note Locks hdf5libMutex.
     *
     */
    VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice, size_t channel = 0) const;

    /**
     * Write the data of vol as slices to this file volume.
     * @param vol The volume containing the slices (slices in z-Dimension).
     * @param firstSlice z-position the first of the slices will be written to.
     *
     * @note vol and this file volume must have the same x and y dimensions and number of channels as this file volume.
     * @note The slice must fit into filevolume, i.e.: vol.dimensions.z+firstlice <= filevolume.dimensions.z
     * @note vol's base type must match the one reported by getBaseType()
     * @note Locks hdf5libMutex.
     */
    void writeSlices(const VolumeRAM* vol, const size_t firstSlice) const;

    /**
     * Write the data of vol as slices to a range channels of this file volume.
     * @param vol The volume containing the slices (slices in z-Dimension).
     * @param firstSlice z-position the first of the slices will be written to.
     * @param firstChannel The first channel the slices will be written to.
     * @param numberOfChannels number of channels to write to.
     *
     * @note vol and this file volume must have the same x and y dimensions.
     * @note The slice must fit into filevolume, i.e.: vol.dimensions.z+firstlice <= filevolume.dimensions.z
     * @note vol must have exactly numberOfChannels channels.
     * @note firstChannel must be in [0, getNumberOfChannels()-numberOfChannels].
     * @note numberOfChannels must be in [1, getNumberOfChannels()].
     * @note vol's base type must match the one reported by getBaseType()
     *
     * @note Locks hdf5libMutex.
     */
    void writeSlices(const VolumeRAM* vol, const size_t firstSlice, const size_t firstChannel, const size_t numberOfChannels = 1) const;

    /**
     * Loads a brick of the HDF5 channel/timestep volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param offset lower-left-front corner voxel of the brick to load
     * @param dimensions dimension of the brick to load
     * @note Locks hdf5libMutex.
     */
    VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions, size_t channel = 0) const;

    /**
     * Write the data of vol as a brick to this file volume.
     * @param vol The volume the brick will be read from.
     * @param offset The offset inside this volume for writing the brick.
     *
     * @note vol and this file volume must have the same number of channels.
     * @note The brick must fit into filevolume, i.e.:
     *      vol.dimensions.d+offset.d <= filevolume.dimensions.d for all dimensions d.
     * @note vol's base type must match the one reported by getBaseType()
     * @note Locks hdf5libMutex.
     */
    void writeBrick(const VolumeRAM* vol, const tgt::svec3& offset) const;

    /**
     * Write the data of vol as a brick to a range of channels of this file volume.
     * @param vol The volume the brick will be read from.
     * @param offset The offset inside this volume for writing the brick.
     * @param firstChannel The first channel the brick will be written to.
     * @param numberOfChannels number of channels to write to.
     *
     * @note The brick must fit into filevolume, i.e.:
     *      vol.dimensions.d+offset.d <= filevolume.dimensions.d for all dimensions d.
     * @note vol must have exactly numberOfChannels channels.
     * @note firstChannel must be in [0, getNumberOfChannels()-numberOfChannels].
     * @note numberOfChannels must be in [1, getNumberOfChannels()].
     * @note vol's base type must match the one reported by getBaseType()
     * @note Locks hdf5libMutex.
     */
    void writeBrick(const VolumeRAM* vol, const tgt::svec3& offset, const size_t firstChannel, const size_t numberOfChannels = 1) const;

private:

    /**
     * Constructor.
     * @param fileName Path to the HDF5 file that contains the volume.
     * @param volumeLocation Location of the DataSet in the HDF5 file that contains the volume.
     * @param channel Channel (i.e. coordinate of the fourth dimension) of the DataSet this volume represents.
     * @note Locks hdf5libMutex.
     */
    HDF5FileVolume(std::unique_ptr<H5::H5File> file, std::unique_ptr<H5::DataSet> dataSet, const std::string& fileName, const std::string& volumeLocation);

    /**
     * Disable copy constructor.
     */
    HDF5FileVolume(const HDF5FileVolume&);


    /**
     * Read the dimensions of this volume from file.
     * @note Locks hdf5libMutex.
     */
    tgt::svec3 readDimensions() const;

    /// The HDF5-library file object
    std::unique_ptr<H5::H5File> file_;
    /// The HDF5-library DataSet object. Contains the volume.
    std::unique_ptr<H5::DataSet> dataSet_;
    /// The file's name.
    std::string fileName_;
    /// The dataSet's location.
    std::string volumeLocation_;
    /// Dimensions are cashed in ram for faster access:
    tgt::svec3 dimensions_;

    static const std::string SPACING_ATTRIBUTE_NAME;
    static const std::string OFFSET_ATTRIBUTE_NAME;
    static const std::string PHYSICAL_TO_WORLD_TRANSFORMATION_ATTRIBUTE_NAME;
    static const std::string METADATA_STRING_NAME;
    static const std::string REALWORLDMAPPING_SCALE_ATTRIBUTE_NAME;
    static const std::string REALWORLDMAPPING_OFFSET_ATTRIBUTE_NAME;
    static const std::string REALWORLDMAPPING_UNIT_ATTRIBUTE_NAME;
    static const std::string REPRESENTATION_MINMAX_ATTRIBUTE_NAME;
    static const std::string REPRESENTATION_HISTOGRAMINTENSITY_ATTRIBUTE_NAME;
    static const std::string REPRESENTATION_HISTOGRAMINTENSITY_METADATA_ATTRIBUTE_NAME;
    static const std::string REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_ATTRIBUTE_NAME;
    static const std::string REPRESENTATION_HISTOGRAMINTENSITYGRADIENT_METADATA_ATTRIBUTE_NAME;
    static const std::string REPRESENTATION_PREVIEW_ATTRIBUTE_NAME;
    static const std::string REPRESENTATION_PREVIEW_METADATA_ATTRIBUTE_NAME;

    static const std::string loggerCat_;
};

} // namespace voreen
#endif //VRN_HDF5FILEVOLUME_H
