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

#include "hdf5volumewriter.h"
#include "../utils/hdf5utils.h"
#include "hdf5filevolume.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/io/progressbar.h"

#include <regex> //For temp hack

#include "tgt/filesystem.h"

namespace { // anonymous helper functions

/// Write the volume data given by the VolumeDisk representation to the specified dataSet slice by slice.
void writeToDataSet(const voreen::HDF5FileVolume& fileVol, const voreen::VolumeBase* volume, voreen::ProgressBar* progressBar) {
    tgtAssert(volume, "No Volume representation.");

    if(progressBar) {
        progressBar->setProgress(0.0f);
    }
    unsigned int numberOfSlices = volume->getDimensions().z;
    for(unsigned int z = 0; z<volume->getDimensions().z; ++z) {
        // Get the slice next to be written to the data set.
        const voreen::VolumeRAM* slice = volume->getSlice(z);

        // Write to the dataSet using the dataSpaces constructed and set up previously
        fileVol.writeSlices(slice, z);
        delete slice;

        if(progressBar) {
            progressBar->setProgress(static_cast<float>(z)/static_cast<float>(numberOfSlices));
        }
    }
}

// Temporary hack: read chunking and deflate value from file name.
void readParameters(const std::string& filename, tgt::svec3& chunkdim, int& deflateLevel, bool& shuffleEnabled) {
    std::smatch m;
    std::regex chunkRegex("_cnk(\\d+)x(\\d+)x(\\d+)");
    std::regex deflateRegex("_zip(\\d)");
    std::regex shuffleRegex("_shf");
    chunkdim = tgt::svec3::zero; //symbolic for maximum size
    deflateLevel = 1; //Deflating is relatively cheap (on level 1) and saves a lot of space for binary volumes
    shuffleEnabled = false;
    if(std::regex_search(filename, m, chunkRegex) && !m.empty()) {
        chunkdim.x = std::stoul(m[1]);
        chunkdim.y = std::stoul(m[2]);
        chunkdim.z = std::stoul(m[3]);
    }
    if(std::regex_search(filename, m, deflateRegex) && !m.empty()) {
        deflateLevel = std::stoi(m[1]);
    }

    // Set Shuffle if wanted.
    if(std::regex_search(filename, m, shuffleRegex) && !m.empty()) {
        shuffleEnabled = true;
    }
}

} // anonymous helper functions

// ---------------------------------------------------------------------------------------------------------------------

namespace voreen {

const std::string HDF5VolumeWriter::loggerCat_("voreen.hdf5.HDF5VolumeWriter");
const std::string HDF5VolumeWriter::VOLUME_DATASET_NAME("volume");

HDF5VolumeWriter::HDF5VolumeWriter(ProgressBar* progress)
    : VolumeWriter(progress)
{
    extensions_.push_back("h5");
}

void HDF5VolumeWriter::write(const std::string& filename, const VolumeBase* volumeHandle, const std::string& dataSetLocation, bool truncateFile, int deflateLevel, tgt::svec3 chunkDim, bool shuffleEnabled) {
    std::string hdf5name = tgt::FileSystem::cleanupPath(filename);

    // Create the file volume
    std::unique_ptr<HDF5FileVolume> fileVolume = HDF5FileVolume::createVolume(hdf5name, dataSetLocation, volumeHandle->getBaseType(), volumeHandle->getDimensions(), volumeHandle->getNumChannels(), truncateFile, deflateLevel, chunkDim, shuffleEnabled);

    // Set up progressbar
    if(getProgressBar()) {
        getProgressBar()->setTitle("Writing to HDF5 data set...");
        getProgressBar()->show();
    }

    // Write to the file volume using an appropriate representation
    writeToDataSet(*fileVolume, volumeHandle, getProgressBar());

    // We are done: Hide progressbar
    if(getProgressBar()) {
        getProgressBar()->hide();
    }

    // Add spacing information about the volume to the data set.
    fileVolume->writeSpacing(volumeHandle->getSpacing());

    // Add offset information about the volume to the data set.
    fileVolume->writeOffset(volumeHandle->getOffset());

    // Add the physical-to-world transformation of the volume to the data set.
    fileVolume->writePhysicalToWorldTransformation(volumeHandle->getPhysicalToWorldMatrix());

    // Write RealWorldMapping information about the volume to the data set, if available.
    if(volumeHandle->hasMetaData("RealWorldMapping")) {
        fileVolume->writeRealWorldMapping(volumeHandle->getRealWorldMapping());
    }

    // Write all available derived data
    fileVolume->writeDerivedData(volumeHandle);

    // Write meta data.
    fileVolume->writeMetaData(volumeHandle);
}
void HDF5VolumeWriter::write(const std::string& filename, const VolumeBase* volumeHandle) {
    tgtAssert(volumeHandle, "No volume");
    std::string hdf5name = tgt::FileSystem::cleanupPath(filename);
    LINFO("Saving " << hdf5name);

    tgt::svec3 chunkDim;
    int deflateLevel;
    bool shuffleEnabled;
    bool truncateFile = true;
    // Read parameters for file creation
    readParameters(hdf5name, chunkDim, deflateLevel, shuffleEnabled);

    // Write the volume using the more general method
    write(hdf5name, volumeHandle, VOLUME_DATASET_NAME, truncateFile, deflateLevel, chunkDim, shuffleEnabled);
}

VolumeWriter* HDF5VolumeWriter::create(ProgressBar* progress) const {
    return new HDF5VolumeWriter(progress);
}

} // namespace voreen
