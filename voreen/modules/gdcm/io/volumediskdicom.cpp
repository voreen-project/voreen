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

#include "volumediskdicom.h"

#include "gdcmvolumereader.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/utils/hashing.h"

namespace voreen {

    const std::string VolumeDiskDicom::loggerCat_("voreen.tiff.VolumeDiskDicom");

VolumeDiskDicom::VolumeDiskDicom(const std::string& format, tgt::svec3 dimensions, const DicomInfo& info, const std::vector<std::string>& files)
    : VolumeDisk(format, dimensions)
    , info_(info)
    , sliceFiles_(files)
{ }

VolumeDiskDicom::~VolumeDiskDicom() { }


std::string VolumeDiskDicom::getHash() const {
    std::string configStr;

    // patient, study and series information
    configStr += info_.getPatientId() + "#";
    configStr += info_.getStudyInstanceUID() + "#";
    configStr += info_.getSeriesInstanceUID() + "#";
    //data format
    configStr += info_.getFormat() + "#";
    configStr += std::to_string(info_.getNumberOfFrames()) + "#";
    //file information
    for (size_t i = 0; i < sliceFiles_.size(); ++i) {
        configStr += sliceFiles_.at(i) + "#";
        configStr += genericToString(tgt::FileSystem::fileTime(sliceFiles_.at(i))) + "#";
        configStr += genericToString(tgt::FileSystem::fileSize(sliceFiles_.at(i))) + "#";
    }

    return VoreenHash::getHash(configStr);
}

VolumeRAM* VolumeDiskDicom::loadVolume() const {

    if (sliceFiles_.empty())
        throw tgt::Exception("VolumeDiskDicom: no files in list!");

    //check if the file is a multiframe DICOM file -> call the GdcmVolumeReader method because loadSlices does not allow that
    if (info_.getNumberOfFrames() > 1) {

        // create GdcmVolumeReader with progress bar
        ProgressBar* progress = VoreenApplication::app()->createProgressDialog();
        GdcmVolumeReader gdcmVolumeReader(progress);

        return gdcmVolumeReader.loadMultiframeDicomFile(info_, sliceFiles_);
    }

    //simply load all files
    return loadSlices(0, sliceFiles_.size() - 1);
}

VolumeRAM* VolumeDiskDicom::loadSlices(const size_t firstSlice, const size_t lastSlice) const {

    if (sliceFiles_.empty())
        throw tgt::Exception("VolumeDiskDicom: no files in list!");

    //check if the data set is a multiframe DICOM file -> throw exception (loading slices of a multiframe data set currently not supported)
    if (info_.getNumberOfFrames() > 1)
        throw tgt::Exception("VolumeDiskDicom: loading slices of a multiframe DICOM data set currently not supported!");

    tgtAssert(firstSlice <= lastSlice && firstSlice < sliceFiles_.size(), "invalid firstZSlice");
    tgtAssert(lastSlice < sliceFiles_.size(), "invalid lastZSlice");
    tgtAssert(!sliceFiles_.empty(), "dicom dataset has no files");

    //create vector with only the slices to be loaded
    std::vector<std::string> slicesToLoad(sliceFiles_.begin() + firstSlice, sliceFiles_.begin() + (lastSlice + 1));

    // create GdcmVolumeReader with progress bar
    ProgressBar* progress = VoreenApplication::app()->createProgressDialog();
    GdcmVolumeReader gdcmVolumeReader(progress);

    //call method of GdcmVolumeReader that gets file list and DicomInfo and return VolumeRAM
    return gdcmVolumeReader.loadDicomSlices(info_, slicesToLoad);
}

VolumeRAM* VolumeDiskDicom::loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
    throw tgt::Exception("VolumeDiskDicom: loading bricks is currently not supported!");
}

} // namespace voreen
