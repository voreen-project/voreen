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

#include "volumediskometiff.h"

#include "ometiffvolumereader.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/utils/hashing.h"

namespace voreen {

OMETiffFile::OMETiffFile(const std::string& filename, size_t numDirectories, size_t firstZ, size_t firstT, size_t firstC)
    : filename_(filename)
    , numDirectories_(numDirectories)
    , firstZ_(firstZ)
    , firstT_(firstT)
    , firstC_(firstC)
{}

OMETiffFile::OMETiffFile()
    : filename_("")
    , numDirectories_(0)
    , firstZ_(0)
    , firstT_(0)
    , firstC_(0)
{}

std::string OMETiffFile::toString() const {
    std::string result = "OMETiffFile[";
    result += "filename=" + filename_ + ", ";
    result += "numDirectories=" + itos(numDirectories_) + ", ";
    result += "firstZ=" + itos(firstZ_) + ", ";
    result += "firstT=" + itos(firstT_) + ", ";
    result += "firstC=" + itos(firstC_) + "]";
    return result;
}

//-------------------------------------------------------------------------------------------------

OMETiffStack::OMETiffStack(const std::vector<OMETiffFile>& files, const std::string& dimensionOrder,
    const std::string& datatype, tgt::svec3 volumeDim, tgt::vec3 voxelSpacing,
    size_t sizeC, size_t sizeT, size_t numSlices)
    : files_(files)
    , dimensionOrder_(dimensionOrder)
    , datatype_(datatype)
    , volumeDim_(volumeDim)
    , voxelSpacing_(voxelSpacing)
    , sizeC_(sizeC)
    , sizeT_(sizeT)
    , numSlices_(numSlices)
{
    tgtAssert(!files_.empty(), "no files passed");
    tgtAssert(!dimensionOrder_.empty(), "no dimension order passed");
    tgtAssert(!datatype.empty(), "no dimension order passed");
    tgtAssert(tgt::hand(tgt::greaterThan(voxelSpacing, tgt::vec3::zero)), "invalid spacing");
    tgtAssert(sizeC_ > 0, "sizeC must be greater 0");
    tgtAssert(sizeT_ > 0, "sizeT must be greater 0");
    tgtAssert(numSlices_ > 0, "num slices must be greater 0")
}

OMETiffStack::OMETiffStack()
    : files_()
    , dimensionOrder_("")
    , datatype_("")
    , volumeDim_(0,0,0)
    , voxelSpacing_(0.f, 0.f, 0.f)
    , sizeC_(0)
    , sizeT_(0)
    , numSlices_(0)
{}

//-------------------------------------------------------------------------------------------------

const std::string VolumeDiskOmeTiff::loggerCat_("voreen.tiff.VolumeDiskOmeTiff");

VolumeDiskOmeTiff::VolumeDiskOmeTiff(const std::string& format, tgt::svec3 dimensions,
        const OMETiffStack& datastack, size_t channel, size_t timestep)
    : VolumeDisk(format, dimensions)
    , datastack_(datastack)
    , channel_(channel)
    , timestep_(timestep)
{
    tgtAssert(!datastack.files_.empty(), "no files passed");
    tgtAssert(!datastack.dimensionOrder_.empty(), "no dimension order passed");
    tgtAssert(!datastack.datatype_.empty(), "no dimension order passed");
    tgtAssert(tgt::hand(tgt::greaterThan(datastack.voxelSpacing_, tgt::vec3::zero)), "invalid spacing");
    tgtAssert(datastack.sizeC_ > 0, "sizeC must be greater 0");
    tgtAssert(datastack.sizeT_ > 0, "sizeT must be greater 0");
    tgtAssert(datastack.numSlices_ > 0, "num slices must be greater 0")

    tgtAssert(channel < datastack.sizeC_, "invalid channel");
    tgtAssert(timestep < datastack.sizeT_, "invalid timestep");
}

VolumeDiskOmeTiff::~VolumeDiskOmeTiff() {
}

const OMETiffStack& VolumeDiskOmeTiff::getDatastack() const {
    return datastack_;
}

std::string VolumeDiskOmeTiff::getHash() const {
    std::string configStr;

    // datastack properties
    for (size_t i=0; i<datastack_.files_.size(); i++) {
        configStr += datastack_.files_.at(i).filename_ + "#";
        configStr += genericToString(tgt::FileSystem::fileTime(datastack_.files_.at(i).filename_)) + "#";
        configStr += genericToString(tgt::FileSystem::fileSize(datastack_.files_.at(i).filename_)) + "#";
    }
    configStr += datastack_.datatype_ + "#";
    configStr += datastack_.dimensionOrder_ + "#";
    configStr += itos(datastack_.sizeC_) + "#";
    configStr += itos(datastack_.sizeT_) + "#";
    configStr += itos(datastack_.numSlices_) + "#";
    configStr += genericToString(datastack_.voxelSpacing_) + "#";
    configStr += genericToString(datastack_.volumeDim_) + "#";

    // channel/timestep
    configStr += itos(channel_) + "#";
    configStr += itos(timestep_) + "#";

    return VoreenHash::getHash(configStr);
}

VolumeRAM* VolumeDiskOmeTiff::loadVolume() const {
    tgtAssert(!datastack_.files_.empty(), "ometiff datastack has no files");
    LDEBUG("Loading OmeTiff volume: filename=" << datastack_.files_.front().filename_ << ", channel=" << channel_ << ", timestep=" << timestep_);

    // create OME Tiff reader with progress bar
    ProgressBar* progress = VoreenApplication::app()->createProgressDialog();
    OMETiffVolumeReader omeTiffReader(progress);

    // load volume
    std::vector<VolumeRAM*> ramVolumes;
    ramVolumes = omeTiffReader.loadVolumesIntoRam(datastack_, static_cast<int>(channel_), static_cast<int>(timestep_));
    tgtAssert(!ramVolumes.empty(), "no ram volume returned (exception expected)");
    if (ramVolumes.size() > 1) {
        LWARNING("OmeTiff volume reader returned more than one volume. Discarding surplus volumes.");
        for (size_t i=1; i<ramVolumes.size(); i++)
            delete ramVolumes.at(i);
    }

    delete progress;

    return ramVolumes.front();
}

VolumeRAM* VolumeDiskOmeTiff::loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
    tgtAssert(firstZSlice <= lastZSlice && firstZSlice < datastack_.volumeDim_.z, "invalid firstZSlice");
    tgtAssert(lastZSlice < datastack_.volumeDim_.z, "invalid lastZSlice");
    tgtAssert(!datastack_.files_.empty(), "ometiff datastack has no files");

    LDEBUG("Loading OmeTiff slices: filename=" << datastack_.files_.front().filename_ << ", channel=" << channel_ << ", timestep=" << timestep_
        << ", firstZSlice=" << firstZSlice << ", lastZSlice=" << lastZSlice);

    // create OME Tiff reader without progress bar
    OMETiffVolumeReader omeTiffReader(0);

    // load volume
    std::vector<VolumeRAM*> ramVolumes;
    ramVolumes = omeTiffReader.loadVolumesIntoRam(datastack_, static_cast<int>(channel_), static_cast<int>(timestep_),
        tgt::ivec3(0, 0, static_cast<int>(firstZSlice)),
        tgt::ivec3(static_cast<int>(getDimensions().x-1), static_cast<int>(getDimensions().y-1), static_cast<int>(lastZSlice)));
    tgtAssert(!ramVolumes.empty(), "no ram volume returned (exception expected)");
    if (ramVolumes.size() > 1) {
        LWARNING("OmeTiff volume reader returned more than one volume. Discarding surplus volumes.");
        for (size_t i=1; i<ramVolumes.size(); i++)
            delete ramVolumes.at(i);
    }

    return ramVolumes.front();
}

VolumeRAM* VolumeDiskOmeTiff::loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
    tgt::ivec3 llf = tgt::ivec3(offset);
    tgt::ivec3 urb = llf + tgt::ivec3(dimensions) - tgt::ivec3(1);

    tgtAssert(tgt::hand(tgt::lessThan(llf, tgt::ivec3(datastack_.volumeDim_))), "invalid offset");
    tgtAssert(tgt::hand(tgt::lessThan(urb, tgt::ivec3(datastack_.volumeDim_))), "invalid dimensions");
    tgtAssert(!datastack_.files_.empty(), "ometiff datastack has no files");

    LDEBUG("Loading OmeTiff brick: filename=" << datastack_.files_.front().filename_ << ", channel=" << channel_ << ", timestep=" << timestep_
        << ", llf=" << llf << ", urb=" << urb);

    // create OME Tiff reader without progress bar
    OMETiffVolumeReader omeTiffReader(0);

    // load volume
    std::vector<VolumeRAM*> ramVolumes;
    ramVolumes = omeTiffReader.loadVolumesIntoRam(datastack_, static_cast<int>(channel_), static_cast<int>(timestep_), llf, urb);
    tgtAssert(!ramVolumes.empty(), "no ram volume returned (exception expected)");
    if (ramVolumes.size() > 1) {
        LWARNING("OmeTiff volume reader returned more than one volume. Discarding surplus volumes.");
        for (size_t i=1; i<ramVolumes.size(); i++)
            delete ramVolumes.at(i);
    }

    return ramVolumes.front();
}

} // namespace voreen
