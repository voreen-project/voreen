/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "ometiffvolumereader.h"

#include "volumediskometiff.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/properties/boolproperty.h"

#include <tinyxml/tinyxml.h>

#include <fstream>
#include <iostream>
#include <assert.h>

#include "tgt/exception.h"
#include "tgt/vector.h"
#include "tgt/texture.h"
#include "tgt/filesystem.h"

using std::string;
using tgt::vec3;
using tgt::ivec3;
using tgt::svec3;
using tgt::Texture;

namespace { // anonymous helper functions

void raiseIOException(const std::string& msg, const std::string& filename, voreen::ProgressBar* progress = 0) {
    LERRORC("voreen.ome.OmeTiffVolumeReader", msg + ": " + filename);
    if (progress)
        progress->hide();
    throw tgt::IOException(msg, filename);
}

void deleteVolumes(const std::vector<std::vector<voreen::VolumeRAM*> >& volumes) {
    for (size_t i=0; i<volumes.size(); i++) {
        for (size_t j=0; j<volumes[i].size(); j++) {
            delete volumes[i][j];
        }
    }
}

} // namespace anonymous

// ------------------------------------------------------------------------------------------------

namespace voreen {

const std::string OMETiffVolumeReader::loggerCat_ = "voreen.ome.OmeTiffVolumeReader";

OMETiffVolumeReader::OMETiffVolumeReader(ProgressBar* progress) : VolumeReader(progress)
{
    extensions_.push_back("ome.tiff");
    extensions_.push_back("ome.tif");

    protocols_.push_back("ome-tiff");
}

VolumeReader* OMETiffVolumeReader::create(ProgressBar* progress) const {
    return new OMETiffVolumeReader(progress);
}

VolumeList* OMETiffVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    const std::string masterFileName = origin.getPath();

    // determine selected channel/timestep
    int requestedChannel = -1;
    if (origin.getSearchParameter("channel") != "")
        requestedChannel = stoi(origin.getSearchParameter("channel"));
    int requestedTimestep = -1;
    if (origin.getSearchParameter("timestep") != "")
        requestedTimestep = stoi(origin.getSearchParameter("timestep"));

    // open master ome tiff file
    TIFF* tiffFile = TIFFOpen(masterFileName.c_str(), "r");
    if (!tiffFile)
        raiseIOException("Failed to open file", masterFileName);

    LINFO("Reading OME XML from master file: " << masterFileName);

    if (getProgressBar()) {
        getProgressBar()->setTitle("Opening OME-TIFF datastack");
        getProgressBar()->setProgressMessage("Scanning files...");
        getProgressBar()->show();
        getProgressBar()->forceUpdate();
    }

    // extract meta data from OME XML, stored in the image description field of the tiff file
    OMETiffStack stack;
    try {
        stack = extractStackInformation(tiffFile, tgt::FileSystem::dirName(masterFileName));
        TIFFClose(tiffFile);
        tiffFile = 0;
    }
    catch (tgt::Exception& e) {
        TIFFClose(tiffFile);
        raiseIOException(e.what(), masterFileName, getProgressBar());
    }

    // determine number of directories stored in each tiff file (requires to open files)
    try {
        determineDirectoryCount(stack);
    }
    catch (tgt::Exception& e) {
        raiseIOException(e.what(), masterFileName, getProgressBar());
    }

    tgtAssert(!stack.files_.empty(), "no files passed");
    tgtAssert(!stack.dimensionOrder_.empty(), "no dimension order passed");
    tgtAssert(!stack.datatype_.empty(), "no dimension order passed");
    tgtAssert(tgt::hand(tgt::greaterThan(stack.voxelSpacing_, tgt::vec3::zero)), "invalid spacing");
    tgtAssert(stack.sizeC_ > 0, "sizeC must be greater 0");
    tgtAssert(stack.sizeT_ > 0, "sizeT must be greater 0");
    tgtAssert(stack.numSlices_ > 0, "num slices must be greater 0")

    // log extracted data
    LINFO("DimensionOrder: " << stack.dimensionOrder_);
    LINFO("Volume dimensions: " << stack.volumeDim_);
    LINFO("Voxel spacing: " << stack.voxelSpacing_ * 1000.f << " micron");
    LINFO("Data type: " << stack.datatype_);
    LINFO("Num channels: " << stack.sizeC_);
    LINFO("Num timesteps: " << stack.sizeT_);
    if (requestedChannel > -1 && requestedTimestep > -1)
        LINFO("Selected channel/timestep: " << requestedChannel << "/" << requestedTimestep);
    else if (requestedChannel > -1)
        LINFO("Selected channel: " << requestedChannel);
    else if (requestedTimestep > -1)
        LINFO("Selected timestep: " << requestedTimestep);
    /*else
        LINFO("Loading " << stack.sizeC_ << " channels with " << stack.sizeT_ << " timesteps each"); */

    /*LDEBUG("TIFF files:");
    for (size_t i=0; i<stack.files_.size(); i++)
        LDEBUG(stack.files_.at(i).toString()); */

    // check requested channel and timestep
    if (requestedChannel >= static_cast<int>(stack.sizeC_)) {
        raiseIOException("Requested channel (" + itos(requestedChannel) + ") is greater/equal than number of channels in stack (" + itos(stack.sizeC_) + ")",
            masterFileName, getProgressBar());
    }
    if (requestedTimestep >= static_cast<int>(stack.sizeT_)) {
        raiseIOException("Requested timestep (" + itos(requestedTimestep) + ") is greater/equal than number of timesteps in stack (" + itos(stack.sizeT_) + ")",
            masterFileName, getProgressBar());
    }

    // create disk volume for each channel/timestep
    VolumeList* volumeList = new VolumeList();
    for (size_t c=0; c<stack.sizeC_; c++) {
        for (size_t t=0; t<stack.sizeT_; t++) {
            if ((c == requestedChannel || requestedChannel == -1) && (t == requestedTimestep || requestedTimestep == -1)) {
                VolumeDiskOmeTiff* diskRep = new VolumeDiskOmeTiff(stack.datatype_, stack.volumeDim_, stack, c, t);
                tgt::vec3 offset = -(stack.voxelSpacing_*(tgt::vec3)stack.volumeDim_) / 2.f;
                Volume* volumeHandle = new Volume(diskRep, stack.voxelSpacing_, offset);
                volumeHandle->setTimestep(static_cast<float>(t));
                volumeHandle->setMetaDataValue<IntMetaData, int>("Channel", (int)c);

                VolumeURL origin("ome-tiff", masterFileName);
                origin.addSearchParameter("channel", itos(c));
                origin.addSearchParameter("timestep", itos(t));
                volumeHandle->setOrigin(origin);

                volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));

                volumeList->add(volumeHandle);
            }
        }
    }
    tgtAssert(!volumeList->empty(), "volume list is empty");

    if (getProgressBar())
        getProgressBar()->hide();

    return volumeList;
}

VolumeBase* OMETiffVolumeReader::read(const VolumeURL& origin) {
    VolumeList* volumeList = read(origin.getURL());
    tgtAssert(!volumeList->empty(), "volume list is empty");

    VolumeBase* volume = volumeList->first();
    tgtAssert(volume, "volume is null");

    if (volumeList->size() > 1) {
        LWARNING("read(origin): more than one volume loaded. Discarding redundant volumes.");
        for (size_t i=1; i<volumeList->size(); i++)
            delete volumeList->at(i);
    }

    return volume;
}

std::vector<VolumeURL> OMETiffVolumeReader::listVolumes(const std::string& urlStr) const {
    VolumeURL url(urlStr);
    std::string filepath = url.getPath();

    TIFF* tiffFile = TIFFOpen(filepath.c_str(), "r");
    if (!tiffFile)
        raiseIOException("Failed to open file", filepath);

    // extract meta data from OME XML, stored in the image description field of the tiff file
    OMETiffStack stack;
    try {
        stack = extractStackInformation(tiffFile, filepath);
    }
    catch (tgt::Exception& e) {
        raiseIOException(e.what(), urlStr);
    }

    // identify volumes by channel and timestep
    std::vector<VolumeURL> volumeURLs;
    for (size_t c=0; c < static_cast<size_t>(stack.sizeC_); c++) {
        for (size_t t=0; t < static_cast<size_t>(stack.sizeT_); t++) {
            VolumeURL subUrl("ome-tiff", filepath);
            subUrl.addSearchParameter("channel", itos(c));
            subUrl.addSearchParameter("timestep", itos(t));
            subUrl.getMetaDataContainer().addMetaData("Channel", new IntMetaData((int)c));
            subUrl.getMetaDataContainer().addMetaData("Timestep", new IntMetaData((int)t));
            subUrl.getMetaDataContainer().addMetaData("Volume Dimensions", new IVec3MetaData(static_cast<ivec3>(stack.volumeDim_)));
            subUrl.getMetaDataContainer().addMetaData("Voxel Spacing", new StringMetaData(genericToString(stack.voxelSpacing_ * 1000.f) + " micron"));
            volumeURLs.push_back(subUrl);
        }
    }

    return volumeURLs;
}

std::vector<VolumeRAM*> OMETiffVolumeReader::loadVolumesIntoRam(const OMETiffStack& stack,
    int requestedChannel /*= -1*/, int requestedTimestep /*= -1*/,
    tgt::ivec3 llf /*= tgt::ivec3(-1)*/, tgt::ivec3 urb /*= tgt::ivec3(-1)*/) const {
    // check parameters
    if (requestedChannel >= (int)stack.sizeC_)
        throw std::invalid_argument("Invalid channel: " + itos(requestedChannel));
    if (requestedTimestep >= (int)stack.sizeT_)
        throw std::invalid_argument("Invalid timestep: " + itos(requestedTimestep));
    if (tgt::hor(tgt::greaterThanEqual(llf, tgt::ivec3(stack.volumeDim_))))
        throw std::invalid_argument("Invalid llf: " + genericToString(llf));
    if (tgt::hor(tgt::greaterThanEqual(urb, tgt::ivec3(stack.volumeDim_))))
        throw std::invalid_argument("Invalid urb: " + genericToString(urb));
    if (tgt::hor(tgt::greaterThan(llf, urb)))
        throw std::invalid_argument("Invalid llf/urb combination: llf=" + genericToString(llf) + ", urb=" + genericToString(urb));

    // stack parameter is created internally, so assertions should suffice
    tgtAssert(!stack.files_.empty(), "no files passed");
    tgtAssert(!stack.dimensionOrder_.empty(), "no dimension order passed");
    tgtAssert(!stack.datatype_.empty(), "no dimension order passed");
    tgtAssert(tgt::hand(tgt::greaterThan(stack.voxelSpacing_, tgt::vec3::zero)), "invalid spacing");
    tgtAssert(stack.sizeC_ > 0, "sizeC must be greater 0");
    tgtAssert(stack.sizeT_ > 0, "sizeT must be greater 0");
    tgtAssert(stack.numSlices_ > 0, "num slices must be greater 0")

    // either both or none of the llf/urb parameters must be given
    if (tgt::hor(tgt::lessThan(llf, tgt::ivec3::zero)) || tgt::hor(tgt::lessThan(urb, tgt::ivec3::zero))) {
        llf = tgt::ivec3(-1);
        urb = tgt::ivec3(-1);
    }

    // do not show progress bar for a single sub-volume (slice or brick)
    bool showProgress = requestedChannel < 0 || requestedTimestep < 0 || tgt::hor(tgt::lessThan(llf, tgt::ivec3(0)));
    if (showProgress && getProgressBar()) {
        getProgressBar()->setTitle("Loading OME-TIFF data set");
        getProgressBar()->show();
        getProgressBar()->forceUpdate();
    }

    // determine dimensions of output volumes: stack volume dim, if no slice range has been selected, or slice range otherwise
    const tgt::svec3 outputVolDim(llf.z >= 0 ? tgt::svec3(urb-llf+1) : stack.volumeDim_);

    // create two nested vectors of output volumes: first coordinate=channel, second coordinate=time
    // (volumes will be created on demand)
    std::vector<std::vector<VolumeRAM*> > volumes;
    for (size_t c=0; c < static_cast<size_t>(stack.sizeC_); c++) {
        std::vector<VolumeRAM*> timeSeries(stack.sizeT_, reinterpret_cast<VolumeRAM*>(0));
        volumes.push_back(timeSeries);
    }
    tgtAssert(volumes.size() == stack.sizeC_, "invalid size of volume vector");

    // insert current and max z,c,t values into helper vectors according to DimensionOrder
    int curZ = 0;
    int curC = 0;
    int curT = 0;
    std::vector<int*> curStackIndices = std::vector<int*>(3, reinterpret_cast<int*>(0));
    std::vector<int> dimSizes = std::vector<int>(3, 0);
    std::string::size_type posZ = stack.dimensionOrder_.find("Z");
    std::string::size_type posC = stack.dimensionOrder_.find("C");
    std::string::size_type posT = stack.dimensionOrder_.find("T");
    if (posZ < 2 || posZ > 4 || posC < 2 || posC > 4 || posT < 2 || posT > 4)
        throw tgt::Exception("Invalid DimensionOrder: " + stack.dimensionOrder_);
    curStackIndices.at(posZ-2) = &curZ;
    dimSizes.at(posZ-2) = static_cast<int>(stack.volumeDim_.z);
    curStackIndices.at(posC-2) = &curC;
    dimSizes.at(posC-2) = static_cast<int>(stack.sizeC_);
    curStackIndices.at(posT-2) = &curT;
    dimSizes.at(posT-2) = static_cast<int>(stack.sizeT_);
    tgtAssert(curStackIndices.at(0) != 0 && dimSizes.at(0) != 0 && curStackIndices.at(1) != 0 && dimSizes.at(1) != 0 && curStackIndices.at(2) != 0 && dimSizes.at(2) != 0,
        "dimension vector not properly initialized");


    //
    // Iterate over TIFF files and copy pixel content to respective positions in the volume stack
    //
    VolumeFactory volumeFac;
    size_t curSlice = 0;
    for (int fileID = 0; fileID < static_cast<int>(stack.files_.size()); fileID++) {
        const OMETiffFile& curFile = stack.files_.at(fileID);
        if (showProgress && getProgressBar())
            getProgressBar()->setProgressMessage("Loading " + curFile.filename_ + " ...");

        // current tiff file (open on demand)
        TIFF* curTiffFile = 0;

        // check firstZ, firstC, firstT parameters of current file against current coordinates
        if (curFile.firstZ_ != curZ || curFile.firstC_ != curC || curFile.firstT_ != curT) {
            LWARNING("First Z/T/C values of Tiff file '" << curFile.filename_ << "' do not match expected values");
            //deleteVolumes(volumes);
            //raiseIOException("First Z/T/C values of Tiff file do not match expected values", curFile.filename_, getProgressBar());
        }

        //
        // iterate over directories of current TIFF file
        //
        for (size_t tiffDir = 0; tiffDir < stack.files_.at(fileID).numDirectories_; tiffDir++) {
            tgtAssert(curC < (int)volumes.size(), "current C value larger than volumes vector");
            tgtAssert(curT < (int)volumes[curC].size(), "current T value larger than volumes vector");

            if (showProgress && getProgressBar()) {
                getProgressBar()->setProgress(static_cast<float>(curSlice) / static_cast<float>(stack.numSlices_-1));
                getProgressBar()->forceUpdate();
            }

            // ignore slice, if a specific channel/timestep/slice is requested, which does not match current channel/timestep/slice
            bool skipSlice = (requestedChannel >= 0  && curC != requestedChannel)  ||
                             (requestedTimestep >= 0 && curT != requestedTimestep) ||
                             (llf.z >= 0 && (llf.z > curZ || urb.z < curZ));
            if (!skipSlice) {

                // open current TIFF file, if not already happened
                if (!curTiffFile) {
                    curTiffFile = TIFFOpen(curFile.filename_.c_str(), "r");
                    if (!curTiffFile) {
                        deleteVolumes(volumes);
                        raiseIOException("Failed to open TIFF file", curFile.filename_, getProgressBar());
                    }
                }
                tgtAssert(curTiffFile, "cur tiff file not opened");

                // set current TIFF directory
                TIFFSetDirectory(curTiffFile,static_cast<tdir_t>(tiffDir));

                // retrieve current volume and create it, if not created yet
                VolumeRAM*& currentVolume = volumes[curC][curT];
                if (!currentVolume) {
                    try {
                        currentVolume = volumeFac.create(stack.datatype_, outputVolDim);
                    }
                    catch (std::exception& e) {
                        LERROR(e.what());
                        deleteVolumes(volumes);
                        if (curTiffFile)
                            TIFFClose(curTiffFile);
                        if (getProgressBar())
                            getProgressBar()->hide();
                        throw e;
                    }
                }
                tgtAssert(currentVolume, "current volume is null");

                // read current directory/slice
                try {
                    tgtAssert(llf.z == -1 || llf.z <= curZ, "invalid curZ (should have skipped this slice)");
                    const size_t bytesPerVoxel = currentVolume->getBytesPerVoxel();
                    const size_t zSliceInOutputVolume = llf.z >= 0 ? curZ-llf.z : curZ;
                    const size_t curSliceByteOffset = tgt::hmul(outputVolDim.xy()) * zSliceInOutputVolume * bytesPerVoxel;
                    char* curSlicePointer = reinterpret_cast<char*>(currentVolume->getData()) + curSliceByteOffset;

                    // should the entire slice be written to the output volume or has a brick been requested?
                    bool brickMode = tgt::hor(tgt::greaterThan(llf.xy(), tgt::ivec2(0)));
                    brickMode |=     urb.x > -1 && tgt::hor(tgt::lessThan(urb.xy(), tgt::ivec2(stack.volumeDim_.xy())-1));
                    if (brickMode) { // brick only => read entire tiff slice into temp buffer and extract sub-slice
                        tgt::ivec2 subsliceDim = urb.xy()-llf.xy() + tgt::ivec2(1);
                        tgtAssert(tgt::hand(tgt::greaterThan(subsliceDim, tgt::ivec2(0))) &&
                                  tgt::hand(tgt::lessThan(subsliceDim, tgt::ivec2(stack.volumeDim_.xy()))), "invalid subsliceDim");
                        size_t subSliceByteSize = tgt::hmul(subsliceDim) * bytesPerVoxel;

                        // read full slice
                        size_t sliceBytesize = tgt::hmul(stack.volumeDim_.xy()) * bytesPerVoxel;
                        char* sliceBuffer = new char[sliceBytesize];
                        readTiffDirectory(curTiffFile, stack.datatype_, stack.volumeDim_.xy(), sliceBuffer);

                        // copy sub-slice to output volume
                        size_t sliceBufferOffset = (llf.y * stack.volumeDim_.x + llf.x) * bytesPerVoxel;
                        size_t subSliceBufferOffset = 0;
                        for (size_t y=0; y<(size_t)subsliceDim.y; y++) {
                            tgtAssert((sliceBufferOffset + subsliceDim.x*bytesPerVoxel) <= sliceBytesize, "invalid sliceBufferOffset");
                            tgtAssert((subSliceBufferOffset + subsliceDim.x*bytesPerVoxel) <= subSliceByteSize, "invalid subSliceBufferOffset");
                            memcpy(curSlicePointer + subSliceBufferOffset, sliceBuffer + sliceBufferOffset, subsliceDim.x*bytesPerVoxel);
                            // advance one line in full slice and sub slice buffer
                            sliceBufferOffset    += stack.volumeDim_.x * bytesPerVoxel;
                            subSliceBufferOffset += subsliceDim.x * bytesPerVoxel;
                        }

                        delete[] sliceBuffer;
                    }
                    else { // full slice => write tiff slice directly into output volume
                        readTiffDirectory(curTiffFile, stack.datatype_, stack.volumeDim_.xy(), curSlicePointer);
                    }
                }
                catch (tgt::Exception& e) {
                    deleteVolumes(volumes);
                    if (curTiffFile)
                        TIFFClose(curTiffFile);
                    raiseIOException("Failed to read tiff slice: " + std::string(e.what()), curFile.filename_, getProgressBar());
                }
            }

            // update stack indices (pointers to curZ, curC, curT)
            curSlice++;
            (*curStackIndices[0])++; //< first index changes fastest
            if (*curStackIndices[0] >= dimSizes[0]) { //< handle overflow
                *curStackIndices[0] = 0;
                (*curStackIndices[1])++;
                if (*curStackIndices[1] >= dimSizes[1]) {
                    *curStackIndices[1] = 0;
                    (*curStackIndices[2])++;
                }
            }
            tgtAssert((*curStackIndices[0] < dimSizes[0] && *curStackIndices[1] < dimSizes[1] && *curStackIndices[2] < dimSizes[2]) ||
                (curSlice == stack.numSlices_), "invalid stack indices");

        } // directory iteration

        // close current TIFF file, if open
        if (curTiffFile) {
            TIFFClose(curTiffFile);
            curTiffFile = 0;
        }

    } // file iteration

    // collect created VolumeRAMs in result vector
    std::vector<VolumeRAM*> result;
    for (size_t c=0; c<volumes.size(); c++) {
        for (size_t t=0; t<volumes[c].size(); t++) {
            if ((requestedChannel == -1 || c == requestedChannel) && (requestedTimestep == -1 || t == requestedTimestep))
                tgtAssert(volumes.at(c).at(t), "missing volume"); //< current channel/timestep is expected to have been created
            if (volumes[c][t]) {
                result.push_back(volumes.at(c).at(t));
            }
        }
    }
    tgtAssert(!result.empty(), "result volume vector is empty");

    if (getProgressBar())
        getProgressBar()->hide();

    return result;
}


// protected/private methods
// -------------------------

void OMETiffVolumeReader::readTiffDirectory(TIFF* tiffFile, const std::string& dataType, const tgt::svec2& sliceDim, void* destBuffer) const {
    tgtAssert(tiffFile, "null pointer passed");

    // determine bits per voxel from pass data type
    int bitsPerVoxel = 0;
    if (dataType == "uint8" || dataType == "int8")
        bitsPerVoxel = 8;
    else if (dataType == "uint16" || dataType == "int16")
        bitsPerVoxel = 16;
    else if (dataType == "uint32" || dataType == "int32" || dataType == "float")
        bitsPerVoxel = 32;
    else {
        tgtAssert(false, "unknown data type"); //< should have been checked before
        throw tgt::Exception("Unknown data type: " + dataType);
    }
    tgtAssert(bitsPerVoxel > 0, "invalid bits per voxel");

    // read properties from tiff file and check against passed parameters
    uint32 width, height;
    uint16 depth, bps;
    TIFFGetField(tiffFile, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiffFile, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tiffFile, TIFFTAG_SAMPLESPERPIXEL, &depth);
    TIFFGetField(tiffFile, TIFFTAG_BITSPERSAMPLE, &bps);
    if ((sliceDim.x != static_cast<size_t>(width)) || (sliceDim.y != static_cast<size_t>(height))) {
        throw tgt::Exception("Tiff image dimensions (" + genericToString(tgt::svec2(width, height)) + ") differ from volume slice dimensions (" + genericToString(sliceDim) + ")");
    }
    else if (depth != 1) {
        throw tgt::Exception("Samples per pixel != 1 (" + itos(depth) + ")");
    }
    else if (bitsPerVoxel != bps) {
        throw tgt::Exception("Bits per sample (" + itos(bitsPerVoxel) + ") does not match data type " + dataType + " (expected: " + itos(bitsPerVoxel) + ")");
    }

    // determine strip parameters
    tsize_t stripCount = TIFFNumberOfStrips(tiffFile);
    tsize_t stripSize = TIFFStripSize(tiffFile);

    // iterate over strips and copy them to dest buffer
    for (tstrip_t stripID=0; stripID<static_cast<tstrip_t>(stripCount); stripID++) {
        if (TIFFReadEncodedStrip(tiffFile, stripID, destBuffer, stripSize) == -1)
            throw tgt::CorruptedFileException("Failed to read strip " + itos(static_cast<int>(stripID)));
        destBuffer = reinterpret_cast<void*>(reinterpret_cast<char*>(destBuffer) + stripSize);
    }

}

// see OME XML schema definition: http://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2012-06/ome.html
OMETiffStack OMETiffVolumeReader::extractStackInformation(TIFF* tiffFile, const std::string& path) const {
    tgtAssert(tiffFile, "no tiff file");

    OMETiffStack stack;

    char* desc = 0;
    TIFFGetField(tiffFile, TIFFTAG_IMAGEDESCRIPTION, &desc);
    if (!desc)
        throw tgt::Exception("Failed to extract image description from Tiff file");

    TiXmlDocument xmlDoc;
    xmlDoc.Parse(desc);
    if (xmlDoc.Error()) {
        LDEBUG("OME XML string: " << desc);
        throw tgt::Exception("Failed to parse OME XML: " + std::string(xmlDoc.ErrorDesc()));
    }

    // log xml string for debugging
    TiXmlPrinter printer;
    xmlDoc.Accept(&printer);
    //LDEBUG("OME XML string: " << printer.Str());

    // Is there a root element?
    const TiXmlElement* rootElem = xmlDoc.RootElement();
    if (!rootElem)
        throw tgt::Exception("OME XML has no root element");

    // get 'Pixels' element
    const TiXmlNode* pixelsNode = 0;
    try {
        pixelsNode = getXMLNode(rootElem, "Image/Pixels");
    }
    catch (tgt::Exception& e) {
        throw tgt::Exception("Node 'Pixels' not not found in OME XML: " + std::string(e.what()));
    }
    tgtAssert(pixelsNode, "no pixels node");
    const TiXmlElement* pixelsElem = pixelsNode->ToElement();
    if (!pixelsElem)
        throw tgt::Exception("Node 'Pixels' is not an element");

    // extract 'DimensionOrder' from 'Pixels' elem (see )
    if (pixelsElem->QueryValueAttribute("DimensionOrder", &stack.dimensionOrder_) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'DimensionOrder' of element 'Pixels' in OME XML");
    stack.dimensionOrder_ = toUpper(trim(stack.dimensionOrder_));
    if (stack.dimensionOrder_ != "XYZCT" &&
        stack.dimensionOrder_ != "XYZTC" &&
        stack.dimensionOrder_ != "XYCTZ" &&
        stack.dimensionOrder_ != "XYCZT" &&
        stack.dimensionOrder_ != "XYTCZ" &&
        stack.dimensionOrder_ != "XYTZC"    )
    {
        throw VoreenException("Unknown DimensionOrder: " + stack.dimensionOrder_);
    }

    // extract dimensions from 'Pixels' elem
    ivec3 intDim;
    int sizeT, sizeC;
    if (pixelsElem->QueryIntAttribute("SizeX", &intDim.x) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'SizeX' of element 'Pixels' in OME XML");
    if (pixelsElem->QueryIntAttribute("SizeY", &intDim.y) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'SizeY' of element 'Pixels' in OME XML");
    if (pixelsElem->QueryIntAttribute("SizeZ", &intDim.z) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'SizeZ' of element 'Pixels' in OME XML");
    if (pixelsElem->QueryIntAttribute("SizeT", &sizeT) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'SizeT' of element 'Pixels' in OME XML");
    if (pixelsElem->QueryIntAttribute("SizeC", &sizeC) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'SizeC' of element 'Pixels' in OME XML");

    // check and convert dimension information
    if (tgt::hor(tgt::lessThanEqual(intDim, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThan(intDim, tgt::ivec3(99999)))) {
        std::ostringstream stream;
        stream << intDim;
        throw VoreenException("Invalid volume dimensions: " + stream.str());
    }
    stack.volumeDim_ = static_cast<tgt::svec3>(intDim);
    if (sizeC <= 0) {
        throw VoreenException("Invalid sizeC: " + itos(sizeC));
    }
    stack.sizeC_ = static_cast<size_t>(sizeC);
    if (sizeT <= 0) {
        throw VoreenException("Invalid sizeT: " + itos(sizeT));
    }
    stack.sizeT_ = static_cast<size_t>(sizeT);

    // extract data type from 'Pixels' elem
    std::string pixelTypeStr;
    if (pixelsElem->QueryValueAttribute("PixelType", &pixelTypeStr) != TIXML_SUCCESS) {
        LWARNING("Failed to read attribute 'PixelType' of element 'Pixels' in OME XML. Trying attribute 'Type' instead...");
        if (pixelsElem->QueryValueAttribute("Type", &pixelTypeStr) != TIXML_SUCCESS)
            throw tgt::Exception("Failed to read attribute 'PixelType'/'Type' of element 'Pixels' in OME XML");
    }
    stack.datatype_ = pixelTypeStr; //< OME pixel type strings equal type strings used by the VolumeFactory
    if (stack.datatype_ != "uint8" && stack.datatype_ != "int8" && stack.datatype_ != "uint16" && stack.datatype_ != "int16" &&
        stack.datatype_ != "uint32" && stack.datatype_ != "int32" && stack.datatype_ != "float" )
        throw tgt::Exception("Unknown/unsupported PixelType: " + stack.datatype_);

    // extract spacing from 'Pixels' elem
    tgt::vec3 physicalsize;
    if (pixelsElem->QueryFloatAttribute("PhysicalSizeX", &physicalsize.x) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'PhysicalSizeX' of element 'Pixels' in OME XML");
    if (pixelsElem->QueryFloatAttribute("PhysicalSizeY", &physicalsize.y) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'PhysicalSizeY' of element 'Pixels' in OME XML");
    if (pixelsElem->QueryFloatAttribute("PhysicalSizeZ", &physicalsize.z) != TIXML_SUCCESS)
        throw tgt::Exception("Failed to read attribute 'PhysicalSizeZ' of element 'Pixels' in OME XML");

    if (tgt::hor(tgt::lessThan(physicalsize, tgt::vec3::zero))) {
        std::ostringstream stream;
        stream << physicalsize;
        LWARNING("negative physicalsize: " + stream.str() + ", using absolute value");
        physicalsize = tgt::abs(physicalsize);
    }
    else if (tgt::hor(tgt::equal(physicalsize, tgt::vec3::zero))) {
        std::ostringstream stream;
        stream << physicalsize;
        LWARNING("physicalsize is zero : " + stream.str() + ", overwriting with 1.0");
        physicalsize = tgt::vec3(1.f);
    }
    // adapt spacing: OME base length is micron, Voreen base length unit is mm
    stack.voxelSpacing_ = physicalsize / 1000.f;

    // collect 'TiffData' nodes
    std::vector<const TiXmlElement*> tiffdataElems;
    try {
        tiffdataElems = getXMLElementList(pixelsNode, "", "TiffData");
    }
    catch (tgt::Exception& e) {
        throw ("Failed to retrieve 'TiffData' nodes from OME XML: " + std::string(e.what()));
    }

    // extract file information from 'TiffData' nodes
    stack.files_.clear();
    std::string currentUUID;
    try {
        for (size_t i=0; i<tiffdataElems.size(); i++) {
            tgtAssert(tiffdataElems.at(i), "elem is null");
            const TiXmlElement* tiffDataElem = tiffdataElems.at(i);

            // extract UUID element and string (mandatory)
            const TiXmlElement* uuidElem = getXMLElement(tiffDataElem, "UUID");
            tgtAssert(uuidElem, "elem is null"); //< exception expected
            std::string uuid = uuidElem->GetText();
            if (uuid.empty())
                throw tgt::Exception("'UUID' element has no text");
            if (uuid == currentUUID) //< file already processed
                continue;
            else
                currentUUID = uuid;

            // extract filename (mandatory) from UUID elem
            const char* c_filename = uuidElem->Attribute("FileName");
            if (!c_filename)
                throw tgt::Exception("failed to read attribute 'FileName' from 'UUID' node");
            std::string filename(c_filename);
            if (filename.empty())
                throw tgt::Exception("'FileName' attribute of 'UUID' node is empty");
            filename = tgt::FileSystem::cleanupPath(path + "/" + filename);

            // extract FirstZ, FirstC, FirstT attributes from TiffData elem (optional)
            int firstZ = 0;
            int firstC = 0;
            int firstT = 0;
            tiffDataElem->QueryIntAttribute("FirstZ", &firstZ);
            tiffDataElem->QueryIntAttribute("FirstC", &firstC);
            tiffDataElem->QueryIntAttribute("FirstT", &firstT);

            // create file descriptor from extracted information
            stack.files_.push_back(OMETiffFile(filename, 0, firstZ, firstT, firstC));
        }
    }
    catch (tgt::Exception& e) {
        throw tgt::Exception("Failed to collect file information from 'TiffData' nodes in OME XML: " + std::string(e.what()));
    }

    tgtAssert(!stack.files_.empty(), "no files passed");
    tgtAssert(!stack.dimensionOrder_.empty(), "no dimension order passed");
    tgtAssert(!stack.datatype_.empty(), "no dimension order passed");
    tgtAssert(tgt::hand(tgt::greaterThan(stack.voxelSpacing_, tgt::vec3::zero)), "invalid spacing");
    tgtAssert(stack.sizeC_ > 0, "sizeC must be greater 0");
    tgtAssert(stack.sizeT_ > 0, "sizeT must be greater 0");
    //tgtAssert(stack.numSlices_ > 0, "num slices must be greater 0")

    return stack;
}

void OMETiffVolumeReader::determineDirectoryCount(OMETiffStack& stack) const {
    bool estimateDirCount = true;

    VoreenModule* tiffModule = 0;
    BoolProperty* estimateProperty = 0;
    if ((tiffModule = VoreenApplication::app()->getModule("TIFF")) &&
        (estimateProperty = dynamic_cast<BoolProperty*>(tiffModule->getProperty("estimateDirectoryCount"))))
    {
        estimateDirCount = estimateProperty->get();
    }
    else {
        LWARNING("BoolProperty 'estimateDirectoryCount' not found in TIFF module");
    }

    // determine number of directories for each TIFF file
    stack.numSlices_ = 0;

    if (estimateDirCount) { // open only first and last file, and assign their directory counts to all files
        size_t firstFileCount = determineDirectoryCount(stack.files_.front().filename_);
        tgtAssert(firstFileCount > 0, "file has no directories");
        size_t lastFileCount = determineDirectoryCount(stack.files_.back().filename_);
        tgtAssert(lastFileCount > 0, "file has no directories");

        if (firstFileCount == lastFileCount) {
            for (size_t i=0; i<stack.files_.size(); i++) {
                stack.files_.at(i).numDirectories_ = firstFileCount;
                stack.numSlices_ += firstFileCount;
            }
        }
        else {
            LDEBUG("determineDirectoryCount(): directory counts of first and last files differ => disabling estimation");
            estimateDirCount = false;
        }
    }

    if (!estimateDirCount) { // detetermine directory count for all files (potentially time-consuming)
        for (size_t i=0; i<stack.files_.size(); i++) {
            OMETiffFile& file = stack.files_.at(i);
            file.numDirectories_ = 0;

            if (getProgressBar()) {
                getProgressBar()->setProgressMessage("Scanning file '" + tgt::FileSystem::fileName(file.filename_) + "' ...");
                getProgressBar()->setProgress((float)i / stack.files_.size());
                //getProgressBar()->forceUpdate();
            }

            file.numDirectories_ = determineDirectoryCount(file.filename_);
            tgtAssert(file.numDirectories_ > 0, "file has no directories");
            stack.numSlices_ += file.numDirectories_;
        }
    }
    tgtAssert(stack.numSlices_ > 0, "num slices must not be 0");

    // check slice count
    if (stack.numSlices_ != stack.volumeDim_.z*stack.sizeC_*stack.sizeT_) {
        raiseIOException("Total number of slices in TIFF files (" + itos(stack.numSlices_) +
            ") does not match stack size (dim.z*sizeC*sizeT = " + itos(stack.volumeDim_.z*stack.sizeC_*stack.sizeT_) + ") ",
            stack.files_.front().filename_, getProgressBar());
    }

}

size_t OMETiffVolumeReader::determineDirectoryCount(const std::string& filename) const {
    int numDirectories = 0;

    TIFF* tiffFile = TIFFOpen(filename.c_str(), "r");
    if (tiffFile) {
        do {
            numDirectories++;
        } while (TIFFSetDirectory(tiffFile, static_cast<tdir_t>(numDirectories)));
        TIFFClose(tiffFile);
        if (numDirectories == 0)
            throw tgt::Exception("No directories found in TIFF file '" + filename + "'");
    }
    else
        throw tgt::Exception("Failed to open TIFF file '" + filename + "'");

    return numDirectories;
}

const TiXmlNode* OMETiffVolumeReader::getXMLNode(const TiXmlNode* parent, const std::string& path) const {
    tgtAssert(parent, "parent node is null");
    tgtAssert(!path.empty(), "path string is empty");

    std::vector<std::string> pathSplit = strSplit(path, '/');
    if (pathSplit.back().empty())
        pathSplit.pop_back();
    std::string subPath;
    const TiXmlNode* currentNode = parent;
    for (size_t i=0; i<pathSplit.size(); i++) {
        subPath += pathSplit.at(i) + "/";
        const TiXmlNode* childNode = 0;
        while ((childNode = currentNode->IterateChildren(childNode)) != 0) {
            if (childNode->Value() == pathSplit.at(i))
                break;
        }
        if (!childNode)
            throw tgt::Exception("Path '" + subPath  + "' does not exist beneath node '" + parent->Value() + "'");
        currentNode = childNode;
    }

    tgtAssert(currentNode, "node is null");
    return currentNode;
}

const TiXmlElement* OMETiffVolumeReader::getXMLElement(const TiXmlNode* parent, const std::string& path) const {
    // retrieve node
    const TiXmlNode* xmlNode = getXMLNode(parent, path);
    tgtAssert(xmlNode, "null pointer returned"); //< exception expected

    // convert to element
    const TiXmlElement* xmlElem = xmlNode->ToElement();
    if (!xmlElem)
        throw tgt::Exception("Element '" + path + "' beneath parent node '" + parent->Value() + "' does not exist");
    return xmlElem;
}

std::vector<const TiXmlNode*> OMETiffVolumeReader::getXMLNodeList(const TiXmlNode* parent,
        const std::string& path, const std::string& nodeName) const {
    tgtAssert(parent, "parent node is null");
    tgtAssert(!nodeName.empty(), "node name is empty");

    // get parent node
    const TiXmlNode* nodeListParent = 0;
    if (path.empty())
        nodeListParent = parent;
    else
        nodeListParent = getXMLNode(parent, path);
    tgtAssert(nodeListParent, "parent node null"); //< otherwise exception expected

    // iterate over children and collect all nodes with matching name
    std::vector<const TiXmlNode*> result;
    const TiXmlNode* childNode = 0;
    while ((childNode = nodeListParent->IterateChildren(childNode)) != 0) {
        if (childNode->Value() == nodeName)
            result.push_back(childNode);
    }
    return result;
}

std::vector<const TiXmlElement*> OMETiffVolumeReader::getXMLElementList(const TiXmlNode* parent, const std::string& path, const std::string& nodeName) const {
    tgtAssert(parent, "parent node is null");
    tgtAssert(!nodeName.empty(), "node name is empty");

    std::vector<const TiXmlNode*> nodeList = getXMLNodeList(parent, path, nodeName);
    std::vector<const TiXmlElement*> elemList;
    for (size_t i=0; i<nodeList.size(); i++) {
        const TiXmlElement* nodeElem = nodeList.at(i)->ToElement();
        if (!nodeElem)
            throw tgt::Exception("Child node '" + path + "/" + nodeName + "' of parent node '" + parent->Value() + "' is not an element");
        elemList.push_back(nodeElem);
    }
    tgtAssert(elemList.size() == nodeList.size(), "list size mis-match");

    return elemList;
}


/* Previous read slice method using Devil (not used anymore)

void OMETiffVolumeReader::readTiffSliceFileIL(const std::string& filename, const std::string& dataType, const tgt::svec3& volumeDim, void* destBuffer) const {
    tgtAssert(destBuffer, "no dest buffer");

    // determine devil data type
    ILenum ilDataType;
    if (dataType == "uint8")
        ilDataType = IL_UNSIGNED_BYTE;
    else if (dataType == "uint16")
        ilDataType = IL_UNSIGNED_SHORT;
    else if (dataType == "uint32")
        ilDataType = IL_UNSIGNED_INT;
    else if (dataType == "int8")
        ilDataType = IL_BYTE;
    else if (dataType == "int16")
        ilDataType = IL_SHORT;
    else if (dataType == "int32")
        ilDataType = IL_INT;
    else if (dataType == "float")
        ilDataType = IL_FLOAT;
    else {
        tgtAssert(false, "unknown data type"); //< should have been checked before
        throw tgt::Exception("Unknown data type: " + dataType);
    }

    // open image
    ILuint ImageName;
    ilGenImages(1, &ImageName);
    ilBindImage(ImageName);
    if (!ilLoad(IL_TIF, filename.c_str())) {
        int errorCode = ilGetError();
        ilDeleteImages(1, &ImageName);
        throw tgt::CorruptedFileException("Failed to open tiff file (error code=" + itos(errorCode) + ")", filename);
    }

    // check image properties
    if (ilGetInteger(IL_IMAGE_FORMAT) != IL_LUMINANCE) {
        ilDeleteImages(1, &ImageName);
        throw tgt::CorruptedFileException("Expected tiff image format 'IL_LUMINANCE'", filename);
    }
    if (ilGetInteger(IL_IMAGE_TYPE) != ilDataType) {
        ilDeleteImages(1, &ImageName);
        throw tgt::CorruptedFileException("Data type of tiff file does not match data type specified in OME XML", filename);
    }
    if (ilGetInteger(IL_IMAGE_WIDTH) != volumeDim.x) {
        ilDeleteImages(1, &ImageName);
        throw tgt::CorruptedFileException("Image width of tiff file does not match image width specified in OME XML", filename);
    }
    if (ilGetInteger(IL_IMAGE_HEIGHT) != volumeDim.y) {
        ilDeleteImages(1, &ImageName);
        throw tgt::CorruptedFileException("Image height of tiff file does not match image width specified in OME XML", filename);
    }
    if (ilGetInteger(IL_IMAGE_DEPTH) != 1) {
        ilDeleteImages(1, &ImageName);
        throw tgt::CorruptedFileException("Image depth of tiff file is greater 1: " + ilGetInteger(IL_IMAGE_DEPTH), filename);
    }

    // copy image data to corresponding position within target data set
    if (!ilCopyPixels(0, 0, 0, static_cast<ILuint>(volumeDim.x), static_cast<ILuint>(volumeDim.y), 1, IL_LUMINANCE, ilDataType, destBuffer)) {
        ilDeleteImages(1, &ImageName);
        throw tgt::CorruptedFileException("Failed to copy pixel data: " + ilGetError(), filename);
    }

    ilDeleteImages(1, &ImageName);
} */

} // namespace voreen
