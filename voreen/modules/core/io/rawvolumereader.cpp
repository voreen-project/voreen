/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "rawvolumereader.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

#include "tgt/exception.h"
#include "tgt/filesystem.h"

#include "voreen/core/io/progressbar.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumefusion.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorswapendianness.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

using tgt::ivec3;
using tgt::vec3;

namespace voreen {

const std::string RawVolumeReader::loggerCat_ = "voreen.RawVolumeReader";

RawVolumeReader::ReadHints::ReadHints(tgt::ivec3 dimensions, tgt::vec3 spacing,
                                      const std::string& objectModel, const std::string& format,
                                      int timeframe, size_t headerskip, bool bigEndian)
    : dimensions_(dimensions), spacing_(spacing),
      objectModel_(objectModel), format_(format), timeframe_(timeframe), headerskip_(headerskip),
      bigEndianByteOrder_(bigEndian),
      transformation_(tgt::mat4::identity),
      modality_(Modality::MODALITY_UNKNOWN),
      hash_(""),
      timeStep_(-1.f),
      sliceOrder_("+z")
{}

RawVolumeReader::RawVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("raw");
    protocols_.push_back("raw");
}

void RawVolumeReader::setReadHints(ivec3 dimensions, vec3 spacing,
                                const std::string& objectModel, const std::string& format,
                                int timeframe, int headerskip, bool bigEndian)
{
    hints_ = ReadHints(dimensions, spacing, objectModel, format, timeframe, headerskip, bigEndian);
}

void RawVolumeReader::setReadHints(const ReadHints& hints) {
    hints_ = hints;
}

VolumeList* RawVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    if (origin.getProtocol() == "raw" && !origin.getSearchString().empty())
        hints_ = extractReadHintsFromOrigin(origin);

    return readSlices(url, 0, 0);
}

Volume* RawVolumeReader::read(const VolumeURL& origin) {
    // read parameters from origin
    std::string filename = origin.getPath();
    hints_ = extractReadHintsFromOrigin(origin);

    // load volume
    Volume* handle = 0;
    VolumeList* collection = read(filename);
    if (!collection->empty())
        handle = static_cast<Volume*>(collection->first());

    delete collection;

    oldVolumePosition(handle);
    return handle;
}

VolumeList* RawVolumeReader::readSlices(const std::string &url, size_t firstSlice, size_t lastSlice ) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    ReadHints& h = hints_;

    // check dimensions
    if (tgt::hor(tgt::lessThan(h.dimensions_, ivec3(0))) || tgt::hor(tgt::greaterThan(h.dimensions_, ivec3(10000)))) {
        LERROR("Invalid volume dimensions: " << h.dimensions_);
        return 0;
    }

    // check if we have to read only some slices instead of the whole volume.
    if ( ! (firstSlice==0 && lastSlice==0)) {
        if (lastSlice > firstSlice) {
            h.dimensions_.z = static_cast<int>(lastSlice - firstSlice);
        }
    }

    std::string info = "Loading raw file " + fileName + " ";

    if (h.dimensions_ == tgt::ivec3::zero)
        throw tgt::CorruptedFileException("No readHints set.", fileName);

    FILE* fin;
    fin = fopen(fileName.c_str(),"rb");

    if (fin == 0)
        throw tgt::IOException("Unable to open raw file for reading", fileName);

    VolumeRAM* volume = NULL;

    if (h.objectModel_ == "I") {
        if (h.format_ == "UCHAR") {
            LINFO(info << "(8 bit dataset)");
            VolumeRAM_UInt8* v = new VolumeRAM_UInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "CHAR") {
            LINFO(info << "(8 bit signed dataset)");
            VolumeRAM_Int8* v = new VolumeRAM_Int8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "USHORT"|| h.format_ == "USHORT_12") {
            LINFO(info << "(16 bit dataset)");
            VolumeRAM_UInt16* v = new VolumeRAM_UInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "SHORT") {
            LINFO(info << "(16 bit signed dataset)");
            VolumeRAM_Int16* v = new VolumeRAM_Int16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT") {
            LINFO(info << "(32 bit dataset)");
            VolumeRAM_UInt32* v = new VolumeRAM_UInt32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT") {
            LINFO(info << "(32 bit signed dataset)");
            VolumeRAM_Int32* v = new VolumeRAM_Int32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT64") {
            LINFO(info << "(64 bit dataset)");
            VolumeRAM_UInt64* v = new VolumeRAM_UInt64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT64") {
            LINFO(info << "(64 bit signed dataset)");
            VolumeRAM_Int64* v = new VolumeRAM_Int64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "FLOAT") {
            LINFO(info << "(32 bit float dataset)");
            VolumeRAM_Float* v = new VolumeRAM_Float(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "DOUBLE") {
            LINFO(info << "(64 bit double dataset)");
            VolumeRAM_Double* v = new VolumeRAM_Double(h.dimensions_);
            volume = v;
        }
        else {
            fclose(fin);
            throw tgt::CorruptedFileException("Format '" + h.format_ + "' not supported", fileName);
        }
    }
    else if (h.objectModel_ == "LA") { // luminance alpha
        if (h.format_ == "UCHAR") {
            LINFO(info << "(2x8 bit dataset)");
            VolumeRAM_2xUInt8* v = new VolumeRAM_2xUInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "CHAR") {
            LINFO(info << "(2x8 bit signed dataset)");
            VolumeRAM_2xInt8* v = new VolumeRAM_2xInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "USHORT" || h.format_ == "USHORT_12") {
            LINFO(info << "(2x16 bit dataset)");
            VolumeRAM_2xUInt16* v = new VolumeRAM_2xUInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "SHORT") {
            LINFO(info << "(2x16 bit signed dataset)");
            VolumeRAM_2xInt16* v = new VolumeRAM_2xInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT") {
            LINFO(info << "(2x32 bit dataset)");
            VolumeRAM_2xUInt32* v = new VolumeRAM_2xUInt32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT") {
            LINFO(info << "(2x32 bit signed dataset)");
            VolumeRAM_2xInt32* v = new VolumeRAM_2xInt32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT64") {
            LINFO(info << "(2x64 bit dataset)");
            VolumeRAM_2xUInt64* v = new VolumeRAM_2xUInt64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT64") {
            LINFO(info << "(2x64 bit signed dataset)");
            VolumeRAM_2xInt64* v = new VolumeRAM_2xInt64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "FLOAT") {
            LINFO(info << "(2x32 bit float dataset)");
            VolumeRAM_2xFloat* v = new VolumeRAM_2xFloat(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "DOUBLE") {
            LINFO(info << "(2x64 bit double dataset)");
            VolumeRAM_2xDouble* v = new VolumeRAM_2xDouble(h.dimensions_);
            volume = v;
        }
        else {
            fclose(fin);
            throw tgt::CorruptedFileException("Format '" + h.format_ + "' not supported", fileName);
        }
    }
    else if (h.objectModel_ == "RGB") {
        if (h.format_ == "UCHAR") {
            LINFO(info << "(3x8 bit dataset)");
            VolumeRAM_3xUInt8* v = new VolumeRAM_3xUInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "CHAR") {
            LINFO(info << "(3x8 bit signed dataset)");
            VolumeRAM_3xInt8* v = new VolumeRAM_3xInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "USHORT" || h.format_ == "USHORT_12") {
            LINFO(info << "(3x16 bit dataset)");
            VolumeRAM_3xUInt16* v = new VolumeRAM_3xUInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "SHORT") {
            LINFO(info << "(3x16 bit signed dataset)");
            VolumeRAM_3xInt16* v = new VolumeRAM_3xInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT") {
            LINFO(info << "(3x32 bit dataset)");
            VolumeRAM_3xUInt32* v = new VolumeRAM_3xUInt32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT") {
            LINFO(info << "(3x32 bit signed dataset)");
            VolumeRAM_3xInt32* v = new VolumeRAM_3xInt32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT64") {
            LINFO(info << "(3x64 bit dataset)");
            VolumeRAM_3xUInt64* v = new VolumeRAM_3xUInt64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT64") {
            LINFO(info << "(3x64 bit signed dataset)");
            VolumeRAM_3xInt64* v = new VolumeRAM_3xInt64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "FLOAT") {
            LINFO(info << "(3x32 bit float dataset)");
            VolumeRAM_3xFloat* v = new VolumeRAM_3xFloat(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "DOUBLE") {
            LINFO(info << "(3x64 bit double dataset)");
            VolumeRAM_3xDouble* v = new VolumeRAM_3xDouble(h.dimensions_);
            volume = v;
        }
        else {
            fclose(fin);
            throw tgt::CorruptedFileException("Format '" + h.format_ + "' not supported", fileName);
        }
    }
    else if (h.objectModel_ == "RGBA") {
        if (h.format_ == "UCHAR") {
            LINFO(info << "(4x8 bit dataset)");
            VolumeRAM_4xUInt8* v = new VolumeRAM_4xUInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "CHAR") {
            LINFO(info << "(3x8 bit signed dataset)");
            VolumeRAM_4xInt8* v = new VolumeRAM_4xInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "USHORT" || h.format_ == "USHORT_12") {
            LINFO(info << "(4x16 bit dataset)");
            VolumeRAM_4xUInt16* v = new VolumeRAM_4xUInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "SHORT") {
            LINFO(info << "(4x16 bit signed dataset)");
            VolumeRAM_4xInt16* v = new VolumeRAM_4xInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT") {
            LINFO(info << "(4x32 bit dataset)");
            VolumeRAM_4xUInt32* v = new VolumeRAM_4xUInt32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT") {
            LINFO(info << "(4x32 bit signed dataset)");
            VolumeRAM_4xInt32* v = new VolumeRAM_4xInt32(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "UINT64") {
            LINFO(info << "(4x64 bit dataset)");
            VolumeRAM_4xUInt64* v = new VolumeRAM_4xUInt64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "INT64") {
            LINFO(info << "(4x64 bit signed dataset)");
            VolumeRAM_4xInt64* v = new VolumeRAM_4xInt64(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "FLOAT") {
            LINFO(info << "(4x32 bit float dataset)");
            VolumeRAM_4xFloat* v = new VolumeRAM_4xFloat(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "DOUBLE") {
            LINFO(info << "(4x64 bit double dataset)");
            VolumeRAM_4xDouble* v = new VolumeRAM_4xDouble(h.dimensions_);
            volume = v;
        }
        else {
            fclose(fin);
            throw tgt::CorruptedFileException("Format '" + h.format_ + "' not supported", fileName);
        }
    }
    else if (h.objectModel_ == "MAT3") { // luminance alpha
        if (h.format_ == "FLOAT") {
            LINFO(info << "(9x32 bit float (rank 3 matrix) dataset");
            VolumeRAM_Mat3Float* v = new VolumeRAM_Mat3Float(h.dimensions_);
            volume = v;
        }
    }
    else if (h.objectModel_.find("TENSOR_") == 0) {
        if (h.format_ == "FLOAT") {
            LINFO(info << "(6x32 bit float (second order tensor) dataset)");
            VolumeRAM_Tensor2Float* v = new VolumeRAM_Tensor2Float(h.dimensions_);
            volume = v;
        }
    }
    else {
        fclose(fin);
        throw tgt::CorruptedFileException("unsupported ObjectModel '" + h.objectModel_ + "'", fileName);
    }

    // Calculate additional skipping if we have to read only slices or not the first time frame
    uint64_t dimx = static_cast<uint64_t>(h.dimensions_.x);
    uint64_t dimy = static_cast<uint64_t>(h.dimensions_.y);
    uint64_t dimz = static_cast<uint64_t>(h.dimensions_.z);
    uint64_t numBytes = static_cast<uint64_t>(volume->getBitsAllocated() / 8);
    uint64_t sliceSkip = dimx * dimy * static_cast<uint64_t>(firstSlice) * numBytes;
    uint64_t frameSkip = dimx * dimy * dimz * static_cast<uint64_t>(h.timeframe_) * numBytes;

    // now add that to the headerskip we might have received
    uint64_t offset = h.headerskip_ + sliceSkip + frameSkip;

    #ifdef _MSC_VER
        _fseeki64(fin, offset, SEEK_SET);
    #else
        fseek(fin, offset, SEEK_SET);
    #endif

    volume->clear();

    if (getProgressBar()) {
        getProgressBar()->setTitle("Loading Volume");
        // getProgress()->setMessage("Loading volume: " + tgt::FileSystem::fileName(fileName));
        getProgressBar()->setProgressMessage("Loading volume: " + fileName);
    }
    VolumeReader::read(volume, fin);

    if (lastSlice == 0) {
        if (feof(fin) ) {
            fclose(fin);
            delete volume;
            if (getProgressBar())
                getProgressBar()->hide();
            // throw exception
            throw tgt::CorruptedFileException("unexpected EOF: raw file truncated or ObjectModel '" +
                                              h.objectModel_ + "' invalid", fileName);
        }
    }

    fclose(fin);

    // correct tensor layout
    if (h.objectModel_.find("TENSOR_") == 0 && h.format_ == "FLOAT") {
        VolumeRAM* shifted = new VolumeRAM_Tensor2Float(h.dimensions_);
        VolumeRAM_Tensor2Float* foo = static_cast<VolumeRAM_Tensor2Float*>(shifted);
        float* plainData = reinterpret_cast<float*>(volume->getData());

        // convert tensors from VolumeFusion like layout to VolumeRAM_Tensor2Float layout
        if (h.objectModel_.find("TENSOR_FUSION_") == 0) {
            for (size_t elem = 0; elem < 6; ++elem) {
                size_t count = h.dimensions_.x * h.dimensions_.y * h.dimensions_.z;
                if (h.objectModel_ == "TENSOR_FUSION_UP") {
                    for (size_t i = 0; i < count; ++i)
                        foo->voxel(i) = Tensor2<float>(plainData[i], plainData[count+i], plainData[2*count + i], plainData[3*count + i], plainData[4*count + i], plainData[5*count + i]);
                }
                else if (h.objectModel_ == "TENSOR_FUSION_LOW") {
                    for (size_t i = 0; i < count; ++i)
                        foo->voxel(i) = Tensor2<float>::createTensorFromLowerDiagonalMatrix(plainData[i], plainData[count+i], plainData[2*count + i], plainData[3*count + i], plainData[4*count + i], plainData[5*count + i]);
                }
                else if (h.objectModel_ == "TENSOR_FUSION_DIAG") {
                    for (size_t i = 0; i < count; ++i)
                        foo->voxel(i) = Tensor2<float>::createTensorFromDiagonalOrder(plainData[i], plainData[count+i], plainData[2*count + i], plainData[3*count + i], plainData[4*count + i], plainData[5*count + i]);
                }
            }
        }
        // convert tensor element ordering
        else {
            if (h.objectModel_ == "TENSOR_UP") {
                for (size_t i = 0; i < volume->getNumVoxels(); ++i)
                    foo->voxel(i) = Tensor2<float>(plainData + 6*i);
            }
            else if (h.objectModel_ == "TENSOR_LOW") {
                for (size_t i = 0; i < volume->getNumVoxels(); ++i)
                    foo->voxel(i) = Tensor2<float>::createTensorFromLowerDiagonalMatrix(plainData + 6*i);
            }
            else if (h.objectModel_ == "TENSOR_DIAG") {
                for (size_t i = 0; i < volume->getNumVoxels(); ++i)
                    foo->voxel(i) = Tensor2<float>::createTensorFromDiagonalOrder(plainData + 6*i);
            }
        }

        std::swap(volume, shifted);
        delete shifted;
    }

    if (h.sliceOrder_ == "-x") {
        LINFO("slice order is -x, reversing order to +x...\n");
        reverseXSliceOrder(volume);
    } else if (h.sliceOrder_ == "-y") {
        LINFO("slice order is -y, reversing order to +y...\n");
        reverseYSliceOrder(volume);
    } else if (h.sliceOrder_ == "-z") {
        LINFO("slice order is -z, reversing order to +z...\n");
        reverseZSliceOrder(volume);
    }

    vec3 offs(0.0f);
    Volume* volumeHandle = new Volume(volume, h.spacing_, offs);
    volumeHandle->setPhysicalToWorldMatrix(h.transformation_);
    volumeHandle->setModality(h.modality_);
    volumeHandle->setTimestep(static_cast<float>(h.timeframe_));

    if (h.bigEndianByteOrder_) {
        VolumeOperatorSwapEndianness::APPLY_OP(volumeHandle);
    }

    if(!h.hash_.empty())
        volumeHandle->setHash(h.hash_);

    volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));
    volumeHandle->setOrigin(VolumeURL("raw", fileName, encodeReadHintsIntoSearchString(hints_)));

    VolumeList* volumeList = new VolumeList();
    volumeList->add(volumeHandle);

    if (getProgressBar())
        getProgressBar()->hide();

    return volumeList;
}

VolumeList* RawVolumeReader::readBrick(const std::string &url, tgt::ivec3 brickStartPos,int brickSize) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    ReadHints& h = hints_;
    tgt::ivec3 datasetDims = h.dimensions_;
    h.dimensions_ = tgt::ivec3(brickSize);

    if (h.dimensions_ == tgt::ivec3::zero) {
        throw tgt::CorruptedFileException("No readHints set.", fileName);
    }

    FILE* fin;
    fin = fopen(fileName.c_str(),"rb");

    if (fin == 0)
        throw tgt::IOException("Unable to open raw file for reading", fileName);

    VolumeRAM* volume;

    if (h.objectModel_ == "I") {
        if (h.format_ == "UCHAR") {
            VolumeRAM_UInt8* v = new VolumeRAM_UInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "USHORT" || h.format_ == "USHORT_12") {
            VolumeRAM_UInt16* v = new VolumeRAM_UInt16(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "FLOAT") {
            VolumeRAM_Float* v = new VolumeRAM_Float(h.dimensions_);
            volume = v;
        }
        else {
            fclose(fin);
            throw tgt::CorruptedFileException("Format '" + h.format_ + "' not supported", fileName);
        }
    }
    else if (h.objectModel_ == "RGBA") {
        if (h.format_ == "UCHAR") {
            //LINFO("Reading 4x8 bit dataset");
            VolumeRAM_4xUInt8* v = new VolumeRAM_4xUInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "USHORT") {
            //LINFO("Reading 4x16 bit dataset");
            VolumeRAM_4xUInt16* v = new VolumeRAM_4xUInt16(h.dimensions_);
            volume = v;
        }
        else {
            fclose(fin);
            throw tgt::CorruptedFileException("Format '" + h.format_ + "' not supported for object model RGBA", fileName);
        }
    }
    else if (h.objectModel_ == "RGB") {
        if (h.format_ == "UCHAR") {
            //LINFO("Reading 3x8 bit dataset");
            VolumeRAM_3xUInt8* v = new VolumeRAM_3xUInt8(h.dimensions_);
            volume = v;
        }
        else if (h.format_ == "USHORT") {
            //LINFO("Reading 3x16 bit dataset");
            VolumeRAM_3xUInt16* v = new VolumeRAM_3xUInt16(h.dimensions_);
            volume = v;
        } else if (h.format_ == "FLOAT") {
            VolumeRAM_3xFloat* v = new VolumeRAM_3xFloat(h.dimensions_);
            volume = v;
        } else {
            fclose(fin);
            throw tgt::CorruptedFileException("Format '" + h.format_ + "' not supported for object model RGB", fileName);
        }
    }
    else if (h.objectModel_ == "LA") { // luminance alpha
        //LINFO("Reading luminance16 alpha16 dataset");
        VolumeRAM_4xUInt8* v = new VolumeRAM_4xUInt8(h.dimensions_);
        volume = v;
    }
    else {
        fclose(fin);
        throw tgt::CorruptedFileException("unsupported ObjectModel '" + h.objectModel_ + "'", fileName);
    }

    volume->clear();

    int voxelSize = static_cast<int>(volume->getBitsAllocated() / 8);

    uint64_t initialSeekPos = (brickStartPos.z *(datasetDims.x*datasetDims.y)*brickSize*voxelSize )+
        (brickStartPos.y * datasetDims.x*brickSize*voxelSize) + (brickStartPos.x*brickSize*voxelSize);

    // add header skip
    initialSeekPos += h.headerskip_;

    #ifdef _MSC_VER
        _fseeki64(fin, initialSeekPos, SEEK_SET);
    #else
        fseek(fin, initialSeekPos, SEEK_SET);
    #endif

        for (int i=0; i < brickSize; i++) {
            for (int j=0; j < brickSize; j++) {
                int volumePos = brickSize*j + (i*brickSize*brickSize);
                volumePos = volumePos*voxelSize;
                if (fread(reinterpret_cast<char*>(volume->getData()) + volumePos, 1, brickSize*voxelSize, fin) == 0)
                    LWARNING("fread() failed");
                #ifdef _MSC_VER
                    _fseeki64(fin, (datasetDims.x - brickSize)*voxelSize, SEEK_CUR);
                #else
                    fseek(fin, (datasetDims.x - brickSize)*voxelSize, SEEK_CUR);
                #endif
            }
            #ifdef _MSC_VER
                _fseeki64(fin, (datasetDims.y - brickSize)*datasetDims.x*voxelSize, SEEK_CUR);
            #else
                fseek(fin, (datasetDims.y - brickSize)*datasetDims.x*voxelSize, SEEK_CUR);
            #endif
        }

    fclose(fin);

    VolumeList* volumeList = new VolumeList();
    Volume* volumeHandle = new Volume(volume, h.spacing_, vec3(0.0f), h.transformation_);
    volumeHandle->setModality(h.modality_);

    // encode raw parameters into search string
    std::ostringstream searchStream;
    searchStream << "objectModel=" << h.objectModel_ << "&";
    searchStream << "format=" << h.format_ << "&";
    searchStream << "headerskip=" << h.headerskip_ << "&";
    if (h.bigEndianByteOrder_)
        searchStream << "bigEndian=" << h.bigEndianByteOrder_ << "&";
    searchStream << "dim_x=" << h.dimensions_.x << "&";
    searchStream << "dim_y=" << h.dimensions_.y << "&";
    searchStream << "dim_z=" << h.dimensions_.z << "&";
    searchStream << "spacing_x=" << h.spacing_.x << "&";
    searchStream << "spacing_y=" << h.spacing_.y << "&";
    searchStream << "spacing_z=" << h.spacing_.z << "&";
    volumeHandle->setOrigin(VolumeURL("raw", fileName, searchStream.str()));

    volumeList->add(volumeHandle);

    return volumeList;
}

Volume* RawVolumeReader::readSliceStack(const std::vector<std::string>& sliceFiles) {

    if (hints_.dimensions_.z > 1) {
        LWARNING("Setting z-dimension of slices to 1.");
        hints_.dimensions_.z = 1;
    }

    std::vector<VolumeRAM*> volumes;
    for (size_t i=0; i<sliceFiles.size(); i++) {
        try {
            VolumeList* collection = read(sliceFiles[i]);
            if (!collection->empty() && collection->first()->getRepresentation<VolumeRAM>()) {
                volumes.push_back(static_cast<Volume*>(collection->first())->getWritableRepresentation<VolumeRAM>());
            }
            delete collection;
        }
        catch (tgt::FileException& e) {
            LWARNING("Reading slice '" << sliceFiles[i] << "' failed: " << e.what());
        }
        catch (std::bad_alloc&) {
            LWARNING("Reading slice '" << sliceFiles[i] << "' failed: bad allocation");
        }
    }

    if (volumes.empty()) {
        LWARNING("No slices");
        return 0;
    }

    LINFO("Constructing volume from " << volumes.size() << " slices of dimensions " << hints_.dimensions_.xy());

    VolumeRAM* result = 0;
    try {
        VolumeRAM* v = volumes.front();
        Volume vh(v, hints_.spacing_, vec3(0.0f));
        Volume* t = VolumeOperatorResize::APPLY_OP(&vh,
            ivec3(hints_.dimensions_.x, hints_.dimensions_.y, static_cast<int>(volumes.size())));
        result = t->getWritableRepresentation<VolumeRAM>();
        t->releaseAllRepresentations();
        delete t;
        vh.releaseAllRepresentations();
        result->clear();
        //result->setTransformation(tgt::mat4::identity);
    }
    catch (std::bad_alloc&) {
        LERROR("Reading slice stack failed: bad allocation");
    }

    // copy slices to output volume
    if (result) {
        size_t sliceDataSize = volumes.front()->getNumBytes();
        tgtAssert(result->getNumBytes() == sliceDataSize*volumes.size(), "Invalid slice sizes");

        for (size_t i=0; i<volumes.size(); i++) {

            if (volumes[i]->getNumBytes() != sliceDataSize) {
                LERROR("Slice '" << sliceFiles[i] << "' differs in size.");
                continue;
            }

            memcpy(((char*)result->getData())+i*sliceDataSize, (char*)volumes[i]->getData(), sliceDataSize);
        }
    }


    // free slices
    for (size_t i=0; i<volumes.size(); i++) {
        delete volumes[i];
    }

    if (result) {
        Volume* outputHandle = new Volume(result, hints_.spacing_, vec3(0.0f));
        outputHandle->setTimestep(hints_.timeStep_);

        // encode raw parameters into search string (currently only first slice)
        std::ostringstream searchStream;
        searchStream << "objectModel=" << hints_.objectModel_ << "&";
        searchStream << "format=" << hints_.format_ << "&";
        searchStream << "headerskip=" << hints_.headerskip_ << "&";
        if (hints_.bigEndianByteOrder_)
            searchStream << "bigEndian=" << hints_.bigEndianByteOrder_ << "&";
        searchStream << "dim_x=" << hints_.dimensions_.x << "&";
        searchStream << "dim_y=" << hints_.dimensions_.y << "&";
        searchStream << "dim_z=" << hints_.dimensions_.z << "&";
        searchStream << "spacing_x=" << hints_.spacing_.x << "&";
        searchStream << "spacing_y=" << hints_.spacing_.y << "&";
        searchStream << "spacing_z=" << hints_.spacing_.z << "&";

        outputHandle->setOrigin(VolumeURL("raw", sliceFiles.front(), searchStream.str()));

        return outputHandle;
    }
    else {
        return 0;
    }
}

VolumeReader* RawVolumeReader::create(ProgressBar* progress) const {
    return new RawVolumeReader(progress);
}

RawVolumeReader::ReadHints RawVolumeReader::extractReadHintsFromOrigin(const VolumeURL& origin) const {
    ReadHints hints;

    hints.objectModel_ = origin.getSearchParameter("objectModel");
    hints.format_ = origin.getSearchParameter("format");

    std::istringstream s0(origin.getSearchParameter("timeframe"));
    s0 >> hints.timeframe_;

    std::istringstream s1(origin.getSearchParameter("dim_x"));
    s1 >> hints.dimensions_.x;
    s1.str(origin.getSearchParameter("dim_y"));
    s1.clear();
    s1 >> hints.dimensions_.y;
    s1.str(origin.getSearchParameter("dim_z"));
    s1.clear();
    s1 >> hints.dimensions_.z;

    std::istringstream s2(origin.getSearchParameter("spacing_x"));
    s2 >> hints.spacing_.x;
    s2.str(origin.getSearchParameter("spacing_y"));
    s2.clear();
    s2 >> hints.spacing_.y;
    s2.str(origin.getSearchParameter("spacing_z"));
    s2.clear();
    s2 >> hints.spacing_.z;

    std::istringstream s3(origin.getSearchParameter("headerskip"));
    s3 >> hints.headerskip_;

    int byteOrder;
    std::istringstream s4(origin.getSearchParameter("bigEndian"));
    s4 >> byteOrder;
    hints.bigEndianByteOrder_ = (byteOrder == 1);

    return hints;
}

std::string RawVolumeReader::encodeReadHintsIntoSearchString(const RawVolumeReader::ReadHints& h) const {
    // encode raw parameters into search string
    std::ostringstream searchStream;
    searchStream << "objectModel=" << h.objectModel_ << "&";
    searchStream << "format=" << h.format_ << "&";
    searchStream << "headerskip=" << h.headerskip_ << "&";
    searchStream << "timeframe=" << h.timeframe_ << "&";
    if (h.bigEndianByteOrder_)
        searchStream << "bigEndian=" << h.bigEndianByteOrder_ << "&";
    searchStream << "dim_x=" << h.dimensions_.x << "&";
    searchStream << "dim_y=" << h.dimensions_.y << "&";
    searchStream << "dim_z=" << h.dimensions_.z << "&";
    searchStream << "spacing_x=" << h.spacing_.x << "&";
    searchStream << "spacing_y=" << h.spacing_.y << "&";
    searchStream << "spacing_z=" << h.spacing_.z << "&";

    return searchStream.str();
}

}   // namespace voreen
