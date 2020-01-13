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

#include "surface.h"
#include "voreen/core/voreenapplication.h"
#include "tgt/filesystem.h"


namespace voreen {
/// SurfaceBuilder ---------------------------------------------
SurfaceBuilder::SurfaceBuilder()
    : filename_(VoreenApplication::app()->getUniqueTmpFilePath())
    , numVoxelsStored_(0)
    , file_(filename_, std::ios::binary | std::ios::trunc)
{
}

SurfaceBuilder::SurfaceBuilder(SurfaceBuilder&& other)
    : filename_(other.filename_)
    , numVoxelsStored_(other.numVoxelsStored_)
    , file_(std::move(other.file_))
{
}

StoredSurface SurfaceBuilder::finalize() && {
    StoredSurface ret(filename_, numVoxelsStored_);
    SurfaceBuilder _ = std::move(*this); //destroy builder and thus flush the file
    return ret;
}

void SurfaceBuilder::push(uint64_t linearVoxelPos) {
    ++numVoxelsStored_;
    file_.write(reinterpret_cast<char*>(&linearVoxelPos), sizeof(linearVoxelPos));
}
void SurfaceBuilder::push_all(SurfaceSlice linearVoxelPositions) {
    std::sort(linearVoxelPositions.begin(), linearVoxelPositions.end());
    auto last = unique(linearVoxelPositions.begin(), linearVoxelPositions.end());
    linearVoxelPositions.erase(last, linearVoxelPositions.end());
    numVoxelsStored_ += linearVoxelPositions.size();
    file_.write(reinterpret_cast<char*>(linearVoxelPositions.data()), sizeof(uint64_t) * linearVoxelPositions.size());
}

/// NoFileSurfaceBuilder ---------------------------------------------

#define MAX_UNWRITTEN_SURFACE_VOXELS 0xffff
NoFileSurfaceBuilder::NoFileSurfaceBuilder()
    : filename_(VoreenApplication::app()->getUniqueTmpFilePath())
    , numVoxelsStored_(0)
    , unwrittenVoxels_()
{
    // Precreate and truncate surface file
    std::ofstream file_(filename_, std::ios::binary | std::ios::trunc);
}

NoFileSurfaceBuilder::NoFileSurfaceBuilder(NoFileSurfaceBuilder&& other)
    : filename_(other.filename_)
    , numVoxelsStored_(other.numVoxelsStored_)
    , unwrittenVoxels_(std::move(other.unwrittenVoxels_))
{
}

NoFileSurfaceBuilder::~NoFileSurfaceBuilder() {
    flush();
}

StoredSurface NoFileSurfaceBuilder::finalize() && {
    StoredSurface ret(filename_, numVoxelsStored_);
    NoFileSurfaceBuilder _ = std::move(*this); //destroy builder and thus flush the file
    return ret;
}

void NoFileSurfaceBuilder::flush() {
    if(unwrittenVoxels_.empty()) {
        return;
    }
    std::ofstream file_(filename_, std::ios::binary | std::ios::app);
    file_.write(reinterpret_cast<char*>(unwrittenVoxels_.data()), sizeof(uint64_t) * unwrittenVoxels_.size());

    unwrittenVoxels_.clear();
}

void NoFileSurfaceBuilder::push(uint64_t linearVoxelPos) {
    ++numVoxelsStored_;
    unwrittenVoxels_.push_back(linearVoxelPos);
    if(unwrittenVoxels_.size() >= MAX_UNWRITTEN_SURFACE_VOXELS) {
        flush();
    }
}


/// SurfaceReader ----------------------------------------------
SurfaceReader::SurfaceReader(StoredSurface surface)
    : surface_(surface)
    , file_(surface.filename_, std::ifstream::binary)
{
    tgtAssert(file_.good(), "Opened bad surface file");
}

SurfaceReader::~SurfaceReader() {
    file_.close();
    tgt::FileSystem::deleteFile(surface_.filename_);
}

bool /* success or not */ SurfaceReader::read(uint64_t& val) {
    file_.read(reinterpret_cast<char*>(&val), sizeof(val));
    return !file_.eof();
}

size_t SurfaceReader::numVoxels() const {
    return surface_.numVoxels_;
}

}
