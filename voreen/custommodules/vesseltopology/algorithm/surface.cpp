/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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
#include "surface.h"
#include "voreen/core/voreenapplication.h"
#include "tgt/filesystem.h"


namespace voreen {
/// SurfaceBuilder ---------------------------------------------
SurfaceBuilder::SurfaceBuilder()
    : filename_(VoreenApplication::app()->getUniqueTmpFilePath())
    , file_(filename_, std::ios::binary | std::ios::trunc)
    , numVoxelsStored_(0)
{
}

SurfaceBuilder::SurfaceBuilder(SurfaceBuilder&& other)
    : filename_(other.filename_)
    , file_(std::move(other.file_))
    , numVoxelsStored_(other.numVoxelsStored_)
{
}

StoredSurface SurfaceBuilder::finalize(SurfaceBuilder&& builder) {
    StoredSurface ret(builder.filename_, builder.numVoxelsStored_);
    SurfaceBuilder _ = std::move(builder); //destroy builder and thus flush the file
    return ret;
}

void SurfaceBuilder::push(uint64_t linearVoxelPos) {
    ++numVoxelsStored_;
    file_.write(reinterpret_cast<char*>(&linearVoxelPos), sizeof(linearVoxelPos));
}
void SurfaceBuilder::push_all(std::set<uint64_t> linearVoxelPositions) {
    for(uint64_t linearpos : linearVoxelPositions) {
        push(linearpos);
    }
}
/// SurfaceReader ----------------------------------------------
SurfaceReader::SurfaceReader(StoredSurface surface)
    : surface_(surface)
    , file_(surface.filename_, std::ifstream::binary)
{
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
