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

#include "volumedisklz4.h"

#include "voreen/core/utils/hashing.h"


namespace voreen {

const std::string VolumeDiskLZ4::loggerCat_("voreen.hdf5.VolumeDiskLZ4");


VolumeDiskLZ4::VolumeDiskLZ4(std::unique_ptr<LZ4SliceVolumeBase> volume)
    : VolumeDisk(volume->getMetaData().getBaseType(), volume->getDimensions())
    , volume_(std::move(volume))
{
}

VolumeDiskLZ4::~VolumeDiskLZ4() {
}

std::string VolumeDiskLZ4::getHash() const {
    std::string configStr;

    configStr += volume_->getFilePath() + "#";

    return VoreenHash::getHash(configStr);
}

VolumeRAM* VolumeDiskLZ4::loadVolume() const {
    return volume_->loadBaseSlab(0, volume_->getNumSlices()).release();
}

VolumeRAM* VolumeDiskLZ4::loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
    return volume_->loadBaseSlab(firstZSlice, lastZSlice+1).release();
}

VolumeRAM* VolumeDiskLZ4::loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
    throw tgt::Exception("loadBrick is not supported");
}

} // namespace voreen