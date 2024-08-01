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

#ifndef VRN_ALGDISTANCETRANSFORM_H
#define VRN_ALGDISTANCETRANSFORM_H

#include "voreen/core/datastructures/volume/volumebase.h"
#include "../datastructures/lz4slicevolume.h"

#include <string>

namespace voreen {

LZ4SliceVolume<float> compute_distance_transform(const VolumeBase& vol, float binarizationThreshold, std::string outputPath, ProgressReporter& progressReporter);

enum MedialStructureType {
    MedialStructureType_Point = 3,
    MedialStructureType_Line = 2,
    MedialStructureType_Surface = 1,
};

// Warning: This is only somewhat working and creates a lot of noisy and disconnected structures.
LZ4SliceVolume<uint8_t> compute_medial_structures(const VolumeBase& vol, float binarizationThreshold, MedialStructureType structureType, std::string outputPath, ProgressReporter& progressReporter);

}   //namespace

#endif
