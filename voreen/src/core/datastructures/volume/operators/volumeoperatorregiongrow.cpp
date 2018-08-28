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

#include "voreen/core/datastructures/volume/operators/volumeoperatorregiongrow.h"

namespace voreen {

const std::string VolumeOperatorRegionGrowBase::loggerCat_("voreen.VolumeOperatorRegionGrow");

float neighborStandardDeviation(const tgt::ivec3& pos, const VolumeRAM* dataset, Statistics& stats) {
    for (int z=-1; z <= 1; z++)
        for (int y=-1; y <= 1; y++)
            for (int x=-1; x <= 1; x++) {
                tgt::ivec3 p(x, y, z);
                if (p != tgt::ivec3::zero)
                    stats.addSample(dataset->getVoxelNormalized(pos + p));
            }

    return stats.getStdDev();
}

float gradientMagnitude(const tgt::ivec3& pos, const VolumeRAM* dataset, const tgt::vec3 spacing) {
    float v0, v1, v2, v3, v4, v5;
    tgt::vec3 gradient;
    float vCenter = dataset->getVoxelNormalized(pos);

    if (static_cast<size_t>(pos.x) != dataset->getDimensions().x-1)
        v0 = dataset->getVoxelNormalized(pos + tgt::ivec3(1, 0, 0));
    else
        v0 = vCenter; // Use center if voxel is on border

    if (static_cast<size_t>(pos.y) != dataset->getDimensions().y-1)
        v1 = dataset->getVoxelNormalized(pos + tgt::ivec3(0, 1, 0));
    else
        v1 = vCenter;

    if (static_cast<size_t>(pos.z) != dataset->getDimensions().z-1)
        v2 = dataset->getVoxelNormalized(pos + tgt::ivec3(0, 0, 1));
    else
        v2 = vCenter;

    if (pos.x != 0)
        v3 = dataset->getVoxelNormalized(pos + tgt::ivec3(-1, 0, 0));
    else
        v3 = vCenter;

    if (pos.y != 0)
        v4 = dataset->getVoxelNormalized(pos + tgt::ivec3(0, -1, 0));
    else
        v4 = vCenter;

    if (pos.z != 0)
        v5 = dataset->getVoxelNormalized(pos + tgt::ivec3(0, 0, -1));
    else
        v5 = vCenter;

    gradient = tgt::vec3(static_cast<float>(v3 - v0), static_cast<float>(v4 - v1), static_cast<float>(v5 - v2));
    gradient /= (spacing * 2.0f);

    return length(gradient);
}

float neighborStandardDeviationGradients(const tgt::ivec3& pos, const VolumeRAM* dataset, const tgt::vec3 spacing, Statistics& stats) {
    for (int z=-1; z <= 1; z++)
        for (int y=-1; y <= 1; y++)
            for (int x=-1; x <= 1; x++) {
                tgt::ivec3 p(x, y, z);
                if (p != tgt::ivec3::zero)
                    stats.addSample(gradientMagnitude(pos + p, dataset, spacing));
            }

    return stats.getStdDev();
}

} // namespace voreen

