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

#ifndef VRN_VOLUMEOPERATORREGIONGROW_H
#define VRN_VOLUMEOPERATORREGIONGROW_H

#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/utils/statistics.h"

namespace voreen {

// Calculate standard deviation of the intensity values from the 26 neighbors of the voxel at
// position pos. Note: position must not be a border voxel.
float neighborStandardDeviation(const tgt::ivec3& pos, const VolumeRAM* dataset, Statistics& stats);

// Calculate gradient magnitude. Note: position must not be a border voxel.
float gradientMagnitude(const tgt::ivec3& pos, const VolumeRAM* dataset, const tgt::vec3 spacing);

// Calculate standard deviation from the gradient magnitudes of the 26 neighbors of the voxel at
// position pos. Note: position must be at least 2 voxels away from the border.
float neighborStandardDeviationGradients(const tgt::ivec3& pos, const VolumeRAM* dataset, const tgt::vec3 spacing, Statistics& stats);

//---------------------------------------------------------------------------------------------------------------------

class VRN_CORE_API VolumeOperatorRegionGrowBase : public UnaryVolumeOperatorBase {
public:
    enum FloodFillMode { FLOODFILL_INTENSITY, FLOODFILL_GRADMAG, FLOODFILL_WEIGHTED };
    enum FloodFillNeighborhood { c6 = 1, c18 = 2, c26 = 3};

    virtual Volume* apply(const VolumeBase* volume, const tgt::ivec3& seed_pos, int segment, float lowerThreshold, float upperThreshold, float strictness, FloodFillMode mode, bool useThresholds, bool adaptive, float maxSeedDistance = 0, FloodFillNeighborhood nb = c26) const = 0;

protected:
    static const std::string loggerCat_;
};

// Generic implementation:
template<typename T>
class VolumeOperatorRegionGrowGeneric : public VolumeOperatorRegionGrowBase {
public:
    virtual Volume* apply(const VolumeBase* volume, const tgt::ivec3& seed_pos, int segment, float lowerThreshold, float upperThreshold, float strictness, FloodFillMode mode, bool useThresholds, bool adaptive, float maxSeedDistance = 0, FloodFillNeighborhood nb = c26) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorRegionGrowGeneric<T>::apply(const VolumeBase* vh, const tgt::ivec3& seed_pos, int segment, float lowerThreshold, float upperThreshold, float strictness, FloodFillMode mode, bool useThresholds, bool adaptive, float maxSeedDistance, FloodFillNeighborhood nb) const {
    const VolumeRAM* vol = vh->getRepresentation<VolumeRAM>();
    if(!vol)
        return 0;

    const VolumeAtomic<T>* dataset = dynamic_cast<const VolumeAtomic<T>*>(vol);
    if(!dataset)
        return 0;

    tgt::ivec3 dims = dataset->getDimensions();
    tgt::vec3 sp = vh->getSpacing();

    VolumeAtomic<bool> markedVoxels(dims);
    markedVoxels.clear();

    std::vector<tgt::ivec3> neighbors;
    for (int z=-1; z <= 1; z++)
        for (int y=-1; y <= 1; y++)
            for (int x=-1; x <= 1; x++) {
               tgt::ivec3 p(x, y, z);
                if (p !=tgt::ivec3::zero) {
                    if((abs(p.x) + abs(p.y) + abs(p.z)) <= nb)
                        neighbors.push_back(p);
                }
            }


    std::stack<tgt::ivec3> voxelStack;
    voxelStack.push(seed_pos);
    for (size_t i=0; i < neighbors.size(); i++)
        voxelStack.push(seed_pos + neighbors[i]);

    Statistics stats_value(false);
    Statistics stats_gradmag(false);

    float seed_value = dataset->getVoxelNormalized(seed_pos);
    float seed_stddev26 = neighborStandardDeviation(seed_pos, dataset, stats_value);

    float seed_gradmag = gradientMagnitude(seed_pos, dataset, sp);
    float seed_gradmag_stddev26 = neighborStandardDeviationGradients(seed_pos, dataset, sp, stats_gradmag);

    if (mode == FLOODFILL_INTENSITY || mode == FLOODFILL_WEIGHTED)
        LINFO("seed value: " << seed_value << ", " << "stddev 26: " << seed_stddev26);

    if (mode == FLOODFILL_GRADMAG || mode == FLOODFILL_WEIGHTED)
        LINFO("seed_gradmag: " << seed_gradmag << ", " << "stddev gradmag 26: " << seed_gradmag_stddev26);

    while (!voxelStack.empty()) {
        tgt::ivec3 pos = voxelStack.top();
        voxelStack.pop();

        if (pos.x < 2 || pos.x > dims.x - 3 ||
            pos.y < 2 || pos.y > dims.y - 3 ||
            pos.z < 2 || pos.z > dims.z - 3 ||
            markedVoxels.voxel(pos))
        {
            // on border or already visited
            continue;
        }

        float value = dataset->getVoxelNormalized(pos);

        if (adaptive /* && not in initial neighbor */) {
            stats_value.addSample(value);
            seed_stddev26 = stats_value.getStdDev();

            if (mode == FLOODFILL_GRADMAG || mode == FLOODFILL_WEIGHTED) {
                stats_gradmag.addSample(gradientMagnitude(pos, dataset, sp));
                seed_gradmag_stddev26 = stats_gradmag.getStdDev();
            }
        }


        if (useThresholds && (value < lowerThreshold || value > upperThreshold)) {
            // invalid value
            continue;
        }

        // Cost function: if less than 1 then voxel is within the region.
        //
        // Based on: Runzhen Huang, Kwan-Liu Ma. RGVis: Region growing based techniques for
        // volume visualization, 2003.
        float cost = 0.f;

        if (mode == FLOODFILL_INTENSITY)
            cost = fabs(value - seed_value) / (strictness * seed_stddev26);
        else if (mode == FLOODFILL_GRADMAG) {
            float gradmag = gradientMagnitude(pos, dataset, sp);
            cost = fabs(gradmag - seed_gradmag) / (strictness * seed_gradmag_stddev26);
        }
        else if (mode == FLOODFILL_WEIGHTED) {
            float cost_a = fabs(value - seed_value) / (strictness * seed_stddev26);

            float gradmag = gradientMagnitude(pos, dataset, sp);
            float cost_b = fabs(gradmag - seed_gradmag) / (strictness * seed_gradmag_stddev26);

            // weight p
            float p = (seed_gradmag_stddev26 / (seed_stddev26 + seed_gradmag_stddev26));
            cost = cost_a * p + cost_b * (1.f - p);
        }

        if (cost >= 1.f)
            continue;

        if (maxSeedDistance > 0 && tgt::distance(tgt::vec3(pos), tgt::vec3(seed_pos)) > maxSeedDistance)
            continue;

        // voxel is valid
        markedVoxels.voxel(pos) = true;

        // add neighbors to stack if not already visited
        for (size_t i=0; i < neighbors.size(); i++)
            if (!markedVoxels.voxel(pos + neighbors[i]))
                voxelStack.push(pos + neighbors[i]);
    }

    // now fill segmentation volume with all marked voxels
    float count = 0;
    VolumeRAM_UInt8* segvol = new VolumeRAM_UInt8(dataset->getDimensions());
    for (size_t z=0; z < markedVoxels.getDimensions().z; z++) {
        for (size_t y=0; y < markedVoxels.getDimensions().y; y++) {
            for (size_t x=0; x < markedVoxels.getDimensions().x; x++) {
                if (markedVoxels.voxel(x, y, z)) {
                    uint8_t& v = segvol->voxel(x, y, z);
                    if (v == 0 || segment == 0) {
                        v = segment;
                        count++;
                    }
                }
            }
        }
    }

    return new Volume(segvol, vh);
    //return static_cast<int>(count);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorRegionGrowBase> VolumeOperatorRegionGrow;

} // namespace

#endif
