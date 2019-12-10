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

#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"

namespace voreen {

tgt::vec3 VolumeOperatorGradient::calcGradientCentralDifferences(const VolumeRAM* input, const tgt::vec3& spacing, const tgt::svec3& pos) {
    float v0, v1, v2, v3, v4, v5;
    tgt::vec3 gradient;

    if (pos.x != input->getDimensions().x-1)
        v0 = input->getVoxelNormalized(pos + tgt::svec3(1, 0, 0));
    else
        v0 = 0;
    if (pos.y != input->getDimensions().y-1)
        v1 = input->getVoxelNormalized(pos + tgt::svec3(0, 1, 0));
    else
        v1 = 0;
    if (pos.z != input->getDimensions().z-1)
        v2 = input->getVoxelNormalized(pos + tgt::svec3(0, 0, 1));
    else
        v2 = 0;

    if (pos.x != 0)
        v3 = input->getVoxelNormalized(pos + tgt::svec3(-1, 0, 0));
    else
        v3 = 0;
    if (pos.y != 0)
        v4 = input->getVoxelNormalized(pos + tgt::svec3(0, -1, 0));
    else
        v4 = 0;
    if (pos.z != 0)
        v5 = input->getVoxelNormalized(pos + tgt::svec3(0, 0, -1));
    else
        v5 = 0;

    gradient = tgt::vec3((v3 - v0), (v4 - v1), (v5 - v2));
    gradient /= (spacing * 2.0f);

    return gradient;
}

tgt::vec3 VolumeOperatorGradient::calcGradientSobel(const VolumeRAM* input, const tgt::vec3& spacing, const tgt::ivec3& pos, const size_t channel) {
    tgt::vec3 gradient = tgt::vec3(0.f);

    if (pos.x >= 1 && pos.x < tgt::ivec3(input->getDimensions()).x-1 &&
        pos.y >= 1 && pos.y < tgt::ivec3(input->getDimensions()).y-1 &&
        pos.z >= 1 && pos.z < tgt::ivec3(input->getDimensions()).z-1)
    {
        //left plane
        float v000 = input->getVoxelNormalized(pos + tgt::ivec3(-1, -1, -1),channel);
        float v001 = input->getVoxelNormalized(pos + tgt::ivec3(-1, -1, 0),channel);
        float v002 = input->getVoxelNormalized(pos + tgt::ivec3(-1, -1, 1),channel);
        float v010 = input->getVoxelNormalized(pos + tgt::ivec3(-1, 0, -1),channel);
        float v011 = input->getVoxelNormalized(pos + tgt::ivec3(-1, 0, 0),channel);
        float v012 = input->getVoxelNormalized(pos + tgt::ivec3(-1, 0, 1),channel);
        float v020 = input->getVoxelNormalized(pos + tgt::ivec3(-1, 1, -1),channel);
        float v021 = input->getVoxelNormalized(pos + tgt::ivec3(-1, 1, 0),channel);
        float v022 = input->getVoxelNormalized(pos + tgt::ivec3(-1, 1, 1),channel);
        //mid plane
        float v100 = input->getVoxelNormalized(pos + tgt::ivec3(0, -1, -1),channel);
        float v101 = input->getVoxelNormalized(pos + tgt::ivec3(0, -1, 0),channel);
        float v102 = input->getVoxelNormalized(pos + tgt::ivec3(0, -1, 1),channel);
        float v110 = input->getVoxelNormalized(pos + tgt::ivec3(0, 0, -1),channel);
        //float v111 = input->getVoxelNormalized(pos + ivec3(0, 0, 0)); //not needed for calculation
        float v112 = input->getVoxelNormalized(pos + tgt::ivec3(0, 0, 1),channel);
        float v120 = input->getVoxelNormalized(pos + tgt::ivec3(0, -1, -1),channel);
        float v121 = input->getVoxelNormalized(pos + tgt::ivec3(0, -1, 0),channel);
        float v122 = input->getVoxelNormalized(pos + tgt::ivec3(0, -1, 1),channel);
        //right plane
        float v200 = input->getVoxelNormalized(pos + tgt::ivec3(1, -1, -1),channel);
        float v201 = input->getVoxelNormalized(pos + tgt::ivec3(1, -1, 0),channel);
        float v202 = input->getVoxelNormalized(pos + tgt::ivec3(1, -1, 1),channel);
        float v210 = input->getVoxelNormalized(pos + tgt::ivec3(1, 0, -1),channel);
        float v211 = input->getVoxelNormalized(pos + tgt::ivec3(1, 0, 0),channel);
        float v212 = input->getVoxelNormalized(pos + tgt::ivec3(1, 0, 1),channel);
        float v220 = input->getVoxelNormalized(pos + tgt::ivec3(1, 1, -1),channel);
        float v221 = input->getVoxelNormalized(pos + tgt::ivec3(1, 1, 0),channel);
        float v222 = input->getVoxelNormalized(pos + tgt::ivec3(1, 1, 1),channel);

        //filter x-direction
        gradient.x += -1 * v000;
        gradient.x += -2 * v010;
        gradient.x += -1 * v020;
        gradient.x += 1 * v200;
        gradient.x += 2 * v210;
        gradient.x += 1 * v220;
        gradient.x += -2 * v001;
        gradient.x += -4 * v011;
        gradient.x += -2 * v021;
        gradient.x += +2 * v201;
        gradient.x += +4 * v211;
        gradient.x += +2 * v221;
        gradient.x += -1 * v002;
        gradient.x += -2 * v012;
        gradient.x += -1 * v022;
        gradient.x += +1 * v202;
        gradient.x += +2 * v212;
        gradient.x += +1 * v222;

        //filter y-direction
        gradient.y += -1 * v000;
        gradient.y += -2 * v100;
        gradient.y += -1 * v200;
        gradient.y += +1 * v020;
        gradient.y += +2 * v120;
        gradient.y += +1 * v220;
        gradient.y += -2 * v001;
        gradient.y += -4 * v101;
        gradient.y += -2 * v201;
        gradient.y += +2 * v021;
        gradient.y += +4 * v121;
        gradient.y += +2 * v221;
        gradient.y += -1 * v002;
        gradient.y += -2 * v102;
        gradient.y += -1 * v202;
        gradient.y += +1 * v022;
        gradient.y += +2 * v122;
        gradient.y += +1 * v222;

        //filter z-direction
        gradient.z += -1 * v000;
        gradient.z += -2 * v100;
        gradient.z += -1 * v200;
        gradient.z += +1 * v002;
        gradient.z += +2 * v102;
        gradient.z += +1 * v202;
        gradient.z += -2 * v010;
        gradient.z += -4 * v110;
        gradient.z += -2 * v210;
        gradient.z += +2 * v012;
        gradient.z += +4 * v112;
        gradient.z += +2 * v212;
        gradient.z += -1 * v020;
        gradient.z += -2 * v120;
        gradient.z += -1 * v220;
        gradient.z += +1 * v022;
        gradient.z += +2 * v122;
        gradient.z += +1 * v222;

        gradient /= 16.f;   // sum of all positive weights
        gradient /= 2.f;    // this mask has a step length of 2 voxels
        gradient /= spacing;
        gradient *= -1.f;
    }

    return gradient;
}

} // namespace

