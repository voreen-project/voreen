/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "colorpool.h"

#include "tgt/assert.h"

#include "../datastructures/ensembledataset.h"

namespace {
    tgt::vec3 generateDistinctColor(size_t index) {
        static std::vector<tgt::vec3> colors;
        if (colors.empty()) {
            colors.push_back(tgt::vec3(230,  25,  75));
            colors.push_back(tgt::vec3( 60, 180,  75));
            colors.push_back(tgt::vec3(255, 225,  25));
            colors.push_back(tgt::vec3(  0, 130, 200));
            colors.push_back(tgt::vec3(245, 130,  48));
            colors.push_back(tgt::vec3(145,  30, 180));
            colors.push_back(tgt::vec3( 70, 240, 240));
            colors.push_back(tgt::vec3(240,  50, 230));
            colors.push_back(tgt::vec3(210, 245,  60));
            colors.push_back(tgt::vec3(250, 190, 190));
            colors.push_back(tgt::vec3(  0, 128, 128));
            colors.push_back(tgt::vec3(230, 190, 255));
            colors.push_back(tgt::vec3(170, 110,  40));
            colors.push_back(tgt::vec3(255, 250, 200));
            colors.push_back(tgt::vec3(170, 255, 195));
            colors.push_back(tgt::vec3(128, 128,   0));
            colors.push_back(tgt::vec3(255, 215, 180));
            colors.push_back(tgt::vec3(  0,   0, 128));

            // Normalize to range [0,1]
            for (tgt::vec3& color : colors)
                color = color / tgt::vec3(255.0f);
        }

        if (index < colors.size())
            return colors[index];

        return tgt::vec3::one;
    }
}

namespace voreen {

const tgt::vec3 ColorPool::FADE_OUT_COLOR = tgt::vec3::one;

ColorPool::ColorPool(const EnsembleDataset& dataset)
    : dataset_(dataset)
{
}

ColorPool::~ColorPool() {
}

tgt::vec3 ColorPool::getColor(const std::string& run, size_t timeStep) {
    for(size_t i=0; i < dataset_.getRuns().size(); i++) {
        if(dataset_.getRuns()[i].name_ == run)
            return getColor(i, timeStep);
    }

    tgtAssert(false, "run not contained in ensemble");
    return tgt::vec3::zero;
}

tgt::vec3 ColorPool::getColor(size_t run, size_t timeStep) {
    tgtAssert(run < dataset_.getRuns().size(), "run index too large");
    size_t numTimeSteps = dataset_.getRuns()[run].timeSteps_.size();
    tgtAssert(timeStep == static_cast<size_t>(-1) || timeStep < numTimeSteps, "time step index too large");

    tgt::vec3 color = getDistinctColor(run);
    if(timeStep == static_cast<size_t>(-1))
        return color;

    float t = static_cast<float>(timeStep) / numTimeSteps;
    return (1.0f - t) * color + t * FADE_OUT_COLOR;
}

tgt::vec3 ColorPool::getDistinctColor(size_t index) {
    return generateDistinctColor(index);
}

}
