/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "utils.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include <random>

namespace voreen {
namespace utils {

tgt::mat4 createTransformationMatrix(const tgt::vec3& position, const tgt::vec3& velocity) {

    tgt::vec3 tangent(tgt::normalize(velocity));

    tgt::vec3 temp(0.0f, 0.0f, 1.0f);
    if (1.0f - std::abs(tgt::dot(temp, tangent)) <= std::numeric_limits<float>::epsilon())
        temp = tgt::vec3(0.0f, 1.0f, 0.0f);

    tgt::vec3 binormal(tgt::normalize(tgt::cross(temp, tangent)));
    tgt::vec3 normal(tgt::normalize(tgt::cross(tangent, binormal)));

    return tgt::mat4(normal.x, binormal.x, tangent.x, position.x,
                     normal.y, binormal.y, tangent.y, position.y,
                     normal.z, binormal.z, tangent.z, position.z,
                     0.0f, 0.0f, 0.0f, 1.0f);
}

std::vector<tgt::vec3> sampleDisk(const VolumeBase* volume, const tgt::vec3& origin, const tgt::vec3& normal, float radius, bool transformSamples, size_t numSamples) {
    return sampleCylinder(volume, origin, normal, radius, 0.0f, transformSamples, numSamples);
}

std::vector<tgt::vec3> sampleCylinder(const VolumeBase* volume, const tgt::vec3& origin, const tgt::vec3& normal, float radius, float length, bool transformSamples, size_t numSamples) {

    using T = tgt::vec3; // Could be quite easy to template..

    if(!volume) {
        return std::vector<T>();
    }

    RealWorldMapping rwm = volume->getRealWorldMapping();
    tgt::mat4 indicatorSpaceMatrix = utils::createTransformationMatrix(origin, normal);
    tgt::mat4 worldToIndicatorSpaceMatrix = volume->getWorldToVoxelMatrix();
    worldToIndicatorSpaceMatrix *= indicatorSpaceMatrix;

    // Estimate number of samples
    if(numSamples == 0) {
        length = length > 0.0f ? length : tgt::length(volume->getSpacing());
        float voxelsPerLength = length / tgt::length(volume->getSpacing());
        float voxelsPerRadius = radius / tgt::length(volume->getSpacing());
        numSamples = 2.0f * tgt::PIf * voxelsPerRadius * voxelsPerRadius * voxelsPerLength; // Use twice as many samples as minimally required.
        if(numSamples < 1) {
            LWARNINGC("SampleCylinder", "radius might be too small for proper sampling");
            numSamples = 10;
        }
    }

    VolumeRAMRepresentationLock lock(volume);
    const size_t numChannels = std::min<size_t>(lock->getNumChannels(), T::size);

    std::function<T(const tgt::vec3&)> sample = [&lock, &rwm, numChannels](const tgt::vec3& pos) {
        T v = T::zero;
        for(size_t channel = 0; channel < numChannels; channel++) {
            v[channel] = rwm.normalizedToRealWorld(lock->getVoxelNormalizedLinear(pos, channel));
        }
        return v;
    };

    // Apply basis transform of vector values into disk space (pun intended).
    if(transformSamples) {
        if(lock->getNumChannels() == 3) {
            tgt::mat4 sampleTransformationMatrix;
            indicatorSpaceMatrix.getRotationalPart().invert(sampleTransformationMatrix);

            sample = [sampleTransformationMatrix, sample](const tgt::vec3& pos) {
                return sampleTransformationMatrix * sample(pos);
            };
        }
        else {
            LWARNINGC("SampleCylinder", "Basis transform only applicable for Vector3 volumes");
        }
    }

    // Set up random generator (predictable!).
    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(0)));

    // Sample the cross section:
    std::vector<T> samples;
    samples.reserve(numSamples);

    for(size_t i=0; i<numSamples; i++) {

        float r   = radius * std::sqrt(rnd());
        float phi = rnd() * 2.0f * tgt::PIf;
        float z   = rnd() * length - length * 0.5f;

        // Calculate sample point in voxel space.
        tgt::vec3 pos = worldToIndicatorSpaceMatrix * tgt::vec3(r * std::cos(phi), r * std::sin(phi), z);
        samples.push_back(sample(pos));
    }

    return samples;
}

}
}