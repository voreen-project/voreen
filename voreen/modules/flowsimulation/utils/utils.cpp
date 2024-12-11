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

std::vector<tgt::vec3> sampleDisk(const VolumeBase* volume, const tgt::vec3& origin, const tgt::vec3& normal, float radius, bool transformSamples, size_t numSamples, const VolumeBase* sampleMask) {
    return sampleCylinder(volume, origin, normal, radius, 0.0f, transformSamples, numSamples, sampleMask);
}

std::vector<tgt::vec3> sampleCylinder(const VolumeBase* volume, const tgt::vec3& origin, const tgt::vec3& normal, float radius, float length, bool transformSamples, size_t numSamples, const VolumeBase* sampleMask) {

    using T = tgt::vec3; // Could be quite easy to template..

    if(!volume) {
        return std::vector<T>();
    }

    RealWorldMapping rwm = volume->getRealWorldMapping();
    tgt::mat4 toIndicatorSpaceMatrix = utils::createTransformationMatrix(origin, normal);
    tgt::mat4 worldToIndicatorSpaceMatrix = volume->getWorldToVoxelMatrix();
    worldToIndicatorSpaceMatrix *= toIndicatorSpaceMatrix;

    tgt::mat4 fromIndicatorSpaceMatrix;
    if(!toIndicatorSpaceMatrix.invert(fromIndicatorSpaceMatrix)) {
        return std::vector<T>();
    }

    // Estimate number of samples
    if(numSamples <= 0) {
        float diag = tgt::length(volume->getSpacing());
        length = length > 0.0f ? length : diag;
        float voxelsPerLength = length / diag;
        float voxelsPerRadius = radius / diag;
        numSamples = 2.0f * tgt::PIf * voxelsPerRadius * voxelsPerRadius * voxelsPerLength; // Use twice as many samples as minimally required.
        if(numSamples < 1) {
            LWARNINGC("SampleCylinder", "radius might be too small for proper sampling");
            numSamples = 10;
        }
    }

    VolumeRAMRepresentationLock volumeLock(volume);
    const size_t numChannels = std::min<size_t>(volumeLock->getNumChannels(), T::size);

    std::function<T(const tgt::vec3&)> sample = [&volumeLock, &rwm, numChannels](const tgt::vec3& pos) {
        T v = T::zero;
        for(size_t channel = 0; channel < numChannels; channel++) {
            v[channel] = rwm.normalizedToRealWorld(volumeLock->getVoxelNormalizedLinear(pos, channel));
        }
        return v;
    };

    // Apply basis transform of vector values into disk space (pun intended).
    if(transformSamples) {
        if(volumeLock->getNumChannels() == 3) {
            auto sampleTransformationMatrix = fromIndicatorSpaceMatrix.getRotationalPartMat3();

            sample = [sampleTransformationMatrix, sample](const tgt::vec3& pos) {
                return sampleTransformationMatrix * sample(pos);
            };
        }
        else {
            LWARNINGC("SampleCylinder", "Basis transform only applicable for Vector3 volumes");
        }
    }

    // Set up sample mask.
    std::vector<tgt::vec3> samplePositions;

    if(sampleMask) {
        tgt::mat4 maskSpaceToWorldMatrix = sampleMask->getVoxelToWorldMatrix();
        VolumeRAMRepresentationLock maskLock(sampleMask);
        auto dim = maskLock->getDimensions();

        auto insideCylinder = [&](tgt::vec3 point) {
            float halfLength = length * 0.5f;
            if(point.z < -halfLength || point.z > halfLength) {
                return false;
            }

            float dx = point.x - origin.x;
            float dy = point.y - origin.y;
            if(dx * dx + dy * dy > radius * radius) {
                return false;
            }

            return true;
        };

        for(size_t z=0; z<dim.z; z++) {
            for(size_t y=0; y<dim.y;  y++) {
                for(size_t x=0; x<dim.x; x++) {
                    tgt::vec3 pos = tgt::vec3(x, y, z);
                    if(maskLock->getVoxelNormalizedLinear(pos, 0) > 0.0f) {
                        auto worldPos = maskSpaceToWorldMatrix * pos;
                        if(insideCylinder(worldPos)) {
                            samplePositions.push_back(worldPos);
                        }
                    }
                }
            }
        }

        if(samplePositions.empty()) {
            //LWARNINGC("SampleCylinder", "No valid samples found in mask volume");
            return std::vector<T>();
        }
    }
    else {
        // Set up random generator (predictable!).
        std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(0)));

        // Sample the cylinder:
        samplePositions.reserve(numSamples);

        for(size_t i=0; i<numSamples; i++) {
            // Generate a random position.
            float r   = radius * std::sqrt(rnd());
            float phi = rnd() * 2.0f * tgt::PIf;
            float z   = rnd() * length - length * 0.5f;

            tgt::vec3 pos(r * std::cos(phi), r * std::sin(phi), z);
            samplePositions.push_back(pos);
        }
    }

    // Sample the cross-section:
    std::vector<T> samples;
    samples.reserve(numSamples);

    for(const tgt::vec3& samplePos : samplePositions) {
        tgt::vec3 pos = worldToIndicatorSpaceMatrix * samplePos;
        samples.push_back(sample(pos));
    }

    return samples;
}


std::vector<tgt::vec3> sampleSphere(const VolumeBase* volume, const tgt::vec3& origin, float radius, size_t numSamples, const VolumeBase* sampleMask) {

    using T = tgt::vec3; // Could be quite easy to template..

    if(!volume) {
        return std::vector<T>();
    }

    RealWorldMapping rwm = volume->getRealWorldMapping();
    tgt::mat4 worldToVoxelMatrix = volume->getWorldToVoxelMatrix();

    // Estimate number of samples
    if(numSamples <= 0) {
        float diag = tgt::length(volume->getSpacing());
        float voxelsPerRadius = radius / diag;
        numSamples = 2.0f * 4.0f / 3.0f * tgt::PIf * voxelsPerRadius * voxelsPerRadius; // Use twice as many samples as minimally required.
        if(numSamples < 1) {
            LWARNINGC("SampleSphere", "radius might be too small for proper sampling");
            numSamples = 10;
        }
    }

    VolumeRAMRepresentationLock volumeLock(volume);
    const size_t numChannels = std::min<size_t>(volumeLock->getNumChannels(), T::size);

    std::function<T(const tgt::vec3&)> sample = [&volumeLock, &rwm, numChannels](const tgt::vec3& pos) {
        T v = T::zero;
        for(size_t channel = 0; channel < numChannels; channel++) {
            v[channel] = rwm.normalizedToRealWorld(volumeLock->getVoxelNormalizedLinear(pos, channel));
        }
        return v;
    };

    // Set up sample mask.
    std::vector<tgt::vec3> samplePositions;

    if(sampleMask) {
        tgt::mat4 maskSpaceToWorldMatrix = sampleMask->getVoxelToWorldMatrix();
        VolumeRAMRepresentationLock maskLock(sampleMask);
        auto dim = maskLock->getDimensions();

        auto insideSphere = [&](tgt::vec3 point) {
            return tgt::distanceSq(point, origin) < radius * radius;
        };

        for(size_t z=0; z<dim.z; z++) {
            for(size_t y=0; y<dim.y;  y++) {
                for(size_t x=0; x<dim.x; x++) {
                    tgt::vec3 pos = tgt::vec3(x, y, z);
                    if(maskLock->getVoxelNormalizedLinear(pos, 0) > 0.0f) {
                        auto worldPos = maskSpaceToWorldMatrix * pos;
                        if(insideSphere(worldPos)) {
                            samplePositions.push_back(worldPos);
                        }
                    }
                }
            }
        }

        if(samplePositions.empty()) {
            //LWARNINGC("SampleSphere", "No valid samples found in mask volume");
            return std::vector<T>();
        }
    }
    else {
        // Set up random generator (predictable!).
        std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(0)));

        // Sample the cylinder:
        samplePositions.reserve(numSamples);

        auto bounds = volume->getBoundingBox().getBoundingBox();

        for(size_t i=0; i<numSamples; i++) {
            // Generate a random position.
            float r     = radius * std::sqrt(rnd());
            float phi   = rnd() * 2.0f * tgt::PIf;
            float theta = rnd() * tgt::PIf;

            tgt::vec3 pos(r * std::sin(theta) * std::cos(phi), r * std::sin(theta) * std::sin(phi), r * std::cos(theta));
            pos += origin; // Transform to world position.

            if(bounds.containsPoint(pos)) {
                samplePositions.push_back(pos);
            }
        }
    }

    // Sample the sphere:
    std::vector<T> samples;
    samples.reserve(numSamples);

    for(const tgt::vec3& samplePos : samplePositions) {
        tgt::vec3 pos = worldToVoxelMatrix * samplePos;
        samples.push_back(sample(pos));
    }

    return samples;
}

}
}