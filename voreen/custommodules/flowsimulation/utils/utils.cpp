#include "utils.h"

#include "voreen/core/datastructures/volume/volumebase.h"

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

tgt::vec3 sampleDisk(const VolumeBase* volume, const tgt::vec3& origin, const tgt::vec3& normal, float radius, size_t div) {

    if(!volume) {
        return tgt::vec3::zero;
    }

    tgt::mat4 worldToIndicatorSpaceMatrix = volume->getWorldToVoxelMatrix();
    worldToIndicatorSpaceMatrix *= utils::createTransformationMatrix(origin, normal);

    VolumeRAMRepresentationLock representation(volume);

    tgt::vec3 meanVelocity = tgt::vec3::zero;
    float maxVelocityMagnitude = 0.0f;

    // Sample the cross section:
    const float da = 2.0f * tgt::PIf / div;
    const float dr = radius / div;
    for(size_t ai=0; ai<div; ai++) {
        float a = da * ai;
        for(size_t ri=0; ri<div; ri++) {
            float r = std::sqrt(dr * (ri + 1)); // Uniform disk sampling.

            // Calculate sample point in voxel space.
            tgt::vec3 voxel = worldToIndicatorSpaceMatrix * tgt::vec3(r * std::cos(a), r * std::sin(a), 0.0f);

            // Determine velocity at the calculated location.
            tgt::vec3 velocity = tgt::vec3::zero;
            for(size_t channel=0; channel < representation->getNumChannels(); channel++) {
                velocity[channel] = representation->getVoxelNormalizedLinear(voxel, channel);
            }

            // Gather statistics.
            float magnitude = tgt::length(velocity);
            if(magnitude > maxVelocityMagnitude) {
                maxVelocityMagnitude = magnitude;
            }

            meanVelocity += velocity;
        }
    }

    if (meanVelocity != tgt::vec3::zero) {
        meanVelocity = tgt::normalize(meanVelocity);
        meanVelocity *= maxVelocityMagnitude;
    }

    return meanVelocity;
}

}
}