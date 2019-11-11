#ifndef VRN_FLOW_SIMULATION_UTILS_H
#define VRN_FLOW_SIMULATION_UTILS_H

#include "tgt/matrix.h"
#include "tgt/vector.h"

namespace voreen {

class VolumeBase;

namespace utils {

tgt::mat4 createTransformationMatrix(const tgt::vec3& position, const tgt::vec3& velocity);

/**
 * Samples a disk inside a 3D-Vectorfield volume.
 * The returned vector is the normalized average of all sampled velocities
 * in the indicator space multiplied by the maximum magnitude.
 */
tgt::vec3 sampleDisk(const VolumeBase* volume, const tgt::vec3& origin, const tgt::vec3& normal, float radius, size_t div = 5);

}
}

#endif