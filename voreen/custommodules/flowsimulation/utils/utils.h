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

#ifndef VRN_FLOW_SIMULATION_UTILS_H
#define VRN_FLOW_SIMULATION_UTILS_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumebase.h"

namespace voreen {

namespace utils {

/**
 * Creates a transformation matrix for transforming the specified position and velocity into
 * an axis aligned coordinate system where position points in x- and velocity in z-direction.
 */
tgt::mat4 createTransformationMatrix(const tgt::vec3& position, const tgt::vec3& velocity);

/**
 * Samples a disk randomly inside a 3D vector field volume.
 * Returns all sampled values.
 * This function can also be used for single channel volumes, however, transformSamples will have no effect.
 * @param volume The vector field volume
 * @param origin The disk's origin
 * @param normal The disk's normal
 * @param radius The disk's radius
 * @param transformSamples determines if the samples should be transformed into a plane defined by the disk
 *                         Note that that you can only rely on the through-plane component!
 * @param numSamples Number of samples to distribute uniformly across the disk surface.
 *                   If set to 0, it will be estimated automatically.
 */
std::vector<tgt::vec3> sampleDisk(const VolumeBase* volume, const tgt::vec3& origin, const tgt::vec3& normal, float radius, bool transformSamples = false, size_t numSamples = 0);

}

}

#endif