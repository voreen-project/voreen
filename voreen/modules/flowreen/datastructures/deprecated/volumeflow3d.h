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

#ifndef VRN_FLOWVOLUME_H
#define VRN_FLOWVOLUME_H

#include "flow3d.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class VolumeFlow3D : public VolumeAtomic<tgt::vec3>
{
public:
    VolumeFlow3D(tgt::vec3* const voxels, const tgt::ivec3& dimensions,
        const float min = 0.0f, const float max = 1.0f, const float maxMagnitude = 1.0f);

    virtual ~VolumeFlow3D() {}

    const Flow3D& getFlow3D() const { return flow3D_; }

    float getMinValue() const { return flow3D_.minValue_; }
    float getMaxValue() const { return flow3D_.maxValue_; }
    float getMaxMagnitude() const { return flow3D_.maxMagnitude_; }

private:
    const Flow3D flow3D_;
};

}   // namespace

#endif  // VRN_FLOWVOLUME_H
