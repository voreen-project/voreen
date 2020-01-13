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

#ifndef VRN_VOLUMECREATEBASE_H
#define VRN_VOLUMECREATEBASE_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class VRN_CORE_API VolumeCreateBase : public VolumeProcessor {
protected:
    void fillPlane(VolumeRAM_UInt8* vds, tgt::vec3 center, tgt::vec3 normal, uint8_t value);
    void fillCircle(VolumeRAM_UInt8* vds, tgt::vec3 center, float radius, uint8_t value);
    void fillOrientedCircle(VolumeRAM_UInt8* vds, tgt::vec3 center, tgt::vec3 normal, float radius, uint8_t value);
    void fillSphere(VolumeRAM_UInt8* vds, tgt::vec3 center, float radius, uint8_t value);
    void fillEllipsoid(VolumeRAM_UInt8* vds, tgt::vec3 center, tgt::vec3 radius, uint8_t value);
    void fillBox(VolumeRAM_UInt8* vds, tgt::ivec3 start, tgt::ivec3 end, uint8_t value);
    void fillOrientedBox(VolumeRAM_UInt8* vds, tgt::vec3 center, tgt::vec3 dir, float lengthA, float lengthB, float yStart, float yEnd, uint8_t value);
    void fillCone(VolumeRAM_UInt8* vds, tgt::vec3 center, tgt::vec3 dir, float length, float yStart, float yEnd, float topRadius, float bottomRadius, uint8_t value);

    void applyPerturbation(VolumeRAM* vds, tgt::ivec3 dimensions, tgt::vec3 frequency, tgt::vec3 amplitude);

    tgt::vec3 getRandVec() const;
    float getRandFloat() const;
};

} // namespace

#endif // VRN_VOLUMECREATEBASE_H
