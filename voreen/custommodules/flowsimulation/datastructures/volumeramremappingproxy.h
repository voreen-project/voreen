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

#ifndef VRN_VOLUMERAMREMAPPINGPROXY_H
#define VRN_VOLUMERAMREMAPPINGPROXY_H

#include "include/voreen/core/datastructures/volume/volumeram.h"
#include "include/voreen/core/datastructures/volume/volume.h"

#include <functional>

namespace voreen {

class VRN_CORE_API VolumeRAMRemappingProxy : public VolumeRAM {
public:

    using RemappingFunction = std::function<bool(tgt::vec3&)>;

    // This function is called by all functions that are or will not be implemented for this proxy representation.
    void unimplemented() const;

    // This function is called by all functions that modify the volume.
    void readOnly() const;

    VolumeRAMRemappingProxy(tgt::svec3 resolution,
                            const VolumeBase* originalVolume,
                            RemappingFunction remappingFunction,
                            float missingValue = 0.0f);

    // Unimplemented:
    virtual VolumeRAM* clone() const;
    virtual VolumeRAM* clone(void* data) const;
    virtual VolumeRAM* createNew(const tgt::svec3& dimensions, bool allocMem = false) const;
    virtual VolumeRAM* getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset) const;

    // Read-only:
    virtual void clear();
    virtual void* getData();
    virtual void setVoxelNormalized(float value, const tgt::svec3& pos, size_t channel);
    virtual void setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel);
    virtual void setVoxelNormalized(float value, size_t index, size_t channel);

    // Delegated:
    virtual size_t getNumBytes() const;
    virtual size_t getNumChannels() const;
    virtual size_t getBitsAllocated() const;
    virtual size_t getBytesPerVoxel() const;
    virtual bool isSigned() const;
    virtual bool isInteger() const;
    virtual std::string getFormat() const;
    virtual std::string getBaseType() const;
    virtual tgt::vec2 elementRange() const;

    // Specialization:
    virtual const void* getData() const;
    virtual void* getBrickData(const tgt::svec3& offset, const tgt::svec3& dimension) const;
    virtual void* getSliceData(const size_t firstSlice, const size_t lastSlice) const;

    virtual float minNormalizedValue(size_t channel) const;
    virtual float maxNormalizedValue(size_t channel) const;
    virtual float minNormalizedMagnitude() const;
    virtual float maxNormalizedMagnitude() const;

    virtual void invalidate() const;

    // Mapping:
    virtual float getVoxelNormalized(const tgt::svec3& pos, size_t channel) const;
    virtual float getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel) const;
    virtual float getVoxelNormalized(size_t index, size_t channel) const;

protected:

    tgt::svec3 indexToPos(size_t index) const;
    float getVoxel(tgt::vec3 pos, size_t channel) const;
    VolumeRAM* getRawData() const;

    const VolumeBase* originalVolume_;
    VolumeRAMRepresentationLock data_;
    RemappingFunction remappingFunction_;
    float missingValue_;

    mutable std::unique_ptr<VolumeRAM> rawData_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMERAMREMAPPINGPROXY_H
