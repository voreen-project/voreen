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

#ifndef VRN_VOLUMEOPERATORCONVERT_H
#define VRN_VOLUMEOPERATORCONVERT_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/utils/exception.h"
#include "tgt/vector.h"


namespace voreen {

/**
 * Converts the source volume's data to the destination volume's data type
 * and writes the result to the dest volume.
 *
 * The source volume has to be passed to the constructor,
 * the dest volume has to be passed to apply().
 *
 * @note The dest volume has to be instantiated by the caller
 *       and must match the source volume in dimension and channel count.
 */
class VRN_CORE_API VolumeOperatorConvert  {
public:
    /**
     * @param srcVolume The volume to be converted. Must not be null.
     */
    VolumeOperatorConvert() :
        progressReporter_(0),
        inputIntensityRange_(-1.0)
    {}

    /**
     * Sets the intensity range to be used for normalizing the
     * input volume's intensity values during the conversion.
     * If no range is set, the input volume's min/max intensities
     * are used.
     *
     * @param range the normalization range that is mapped to the
     *  data range of the output volume. Precondition: range.x < range.y
     *
     * @note This setting does only affect the conversion of unbounded
     *  data types (float, double).
     */
    void setInputIntensityRange(const tgt::dvec2& range) {
        if (range.x >= range.y)
            throw VoreenException("Illegal range parameter. Expecting range.x < range.y");
        inputIntensityRange_ = range;
    }

    /**
     * Performs the conversion.
     *
     * Returns a Volume with a VolumeAtomic<T>
     */
    template<class T>
    Volume* apply(const VolumeBase* srcVolume) const;

    /**
     * Assigns a progress bar that should be used by the
     * operator for indicating progress.
     */
    void setProgressReporter(ProgressReporter* progress) {
        progressReporter_ = progress;
    }

protected:
    ProgressReporter* progressReporter_;  ///< to be used by concrete subclasses for indicating progress
private:
    tgt::dvec2 inputIntensityRange_;
};

template<class T>
Volume* VolumeOperatorConvert::apply(const VolumeBase* srcVolumeHandle) const {
    const VolumeRAM* srcVolume = srcVolumeHandle->getRepresentation<VolumeRAM>();
    VolumeAtomic<T>* destVolume = new VolumeAtomic<T>(srcVolume->getDimensions());

    if (!srcVolume)
        throw VoreenException("VolumeOperatorConvert: source volume is null pointer");
    if (!destVolume)
        throw VoreenException("VolumeOperatorConvert: dest volume is null pointer");

    if (destVolume->getDimensions() != srcVolume->getDimensions())
        throw VoreenException("VolumeOperatorConvert: volume dimensions must match");
    if (destVolume->getNumChannels() != srcVolume->getNumChannels())
        throw VoreenException("VolumeOperatorConvert: number of channels must match");

    // check the source volume's type
    const VolumeRAM_UInt8* srcUInt8 = dynamic_cast<const VolumeRAM_UInt8*>(srcVolume);
    const VolumeRAM_UInt16* srcUInt16 = dynamic_cast<const VolumeRAM_UInt16*>(srcVolume);

    // check the dest volume's type
    VolumeRAM_UInt8* destUInt8 = dynamic_cast<VolumeRAM_UInt8*>(destVolume);
    VolumeRAM_UInt16* destUInt16 = dynamic_cast<VolumeRAM_UInt16*>(destVolume);

    RealWorldMapping destMapping;
    if (srcUInt8 && destUInt8) {
        LINFOC("voreen.VolumeOperatorConvert" ,"No conversion necessary: source and dest type equal (VolumeRAM_UInt8)");
        VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::svec3(0, 0, 0), srcUInt8->getDimensions(), progressReporter_)
            destUInt8->voxel(i) = srcUInt8->voxel(i);
        destMapping = srcVolumeHandle->getRealWorldMapping();
    }
    else if (srcUInt16 && destUInt16) {
        LINFOC("voreen.VolumeOperatorConvert" ,"No conversion necessary: source and dest type equal (VolumeRAM_UInt16)");
        VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::svec3(0, 0, 0), srcUInt16->getDimensions(), progressReporter_)
            destUInt16->voxel(i) = srcUInt16->voxel(i);
        destMapping = srcVolumeHandle->getRealWorldMapping();
    }
    else if (srcUInt16 && destUInt8) {
        LINFOC("voreen.VolumeOperatorConvert", "Using accelerated conversion from VolumeRAM_UInt16 -> VolumeRAM_UInt8");
        VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::svec3(0, 0, 0), srcUInt16->getDimensions(), progressReporter_)
            destUInt8->voxel(i) = srcUInt16->voxel(i) >> 8;
        destMapping = srcVolumeHandle->getRealWorldMapping();
    }
    else if (srcUInt8 && destUInt16) {
        LINFOC("voreen.VolumeOperatorConvert", "Using accelerated conversion from VolumeRAM_UInt8 -> VolumeRAM_UInt16");
        VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::svec3(0, 0, 0), srcUInt8->getDimensions(), progressReporter_)
            destUInt16->voxel(i) = srcUInt8->voxel(i) << 8;
        destMapping = srcVolumeHandle->getRealWorldMapping();
    }
    else {
        // fallback using getVoxelFloat/setVoxelFloat

        // map input range to output range
        float scale = 1.f;
        float offset = 0.f;
        if (srcVolume->isInteger() && destVolume->isInteger()) {
            // for int-to-uint or uint-to-int conversion, the normalized value
            // has to be mapped from [-1.0;1.0] to [0.0;1.0] or vice-versa, respectively
            if (srcVolume->isSigned() && !destVolume->isSigned()) {
                scale = 0.5f;
                offset = 0.5f;
            }
            else if (!srcVolume->isSigned() && destVolume->isSigned()) {
                scale = 2.f;
                offset = -1.f;
            }
        }
        else if (!srcVolume->isInteger()) { // float/double input volume
            // determine input mapping range
            float min, max;
            if (inputIntensityRange_ == tgt::dvec2(-1.0)) {
                // conversion range not set => use input volume's intensity range
                min = srcVolume->minNormalizedValue(0);
                max = srcVolume->maxNormalizedValue(0);
                for (size_t i=1; i<srcVolume->getNumChannels(); i++) {
                    min = std::min(min, srcVolume->minNormalizedValue(i));
                    max = std::max(max, srcVolume->maxNormalizedValue(i));
                }
            }
            else {
                // use assigned conversion range
                tgtAssert(inputIntensityRange_.x < inputIntensityRange_.y, "invalid intensity range");
                min = static_cast<float>(inputIntensityRange_.x);
                max = static_cast<float>(inputIntensityRange_.y);
            }

            // map [min:max] to output range
            if (destVolume->isInteger() && destVolume->isSigned()) {
                // map [min:max] to [-1.0:1.0]
                scale = 1.f / (max - min + 1.f);
                offset = -min * scale - 1.f;
            }
            else {
                // map [min:max] to [0.0:1.0]
                scale = 1.f / (max - min);
                offset = -min * scale;
            }
        }
        else {
            // id mapping for all other cases
            scale = 1.f;
            offset = 0.f;
        }

        // differentiate single-channel from multi-channel volumes
        if (srcVolume->getNumChannels() == 1) {
            LINFOC("voreen.VolumeOperatorConvert", "Using fallback with setVoxelNormalized and getVoxelNormalized (single-channel)");
            VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::ivec3(0, 0, 0), srcVolume->getDimensions(), progressReporter_)
                destVolume->setVoxelNormalized(srcVolume->getVoxelNormalized(i)*scale + offset, i);
        }
        else {
            tgtAssert(srcVolume->getNumChannels() == destVolume->getNumChannels(), "channel-count mis-match");
            size_t numChannels = srcVolume->getNumChannels();
            LINFOC("voreen.VolumeOperatorConvert", "Using fallback with setVoxelNormalized and getVoxelNormalized (" << numChannels << " channels)");
            VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::ivec3(0, 0, 0), srcVolume->getDimensions(), progressReporter_) {
                for (size_t channel=0; channel < numChannels; channel++)
                    destVolume->setVoxelNormalized(srcVolume->getVoxelNormalized(i, channel)*scale + offset, i, channel);
            }
        }

        // real world mapping has to revert the applied value transformation
        RealWorldMapping srcMapping = srcVolumeHandle->getRealWorldMapping();
        destMapping = RealWorldMapping::combine(RealWorldMapping(1.f/scale, -offset/scale, ""), srcVolumeHandle->getRealWorldMapping());
    }

    if (progressReporter_)
        progressReporter_->setProgress(1.f);

    Volume* destVolumeHandle = new Volume(destVolume, srcVolumeHandle);
    destVolumeHandle->setRealWorldMapping(destMapping);
    return destVolumeHandle;
}

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
