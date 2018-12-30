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

#ifndef VRN_VOLUMESLICEHELPER_H
#define VRN_VOLUMESLICEHELPER_H

#include "voreen/core/datastructures/volume/slice/slicetexture.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"

namespace voreen {

class TriangleMeshGeometryNormal;

/**
 * Helper class for the extraction of slices from a volume.
 */
class VRN_CORE_API SliceHelper {

public:

/**
 * Extracts an axis-aligned slice of the passed volume.
 * The caller takes ownership of the returned object.
 *
 * @param volumeRAM the volume to create the slice texture from.
 * @param alignment Slice direction to be extracted, must be axis-aligned.
 * @param sliceIndex Index of the slice to be extracted, must be smaller than volumeDim[alignment].
 * @param shiftArray Can be null-pointer. Should contain the shift of each channel in alignment direction.
 * @param levelOfDetail specifies the desired level of detail of the slice, with 0 being the full resolution.
 *        This parameter is currently only relevant, if the slice if constructed from an octree representation.
 * @param timeLimit maximum time in milliseconds the slice construction is allowed to take. After the time-limit
 *        has been reached, all subsequent octree brick accesses are omitted unless the bricks are already present in the CPU RAM.
 *        A value of 0 is interpreted as no time limit. This parameter is currently only relevant for octree representations.
 * @param complete out-parameter returning whether the slice could be constructed in full resolution (or the desired LOD)
 *        or whether the time limit has caused a resolution reduction. If the null pointer is passed, no value is returned.
 *
 * @see VolumeOctreeBase
 *
 * @return the created slice texture, whose dimensions and data type matches the input volume's properties.
 *         If the slice texture could not be created, 0 is returned.
 */
static SliceTexture* getVolumeSlice(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray = 0,
                                     size_t levelOfDetail = 0, clock_t timeLimit = 0, bool* complete = 0);

/**
 * Extracts an arbitrarily aligned slice from the passed volume.
 * TODO: more doc
 */
static SliceTexture* getVolumeSlice(const VolumeBase* volume, tgt::plane pl, float samplingRate);

/**
 * @brief Generates a geometry that represents slice number \p sliceIndex with orientation \p alignment through the volume \p vh
 *
 * @param vh The primary volume.
 * @param applyTransformation Apply the physicalToWorld-Matrix?
 * @param secondaryVolumes You can specify additional volumes to extend the area of the slice geometry to include these volumes. (For multi-volume slicing)
 */
static TriangleMeshGeometryNormal* getSliceGeometry(const VolumeBase* volume, SliceAlignment alignment, float sliceIndex, bool applyTransformation = true, const std::vector<const VolumeBase*> secondaryVolumes = std::vector<const VolumeBase*>());

protected:
/**
 * Extracts an axis-aligned 2D slice from the passed RAM volume and copies the pixel data to the passed buffer.
 * It is a helper function used by the SliceHelper.
 *
 * @param volumeRAM the volume to create the slice from.
 * @param alignment Slice direction to be extracted, must be axis-aligned.
 * @param sliceIndex Index of the slice to be extracted, must be smaller than volumeDim[alignment].
 * @param dataBuffer Out-parameter that will hold the pixel data buffer.
 *                   Note: The buffer will be allocated internally by the function!
 * @param textureFormat Out-parameter specifiying the OpenGL texture format to use when creating a tgt::Texture from the returned pixel data buffer.
 * @param internalFormat Out-parameter specifiying the OpenGL internal texture format to use when creating a tgt::Texture from the returned pixel data buffer.
 * @param textureDataType Out-parameter specifiying the OpenGL texture data type to use when creating a tgt::Texture from the returned pixel data buffer.
 *
 */
static void extractAlignedSlicePixelDataHelper(const VolumeRAM* volumeRAM, SliceAlignment alignment, size_t sliceIndex, int* shiftArray,
    void*& dataBuffer, GLint& textureFormat, GLint& internalFormat, GLenum& textureDataType);

    static const std::string loggerCat_;    ///< meow
};

} // namespace voreen

#endif
