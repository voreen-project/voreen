/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_SLICETEXTURE_H
#define VRN_SLICETEXTURE_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h" //for the realworldmapping

#include "tgt/texture.h"
#include "tgt/matrix.h"

namespace voreen {

/**
 * Enum used to determine the alignment of a slice.
 */
enum SliceAlignment {
    YZ_PLANE = 0,
    XZ_PLANE = 1,
    XY_PLANE = 2,
    UNALIGNED_PLANE = 3,
};

/**
 * Voreen representation of a volume slice.
 * It is basically a tgt::Texture with meta informations like a real world mapping.
 */
class VRN_CORE_API SliceTexture : public tgt::Texture {
public:
    /** Constructor */
    SliceTexture(const tgt::ivec2& sliceDim, SliceAlignment alignment, const std::string& format, const std::string& baseType,
                 tgt::vec3 originInWorldSpace, tgt::vec3 xDirectionInWorldSpace, tgt::vec3 yDirectionInWorldSpace, RealWorldMapping rwm,
                 void* data, GLint textureFormat, GLint internalFormat, GLenum textureDataType);
    /** Destructor */
    virtual ~SliceTexture() {};

    //------------------------------
    //  Getter
    //------------------------------
    /** Returns the alignment according to the original volume. */
    SliceAlignment getAlignment() const;
    /** Returns the Voreen volume format. */
    std::string getFormat() const;
    /** Returns the Voreen volume base type. */
    std::string getBaseType() const;
    /** Returns the dimensions of the slice. */
    tgt::ivec2 getSliceDimensions() const;
    /** Returns the texture to world matrix. */
    tgt::mat4 getTextureToWorldMatrix() const;
    /** Returns the world to texture matrix. */
    tgt::mat4 getWorldToTextureMatrix() const;
    /** Returns the real world mapping. */
    RealWorldMapping getRealWorldMapping() const;

private:
    SliceAlignment alignment_;  ///< alignment of the slice texture

    std::string format_;        ///< Voreen volume format (see VolumeFactory)
    std::string baseType_;      ///< Voreen base type (e.g., "float" for format "Vector3(float)")

    tgt::vec3 originInWorldSpace_;      ///< slice position (0,0) in world space
    tgt::vec3 xDirectionInWorldSpace_;  ///< slice x direction in world space. Is dependent of the alignment
    tgt::vec3 yDirectionInWorldSpace_;  ///< slice y direction in world space. Is dependent of the alignment

    RealWorldMapping rwm_;  ///< the real world mapping
};

} // namespace voreen

#endif
