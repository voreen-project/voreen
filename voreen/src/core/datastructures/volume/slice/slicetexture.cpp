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

#include "voreen/core/datastructures/volume/slice/slicetexture.h"


namespace voreen {

SliceTexture::SliceTexture(const tgt::ivec2& sliceDim, SliceAlignment alignment, const std::string& format, const std::string& baseType,
                           tgt::vec3 originInWorldSpace, tgt::vec3 xDirectionInWorldSpace, tgt::vec3 yDirectionInWorldSpace, const RealWorldMapping& rwm,
                           void* data, GLint textureFormat, GLint internalFormat, GLenum textureDataType)
    : tgt::Texture(tgt::ivec3(sliceDim,1),textureFormat,internalFormat,textureDataType,tgt::Texture::LINEAR,tgt::Texture::REPEAT,reinterpret_cast<GLubyte*>(data),true)
    , alignment_(alignment)
    , format_(format)
    , baseType_(baseType)
    , originInWorldSpace_(originInWorldSpace)
    , xDirectionInWorldSpace_(xDirectionInWorldSpace)
    , yDirectionInWorldSpace_(yDirectionInWorldSpace)
    , rwm_(rwm)
{
    tgtAssert(data, "no data passed");
    //upload the texture
    uploadTexture();
}

std::string SliceTexture::getFormat() const {
    return format_;
}

std::string SliceTexture::getBaseType() const {
    return baseType_;
}

SliceAlignment SliceTexture::getAlignment() const {
    return alignment_;
}

tgt::mat4 SliceTexture::getTextureToWorldMatrix() const {
    tgt::vec3 zVec = normalize(cross(xDirectionInWorldSpace_, yDirectionInWorldSpace_));

    tgt::mat4 m(xDirectionInWorldSpace_.x, yDirectionInWorldSpace_.x, zVec.x, originInWorldSpace_.x,
                xDirectionInWorldSpace_.y, yDirectionInWorldSpace_.y, zVec.y, originInWorldSpace_.y,
                xDirectionInWorldSpace_.z, yDirectionInWorldSpace_.z, zVec.z, originInWorldSpace_.z,
                                     0.0f,                      0.0f,   0.0f,                  1.0f);
    return m;
}

tgt::mat4 SliceTexture::getWorldToTextureMatrix() const {
    tgt::mat4 m = getTextureToWorldMatrix();
    tgt::mat4 inv;
    m.invert(inv);
    return inv;
}

voreen::RealWorldMapping SliceTexture::getRealWorldMapping() const {
    return rwm_;
}


} // namespace voreen
