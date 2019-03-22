/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/texture.h"
#include "tgt/types.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/gpucapabilities.h"
#include "tgt/filesystem.h"

namespace tgt {

//------------------------------------------------------------------------------
// Constructor and Destructor
//------------------------------------------------------------------------------
Texture::Texture(const tgt::ivec3& dimensions, GLint format, GLint internalformat,
                GLenum dataType, Filter filter, Wrapping wrapping, GLubyte* data, bool ownership)
    : dimensions_(dimensions)
    , internalformat_(internalformat)
    , format_(format)
    , dataType_(dataType)
    , filter_(filter)
    , wrapping_(wrapping)
    , cpuTextureData_(data)
    , ownsCpuTextureData_(ownership)
    , useMipMap_(false) // is set by setFilter() method
{
    //generate missing members
    generateId();
    type_ = calcType(dimensions_);
    bpt_ = calcFormatBpt(format_, dataType_);
    //update texture settings
    setFilter(filter_);
    setWrapping(wrapping_);
}

Texture::~Texture()
    {
        if (id_)
            glDeleteTextures(1, &id_);

        if(ownsCpuTextureData_)
            delete[] cpuTextureData_;
    }

//--------------------------------------------------------------------------------------
// Getter and Setter
//--------------------------------------------------------------------------------------
size_t Texture::getSizeOnCPU() const {
    return bpt_ * hmul(dimensions_);
}

int Texture::getSizeOnGPU() const {
    int bpt = calcInternalFormatBpt(internalformat_);
    return bpt * hmul(dimensions_);
}

void Texture::updateDimensions(const ivec3& dimensions, bool reuploadTexture) {
    tgtAssert(type_ == Texture::calcType(dimensions),"Unsupported dimension changes");

    if(ownsCpuTextureData_)
        delete[] cpuTextureData_;
    cpuTextureData_ = 0;

    dimensions_ = dimensions;

    if(reuploadTexture)
        uploadTexture();
}

void Texture::setFilter(Filter filter) {
    filter_ = filter;
    useMipMap_ = (filter == MIPMAP);

    bind();

    switch(filter_) {
        case NEAREST:
            glTexParameteri(type_,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
            glTexParameteri(type_,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
            break;

        case LINEAR:
            glTexParameteri(type_, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(type_, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            break;

        case ANISOTROPIC:
            glTexParameterf(type_, GL_TEXTURE_MAX_ANISOTROPY_EXT, GpuCaps.getMaxTextureAnisotropy());

        case MIPMAP:
            glTexParameteri(type_,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
            glTexParameteri(type_,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_LINEAR);
            //glTexParameteri(type_, GL_SGIS_generate_mipmap, GL_TRUE);
            break;
    }
}

void Texture::setWrapping(Wrapping w) {
    wrapping_ = w;

    bind();

    GLint wrap = wrapping_;

    /* set wrapping for all dimensions */
    switch(type_) {
    case GL_TEXTURE_3D:
        glTexParameteri(type_, GL_TEXTURE_WRAP_R, wrap);
        //no break
    case GL_TEXTURE_2D:
        glTexParameteri(type_, GL_TEXTURE_WRAP_T, wrap);
        //no break
    default: //GL_TEXTURE_1D:
        glTexParameteri(type_, GL_TEXTURE_WRAP_S, wrap);
        break;
    }
}

void Texture::setSwizzle(const GLint* swizzleMask) {
    bind();
    glTexParameteriv(type_, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);
}

void Texture::setCpuTextureData(GLubyte* data, bool ownership) {
    if(ownsCpuTextureData_)
        delete[] cpuTextureData_;

    cpuTextureData_ = data;
    ownsCpuTextureData_ = ownership;
}

//--------------------------------------------------------------------------------------
// Wrapper of OpenGL functions
//--------------------------------------------------------------------------------------
void Texture::generateId() {
    id_ = 0;
    glGenTextures(1, &id_);
}

void Texture::bind() const {
    glBindTexture(type_ , id_);
}

void Texture::enable() const {
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    glEnable(type_);
#else
    ImmediateMode::TexturMode texmode = static_cast<ImmediateMode::TexturMode>(type_);
    IMode.setTextureMode(texmode);
#endif
}

void Texture::disable() const {
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    glDisable(type_);
#else
    IMode.setTextureMode(ImmediateMode::TEXNONE);
#endif
}

GLubyte* Texture::alloc(bool ownership) {
    if(ownsCpuTextureData_)
        delete[] cpuTextureData_;

    cpuTextureData_ = new GLubyte[getSizeOnCPU()];
    ownsCpuTextureData_ = ownership;

    return cpuTextureData_;
}

//--------------------------------------------------------------------------------------
// Upload and download texture
//--------------------------------------------------------------------------------------
void Texture::uploadTexture() {
    bind();

    switch(type_) {
        case GL_TEXTURE_1D:
            glTexImage1D(GL_TEXTURE_1D, 0, internalformat_, dimensions_.x, 0,
                         format_, dataType_, cpuTextureData_);
            break;

        case GL_TEXTURE_2D:
            if(format_ == GL_RGBA) {
                glTexImage2D(GL_TEXTURE_2D, 0, internalformat_, dimensions_.x, dimensions_.y, 0,
                         format_, dataType_, cpuTextureData_);
            } else {
                glPixelStorei(GL_UNPACK_ALIGNMENT, 1); //alignment should not modify the texture (see http://www.opengl.org/wiki/Common_Mistakes#Texture_upload_and_pixel_reads)
                glTexImage2D(GL_TEXTURE_2D, 0, internalformat_, dimensions_.x, dimensions_.y, 0,
                         format_, dataType_, cpuTextureData_);
                glPixelStorei(GL_UNPACK_ALIGNMENT, 4); //switch back to default
            }
            break;

        case GL_TEXTURE_3D:
            glTexImage3D(GL_TEXTURE_3D, 0, internalformat_,
                         dimensions_.x, dimensions_.y, dimensions_.z, 0,
                         format_, dataType_, cpuTextureData_);
            break;

        case GL_TEXTURE_2D_ARRAY:
            glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, internalformat_,
                         dimensions_.x, dimensions_.y, dimensions_.z, 0,
                         format_, dataType_, cpuTextureData_);
            break;

        case GL_TEXTURE_RECTANGLE:
            glTexImage2D(GL_TEXTURE_RECTANGLE, 0, internalformat_, dimensions_.x, dimensions_.y, 0,
                            format_, dataType_, cpuTextureData_);
            break;

    }

    if (cpuTextureData_ && useMipMap_) {
        glGenerateMipmap(type_);
    }
}

void Texture::downloadTexture() const {
    bind();

    if (cpuTextureData_ == 0)
        const_cast<Texture*>(this)->alloc(true);

    glGetTexImage(type_, 0, format_, dataType_, cpuTextureData_);
}

GLubyte* Texture::downloadTextureToBuffer() const {
    return downloadTextureToBuffer(format_, dataType_);
}

void Texture::downloadTextureToBuffer(GLubyte* pixels, size_t numBytesAllocated) const {
    downloadTextureToBuffer(format_, dataType_, pixels, numBytesAllocated);
}

void Texture::downloadTextureToBuffer(GLint format, GLenum dataType,
                                      GLubyte* pixels, size_t numBytesAllocated) const {
    bind();

    size_t arraySize = hmul(dimensions_) * calcFormatBpt(format, dataType);
    if(numBytesAllocated < arraySize) {
        LWARNINGC("tgt.texture", "downloadTextureToBuffer: allocated buffer is too small");
    }
    else {
        glGetTexImage(type_, 0, format, dataType, pixels);
    }
}

GLubyte* Texture::downloadTextureToBuffer(GLint format, GLenum dataType) const {
    bind();

    int arraySize = hmul(dimensions_) * calcFormatBpt(format, dataType);
    GLubyte* pixels = new GLubyte[arraySize];

    glGetTexImage(type_, 0, format, dataType, pixels);
    return pixels;
}

//-----------------------------------------------------------------------------------------
// Texel access (only if downloadTexture has been called or cpuTextureData_ is up to date)
//-----------------------------------------------------------------------------------------
tgt::Color Texture::texelAsFloat(size_t x, size_t y) const {
    tgt::Color ret = tgt::Color(0.0f);
    switch(format_) {
        case GL_RGBA:
            switch(dataType_) {
                case GL_UNSIGNED_BYTE: {
                    tgt::Vector4<uint8_t> t = texel< tgt::Vector4<uint8_t> >(x,y);
                    ret.x = (float )t.x / 0xFF;
                    ret.y = (float )t.y / 0xFF;
                    ret.z = (float )t.z / 0xFF;
                    ret.w = (float )t.w / 0xFF;
                    break;
                }
                case GL_UNSIGNED_SHORT: {
                    tgt::Vector4<uint16_t> t = texel< tgt::Vector4<uint16_t> >(x,y);
                    ret.x = (float )t.x / 0xFFFF;
                    ret.y = (float )t.y / 0xFFFF;
                    ret.z = (float )t.z / 0xFFFF;
                    ret.w = (float )t.w / 0xFFFF;
                    break;
                }
                case GL_FLOAT:
                    ret = texel<tgt::Color>(x,y);
                    break;
                default:
                    LWARNINGC("tgt.texture", "texelAsFloat: Unknown data type!");
            }
            break;
        case GL_RGB:
            switch(dataType_) {
                case GL_UNSIGNED_BYTE: {
                    tgt::Vector3<uint8_t> t = texel< tgt::Vector3<uint8_t> >(x,y);
                    ret.x = (float )t.x / 0xFF;
                    ret.y = (float )t.y / 0xFF;
                    ret.z = (float )t.z / 0xFF;
                    ret.w = 1.0f;
                    break;
                }
                case GL_UNSIGNED_SHORT: {
                    tgt::Vector3<uint16_t> t = texel< tgt::Vector3<uint16_t> >(x,y);
                    ret.x = (float )t.x / 0xFFFF;
                    ret.y = (float )t.y / 0xFFFF;
                    ret.z = (float )t.z / 0xFFFF;
                    ret.w = 1.0f;
                    break;
                }
                case GL_FLOAT: {
                    tgt::Vector3f t = texel<tgt::Vector3f>(x,y);
                    ret.x = t.x;
                    ret.y = t.y;
                    ret.z = t.z;
                    ret.w = 1.0f;
                    break;
                }
                default:
                    LWARNINGC("tgt.texture", "texelAsFloat: Unknown data type!");
            }
            break;
        case GL_RED:
            switch(dataType_) {
                case GL_UNSIGNED_BYTE: {
                    tgt::Vector3<uint8_t> t = tgt::vec3(texel<uint8_t>(x,y));
                    ret.x = (float )t.x / 0xFF;
                    ret.y = (float )t.y / 0xFF;
                    ret.z = (float )t.z / 0xFF;
                    ret.w = 1.0f;
                    break;
                                       }
                case GL_UNSIGNED_SHORT: {
                    tgt::Vector3<uint16_t> t = tgt::vec3(texel<uint16_t>(x,y));
                    ret.x = (float )t.x / 0xFFFF;
                    ret.y = (float )t.y / 0xFFFF;
                    ret.z = (float )t.z / 0xFFFF;
                    ret.w = 1.0f;
                    break;
                                        }
                case GL_FLOAT: {
                    tgt::Vector3f t = tgt::vec3(texel<GLfloat>(x,y));
                    ret.x = t.x;
                    ret.y = t.y;
                    ret.z = t.z;
                    ret.w = 1.0f;
                    break;
                }
                default:
                    LWARNINGC("tgt.texture", "texelAsFloat: Unknown data type!");
        }
        break;

        default:
            LWARNINGC("tgt.texture", "texelAsFloat: Unknown format!");
    }
    return ret;
}

void Texture::texelFromFloat(const tgt::Color& color, size_t x, size_t y) {
    switch(format_) {
        case GL_RGBA:
            switch(dataType_) {
                case GL_UNSIGNED_BYTE: {
                    auto& t = texel< tgt::Vector4<uint8_t> >(x,y);
                    t.x = static_cast<uint8_t>(color.x * 0xFF);
                    t.y = static_cast<uint8_t>(color.y * 0xFF);
                    t.z = static_cast<uint8_t>(color.z * 0xFF);
                    t.w = static_cast<uint8_t>(color.w * 0xFF);
                    break;
                }
                case GL_UNSIGNED_SHORT: {
                    auto& t = texel< tgt::Vector4<uint16_t> >(x,y);
                    t.x = static_cast<uint16_t>(color.x * 0xFFFF);
                    t.y = static_cast<uint16_t>(color.y * 0xFFFF);
                    t.z = static_cast<uint16_t>(color.z * 0xFFFF);
                    t.w = static_cast<uint16_t>(color.w * 0xFFFF);
                    break;
                }
                case GL_FLOAT:
                    texel<tgt::Color>(x,y) = color;
                    break;
                default:
                    LWARNINGC("tgt.texture", "texelAsFloat: Unknown data type!");
            }
            break;
        case GL_RGB:
            switch(dataType_) {
                case GL_UNSIGNED_BYTE: {
                    auto& t = texel< tgt::Vector3<uint8_t> >(x,y);
                    t.x = static_cast<uint8_t>(color.x * 0xFF);
                    t.y = static_cast<uint8_t>(color.y * 0xFF);
                    t.z = static_cast<uint8_t>(color.z * 0xFF);
                    break;
                }
                case GL_UNSIGNED_SHORT: {
                    auto& t = texel< tgt::Vector3<uint16_t> >(x,y);
                    t.x = static_cast<uint16_t>(color.x * 0xFFFF);
                    t.y = static_cast<uint16_t>(color.y * 0xFFFF);
                    t.z = static_cast<uint16_t>(color.z * 0xFFFF);
                    break;
                }
                case GL_FLOAT: {
                    auto& t = texel<tgt::Vector3f>(x,y);
                    t.x = color.x;
                    t.y = color.y;
                    t.z = color.z;
                    break;
                }
                default:
                    LWARNINGC("tgt.texture", "texelAsFloat: Unknown data type!");
            }
            break;
        case GL_RED:
            switch(dataType_) {
                case GL_UNSIGNED_BYTE: {
                    auto& t = texel<uint8_t>(x, y);
                    t = static_cast<uint8_t>(color.x * 0xFF);
                    break;
                }
                case GL_UNSIGNED_SHORT: {
                    auto& t = texel<uint16_t>(x,y);
                    t = static_cast<uint16_t>(color.x * 0xFFFF);
                    break;
                }
                case GL_FLOAT: {
                    texel<GLfloat>(x, y) = color.x;
                    break;
                }
                default:
                    LWARNINGC("tgt.texture", "texelAsFloat: Unknown data type!");
            }
            break;

        default:
            LWARNINGC("tgt.texture", "texelAsFloat: Unknown format!");
    }
}

//--------------------------------------------------------------------------------------
// Helper functions (static)
//--------------------------------------------------------------------------------------
int Texture::calcFormatBpt(GLint format, GLenum dataType) {

    int numChannels = calcNumChannels(format);

    int typeSize = 0;
    switch (dataType) {
        case GL_BYTE:
        case GL_UNSIGNED_BYTE:
            typeSize = 1;
            break;

        case GL_SHORT:
        case GL_UNSIGNED_SHORT:
            typeSize = 2;
            break;

        case GL_INT:
        case GL_UNSIGNED_INT:
        case GL_FLOAT:
        case GL_UNSIGNED_INT_24_8:
            typeSize = 4;
            break;
        case GL_FLOAT_32_UNSIGNED_INT_24_8_REV:
        case GL_DOUBLE:
            typeSize = 8;
            break;

        default:
            LWARNINGC("tgt.Texture", "unknown dataType");
            tgtAssert(false,"unknown dataType");
    }

    return typeSize * numChannels;
}

int Texture::calcInternalFormatBpt(GLint internalformat) {

    int bpt = 0;
    switch (internalformat) {
        case 1:
        case GL_RED:
        case GL_GREEN:
        case GL_BLUE:
        case GL_DEPTH_COMPONENT:
            bpt = 1;
            break;

        case 2:
        case GL_RG:
        case GL_DEPTH_COMPONENT16:
            bpt = 2;
            break;

        case GL_RGB:
        case GL_BGR:
        case GL_DEPTH_COMPONENT24:
            bpt = 3;
            break;

        case GL_RGBA:
        case GL_RGBA8:
        case GL_BGRA:
        case GL_DEPTH_COMPONENT32:
        case GL_DEPTH24_STENCIL8:
        case GL_R32F:
        case GL_R32UI:
        case GL_RG32F:
            bpt = 4;
            break;
#ifdef GL_DEPTH32F_STENCIL8
        case GL_DEPTH32F_STENCIL8:
            bpt = 5;
            break;
#endif
        case GL_RGB16:
        case GL_RGB16F:
            bpt = 6;
            break;

        case GL_RGBA16:
        case GL_RGBA16F:
            bpt = 8;
            break;

        case GL_RGB32F:
            bpt = 12;
            break;

        case GL_RGBA32F:
            bpt = 16;
            break;

        default:
            LWARNINGC("tgt.Texture", "unknown internal format");
            tgtAssert(false, "unknown internal format");
            break;
    }

    return bpt;
}

int Texture::calcNumChannels(GLint format) {

    switch (format) {
    case 1:
    case GL_RED:
    case GL_RED_INTEGER:
    case GL_R32F:
    case GL_GREEN:
    case GL_BLUE:
    case GL_DEPTH_COMPONENT:
    case GL_DEPTH_COMPONENT24:
    case GL_DEPTH_STENCIL:
        return 1;
        break;

    case 2:
    case GL_RG:
        return 2;
        break;

    case 3:
    case GL_RGB:
    case GL_BGR:
        return 3;
        break;

    case 4:
    case GL_RGBA:
    case GL_BGRA:
    case GL_RGBA16:
    case GL_RG32F:
        return 4;
        break;

    default:
        LWARNINGC("tgt.Texture", "unknown format");
        tgtAssert(false, "unknown format");
        return 0;
    }
}

GLenum Texture::calcType(const ivec3& dimensions) {
    GLenum type;
    if (dimensions.z == 1)    {
        if (dimensions.y == 1)
            type = GL_TEXTURE_1D;
        else
            type = GL_TEXTURE_2D;
    }
    else {
        type = GL_TEXTURE_3D;
    }

    return type;
}

} // namespace tgt
