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

#include "voreen/core/datastructures/volume/volumegl.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "tgt/glcontextmanager.h"
#include "tgt/gpucapabilities.h"

#include <algorithm>
#include <typeinfo>
#include <sstream>

using tgt::vec3;
using tgt::bvec3;
using tgt::ivec3;
using tgt::svec3;
using tgt::col4;
using tgt::mat4;

namespace {

// calculates next power of two if it is not already power of two
static int fitPowerOfTwo(int i) {
    int result = 1;

    while (result < i)
        result <<= 1;

    return result;
}

// same as above but for each component of a ivec3
static ivec3 fitPowerOfTwo(const ivec3& v) {
    ivec3 result;

    result.x = fitPowerOfTwo(v.x);
    result.y = fitPowerOfTwo(v.y);
    result.z = fitPowerOfTwo(v.z);

    return result;
}

// returns true, if dims.x, dims.y and dims.z are power of two
bool isPowerOfTwo(const ivec3& dims) {
    return dims == fitPowerOfTwo(dims);
}

} // namespace

namespace voreen {

const std::string VolumeGL::loggerCat_("voreen.VolumeGL");

VolumeGL::VolumeGL(const VolumeRAM* volume) 
  : VolumeRepresentation(volume->getDimensions())
  , texture_(0)
{
    tgtAssert(volume, "No volume");
    VolumeFactory vf;
    format_ = vf.getFormat(volume);
    baseType_ = vf.getBaseType(format_);
    generateTexture(volume);
}

VolumeGL::~VolumeGL() {
    destroy();
}

void VolumeGL::destroy() {
    delete texture_;
}

const VolumeTexture* VolumeGL::getTexture() const {
    return texture_;
}

VolumeTexture* VolumeGL::getTexture() {
    return texture_;
}

void VolumeGL::generateTexture(const VolumeRAM* volume) {
    if (!GpuCaps.is3DTexturingSupported()) {
        std::string message = "3D textures apparently not supported by the OpenGL driver";
        LERROR(message);
        throw VoreenException(message);
    }

    if (!GpuCaps.isNpotSupported() && !isPowerOfTwo(getDimensions())) {
        std::string message = "Non-power-of-two textures apparently not supported by the OpenGL driver";
        LERROR(message);
        throw VoreenException(message);
    }

    //
    // Determine GL format
    //
    GLint format = 0;         ///< The format of textures which will are created.
    GLint internalFormat = 0; ///< The internal format of the textures which are created.
    GLenum dataType = 0;      ///< The data type of the textures which are created.
    // scalar
    if (dynamic_cast<const VolumeAtomic<uint8_t>*>(volume)) {
        format = GL_RED;
        internalFormat = GL_R8;
        dataType = GL_UNSIGNED_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<int8_t>*>(volume)) {
        format = GL_RED;
        internalFormat = GL_R8_SNORM;
        dataType = GL_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(volume)) {
        format = GL_RED;
        internalFormat = GL_R16;
        dataType = GL_UNSIGNED_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<int16_t>*>(volume)) {
        format = GL_RED;
        internalFormat = GL_R16_SNORM;
        dataType = GL_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<uint32_t>*>(volume)) {
        format = GL_RED;
        internalFormat = GL_R32F;
        dataType = GL_UNSIGNED_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<int32_t>*>(volume)) {
        format = GL_RED;
        internalFormat = GL_R32F;
        dataType = GL_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<uint64_t>*>(volume)) {
        std::string message = "VolumeRAM_UInt64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<int64_t>*>(volume)) {
        std::string message = "VolumeRAM_Int64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<float>*>(volume)) {
        format = GL_RED;
        internalFormat = GL_R32F;
        dataType = GL_FLOAT;
    }
    else if (dynamic_cast<const VolumeAtomic<double>*>(volume)) {
        std::string message = "VolumeRAM_Double not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    // vec2
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint8_t> >*>(volume)) {
        format = GL_RG;
        internalFormat = GL_RG8;
        dataType = GL_UNSIGNED_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int8_t> >*>(volume)) {
        format = GL_RG;
        internalFormat = GL_RG8_SNORM;
        dataType = GL_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint16_t> >*>(volume)) {
        format = GL_RG;
        internalFormat = GL_RG16;
        dataType = GL_UNSIGNED_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int16_t> >*>(volume)) {
        format = GL_RG;
        internalFormat = GL_RG16_SNORM;
        dataType = GL_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint32_t> >*>(volume)) {
        format = GL_RG;
        internalFormat = GL_RG32F;
        dataType = GL_UNSIGNED_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int32_t> >*>(volume)) {
        format = GL_RG;
        internalFormat = GL_RG32F;
        dataType = GL_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint64_t> >*>(volume)) {
        std::string message = "VolumeRAM_2xUInt64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int64_t> >*>(volume)) {
        std::string message = "VolumeRAM_2xInt64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<float> >*>(volume)) {
        format = GL_RG;
        internalFormat = GL_RG32F;
        dataType = GL_FLOAT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<double> >*>(volume)) {
        std::string message = "VolumeRAM_2xDouble not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    // vec3
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(volume)) {
        format = GL_RGB;
        internalFormat = GL_RGB8;
        dataType = GL_UNSIGNED_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int8_t> >*>(volume)) {
        format = GL_RGB;
        internalFormat = GL_RGB8_SNORM;
        dataType = GL_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(volume)) {
        format = GL_RGB;
        internalFormat = GL_RGB16;
        dataType = GL_UNSIGNED_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int16_t> >*>(volume)) {
        format = GL_RGB;
        internalFormat = GL_RGB16_SNORM;
        dataType = GL_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint32_t> >*>(volume)) {
        format = GL_RGB;
        internalFormat = GL_RGB;
        dataType = GL_UNSIGNED_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int32_t> >*>(volume)) {
        format = GL_RGB;
        internalFormat = GL_RGB32F;
        dataType = GL_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint64_t> >*>(volume)) {
        std::string message = "VolumeRAM_3xUInt64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int64_t> >*>(volume)) {
        std::string message = "VolumeRAM_3xInt64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<float> >*>(volume)) {
        format = GL_RGB;
        internalFormat = GL_RGB32F;
        dataType = GL_FLOAT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<double> >*>(volume)) {
        std::string message = "VolumeRAM_3xDouble not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    // vec4
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(volume)) {
        format = GL_RGBA;
        internalFormat = GL_RGBA8;
        dataType = GL_UNSIGNED_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int8_t> >*>(volume)) {
        format = GL_RGBA;
        internalFormat = GL_RGBA8_SNORM;
        dataType = GL_BYTE;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint16_t> >*>(volume)) {
        format = GL_RGBA;
        internalFormat = GL_RGBA16;
        dataType = GL_UNSIGNED_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int16_t> >*>(volume)) {
        format = GL_RGBA;
        internalFormat = GL_RGBA16_SNORM;
        dataType = GL_SHORT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint32_t> >*>(volume)) {
        format = GL_RGBA;
        internalFormat = GL_RGBA;
        dataType = GL_UNSIGNED_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int32_t> >*>(volume)) {
        format = GL_RGBA;
        internalFormat = GL_RGBA32F;
        dataType = GL_INT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint64_t> >*>(volume)) {
        std::string message = "VolumeRAM_4xUInt64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int64_t> >*>(volume)) {
        std::string message = "VolumeRAM_4xInt64 not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<float> >*>(volume)) {
        format = GL_RGBA;
        internalFormat = GL_RGBA32F;
        dataType = GL_FLOAT;
    }
    else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<double> >*>(volume)) {
        std::string message = "VolumeRAM_4xDouble not supported as OpenGL volume.";
        LERROR(message);
        throw VoreenException(message);
    }
    else {
        LERROR("unknown or unsupported volume type");
    }

    //
    // Create texture
    //
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    VolumeTexture* vTex = 0;
    char* tempVolumeData = 0;
    if (volume->getDimensions().z > 1) { // multi-slice volume => just create OGL texture from it
        vTex = new VolumeTexture(static_cast<const GLubyte*>(volume->getData()),
            volume->getDimensions(),
            format, internalFormat, dataType, tgt::Texture::LINEAR);
    }
    else { // single-slice volume (not allowed as OGL texture) => double slice
        LWARNING("OpenGL does not allow 3D textures consisting of only one slice: cloning slice");
        try {
            tempVolumeData = new char[2*volume->getNumBytes()];
            memcpy(tempVolumeData, volume->getData(), volume->getNumBytes());
            memcpy(tempVolumeData+volume->getNumBytes(), volume->getData(), volume->getNumBytes());

            vTex = new VolumeTexture(reinterpret_cast<const GLubyte*>(tempVolumeData),
                tgt::svec3(volume->getDimensions().x, volume->getDimensions().y, volume->getDimensions().z * 2),
                format, internalFormat, dataType, tgt::Texture::LINEAR);
        }
        catch (std::bad_alloc&) {
            LERROR("bad allocation while creating OpenGL texture");
            throw VoreenException("bad allocation while creating OpenGL texture");
        }
    }
    tgtAssert(vTex, "volume texture not created");
    LGL_ERROR;

    vTex->bind();

    /// Make 1 and 2 channel texture appear to the shader as LUMINANCE and LUMINANCE_ALPHA for compatibility reasons
    if (format == GL_RED) {
        GLenum swizzleMask[4] = {GL_RED, GL_ZERO, GL_ZERO, GL_RED};
        vTex->setSwizzle((GLint*)swizzleMask);
    } else if (format == GL_RG) {
        //GLenum swizzleMask[4] = {GL_RED, GL_RED, GL_RED, GL_GREEN};
        GLenum swizzleMask[4] = {GL_RED, GL_GREEN, GL_ZERO, GL_GREEN}; // <- bugfix compared to the line above for two-channel volumes in SliceViewer
        vTex->setSwizzle((GLint*)swizzleMask);
    }

    if (volume->getData())
        vTex->uploadTexture();

    // set texture wrap to clamp
    vTex->setWrapping(tgt::Texture::CLAMP);

    LGL_ERROR;

    // prevent deleting twice
    vTex->setCpuTextureData(0,false);

    delete[] tempVolumeData;
    tempVolumeData = 0;

    // append to internal data structure
    texture_ = vTex;

    LGL_ERROR;
}

std::string VolumeGL::getFormat() const {
    return format_;
}

std::string VolumeGL::getBaseType() const {
    return baseType_;
}

size_t VolumeGL::getNumChannels() const {
    switch(getTexture()->getGLFormat()) {
        case GL_RED: return 1;
        case GL_RG: return 2;
        case GL_RGB: return 3;
        case GL_RGBA: return 4;
        default:
            LERROR("Unsupported volume format!");
            tgtAssert(false, "unsupported volume format");
            return 0;
    }
}

size_t VolumeGL::getBytesPerVoxel() const {
    int bytesPerChannel = 0;
    switch(getTexture()->getGLDataType()) {
        case GL_UNSIGNED_BYTE:
        case GL_BYTE:
            bytesPerChannel = 1;
            break;
        case GL_SHORT:
        case GL_UNSIGNED_SHORT:
            bytesPerChannel = 2;
            break;
        case GL_UNSIGNED_INT:
        case GL_INT:
            bytesPerChannel = 4;
            break;
        case GL_FLOAT:
            bytesPerChannel = 4;
            break;
        default:
            LERROR("Unsupported data type!");
    }
    return bytesPerChannel * getNumChannels();
}

//--------------------------------------------------------

bool RepresentationConverterUploadGL::canConvert(const VolumeRepresentation* source) const {
    // We can only convert from RAM volumes.
    if (!dynamic_cast<const VolumeRAM*>(source) )
        return false;

    // Test if a valid context is active.
    // TODO: ideally, we want to perform the conversion anyways.
    //       For that, we'd need it to be executed in the opengl/main thread the moment it's requested.
    if (tgt::GLContextManager::isInited() && !GLContextMgr.hasActiveContext())
        return false;

    return true;
}

VolumeRepresentation* RepresentationConverterUploadGL::convert(const VolumeRepresentation* source) const {
    const VolumeRAM* v = dynamic_cast<const VolumeRAM*>(source);

    if (!v) {
        //should have checked before...
        //LERROR("Failed to convert!");
        return 0;
    }

    VolumeGL* vgl = 0;
    try {
        vgl = new VolumeGL(v);
        LGL_ERROR;
    }
    catch (VoreenException& e) {
        LERRORC("voreen.RepresentationConverterUploadGL", "Failed to create VolumeGL: " << e.what());
    }
    catch (std::bad_alloc& /*e*/) {
        LERRORC("voreen.RepresentationConverterUploadGL", "Bad allocation during creation of VolumeGL");
    }
    return vgl;
}

//--------------------------------------------------------

bool RepresentationConverterDownloadGL::canConvert(const VolumeRepresentation* source) const {
    // We can only convert from RAM volumes.
    if (!dynamic_cast<const VolumeRAM*>(source) )
        return false;

    // Test if a valid context is active.
    // TODO: ideally, we want to perform the conversion anyways.
    //       For that, we'd need it to be executed in the opengl/main thread the moment it's requested.
    if (tgt::GLContextManager::isInited() && !GLContextMgr.hasActiveContext())
        return false;

    return true;
}

VolumeRepresentation* RepresentationConverterDownloadGL::convert(const VolumeRepresentation* source) const {
    const VolumeGL* vgl = dynamic_cast<const VolumeGL*>(source);

    if(vgl) {
        const VolumeTexture* tex = vgl->getTexture();
        tgt::svec3 dims = tex->getDimensions();
        VolumeRAM* v = 0;

        switch(tex->getGLFormat()) {
            case GL_RED:
                switch(tex->getGLDataType()) {
                    // VolumeRAM_UIntX
                    case GL_UNSIGNED_BYTE:
                        v = new VolumeRAM_UInt8(dims);
                        //internalFormat = GL_ALPHA8;
                        break;
                    case GL_UNSIGNED_SHORT:
                        v = new VolumeRAM_UInt16(dims);
                        //internalFormat = GL_ALPHA16;
                        break;
                    case GL_UNSIGNED_INT:
                        v = new VolumeRAM_UInt32(dims);
                        //internalFormat = GL_ALPHA;
                        break;
                        // VolumeRAM_IntX
                    case GL_BYTE:
                        v = new VolumeRAM_Int8(dims);
                        //internalFormat = GL_ALPHA8;
                        break;
                    case GL_SHORT:
                        v = new VolumeRAM_Int16(dims);
                        //internalFormat = GL_ALPHA16;
                        break;
                    case GL_INT:
                        v = new VolumeRAM_Int32(dims);
                        //internalFormat = GL_ALPHA;
                        break;
                        // VolumeRAM_Float and VolumeRAM_Double
                    case GL_FLOAT:
                        v = new VolumeRAM_Float(dims);
                        //internalFormat = GL_ALPHA;
                        break;
                    default:
                        //v = new VolumeRAM_Double(dims);
                        LERRORC("voreen.RepresentationConverterDownloadGL", "Unsupported texture format! (GL_ALPHA)");
                }
                break;
            case GL_RG:
                switch(tex->getGLDataType()) {
                    // VolumeRAM_2x with int16 types
                    case GL_UNSIGNED_SHORT:
                        v = new VolumeRAM_2xUInt16(dims);
                        //internalFormat = GL_LUMINANCE16_ALPHA16;
                        break;
                    case GL_SHORT:
                        v = new VolumeRAM_2xInt16(dims);
                        //internalFormat = GL_LUMINANCE16_ALPHA16;
                        break;
                    default:
                        LERRORC("voreen.RepresentationConverterDownloadGL", "Unsupported texture format! (GL_LUMINANCE_ALPHA)");
                        //v = new VolumeRAM_Double(dims);
                }
                break;
            case GL_RGB:
                switch(tex->getGLDataType()) {
                    // VolumeRAM_3x with int8 types
                    case GL_UNSIGNED_BYTE:
                        v = new VolumeRAM_3xUInt8(dims);
                        //internalFormat = GL_RGB8;
                        break;
                    case GL_BYTE:
                        v = new VolumeRAM_3xInt8(dims);
                        //internalFormat = GL_RGB8;
                        break;
                        // VolumeRAM_3x with int16 types
                    case GL_UNSIGNED_SHORT:
                        v = new VolumeRAM_3xUInt16(dims);
                        //internalFormat = GL_RGB16;
                        break;
                    case GL_SHORT:
                        v = new VolumeRAM_3xInt16(dims);
                        //internalFormat = GL_RGB16;
                        break;
                        // VolumeRAM_3x with real types
                    case GL_FLOAT:
                        v = new VolumeRAM_3xFloat(dims);
                        //internalFormat = GL_RGB;
                        break;
                    default:
                        LERRORC("voreen.RepresentationConverterDownloadGL", "Unsupported texture format! (GL_RGB)");
                        //v = new VolumeRAM_3xDouble(dims);
                }
                break;
            case GL_RGBA:
                switch(tex->getGLDataType()) {
                    // VolumeRAM_4x with int8 types
                    case GL_UNSIGNED_BYTE:
                        v = new VolumeRAM_4xUInt8(dims);
                        //internalFormat = GL_RGBA8;
                        break;
                    case GL_BYTE:
                        v = new VolumeRAM_4xInt8(dims);
                        //internalFormat = GL_RGBA8;
                        break;
                        // VolumeRAM_4x with int16 types
                    case GL_UNSIGNED_SHORT:
                        v = new VolumeRAM_4xUInt16(dims);
                        //internalFormat = GL_RGBA16;
                        break;
                    case GL_SHORT:
                        v = new VolumeRAM_4xInt16(dims);
                        //internalFormat = GL_RGBA16;
                        break;
                        // VolumeRAM_4x with real types
                    case GL_FLOAT:
                        v = new VolumeRAM_4xFloat(dims);
                        //internalFormat = GL_RGBA;
                        break;
                    default:
                        LERRORC("voreen.RepresentationConverterDownloadGL", "Unsupported texture format! (GL_RGBA)");
                        //v = new VolumeRAM_4xDouble(dims);
                }
                break;
            default:
                tgtAssert(false, "unsupported volume format");
        }

        if(v) {
            tex->downloadTextureToBuffer(static_cast<GLubyte*>(v->getData()), v->getNumBytes());
        }
        return v;
    }
    else {
        //should have checked before...
        //LERROR("Failed to convert!");
        return 0;
    }
}

} // namespace voreen
