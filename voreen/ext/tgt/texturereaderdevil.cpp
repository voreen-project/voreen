/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#include "tgt/texturereaderdevil.h"

#ifdef TGT_HAS_DEVIL

#include <IL/il.h>

#include "tgt/logmanager.h"
#include "tgt/filesystem.h"
#include <cstring>

namespace tgt {

//------------------------------------------------------------------------------
// TextureReaderDevil
//------------------------------------------------------------------------------


const std::string TextureReaderDevil::loggerCat_("tgt.TextureReaderDevIl");

TextureReaderDevil::TextureReaderDevil()
    : TextureReader()
{
    readerName_ = "DevIL Reader";
    supportedEndings_.push_back("bmp");
    supportedEndings_.push_back("cut");
    supportedEndings_.push_back("dcx");
    supportedEndings_.push_back("dds");
    supportedEndings_.push_back("ico");
    supportedEndings_.push_back("gif");
    supportedEndings_.push_back("jpg");
    supportedEndings_.push_back("jpeg");
    supportedEndings_.push_back("lbm");
    supportedEndings_.push_back("lif");
    supportedEndings_.push_back("mdl");
    supportedEndings_.push_back("pcd");
    supportedEndings_.push_back("pcx");
    supportedEndings_.push_back("pic");
    supportedEndings_.push_back("png");
    supportedEndings_.push_back("pnm");
    supportedEndings_.push_back("psd");
    supportedEndings_.push_back("psp");
    supportedEndings_.push_back("raw");
    supportedEndings_.push_back("sgi");
    supportedEndings_.push_back("tga");
    supportedEndings_.push_back("tif");
    supportedEndings_.push_back("wal");
    supportedEndings_.push_back("act");
    supportedEndings_.push_back("pal");
    supportedEndings_.push_back("hdr");

    // Initialize DevIL
    ilInit();
    ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
    ilEnable(IL_ORIGIN_SET); // Flip images
}

TextureReaderDevil::~TextureReaderDevil() {
}

Texture* TextureReaderDevil::loadTexture(const std::string& filename, Texture::Filter filter, Texture::Wrapping wrapping, bool compress,
                                         bool keepPixels, bool uploadTexture, bool extend3To4Channels)
{
    File* file = FileSys.open(filename);

    // check if file is open
    if (!file || !file->isOpen()) {
        delete file;
        return 0;
    }

    size_t len = file->size();

    // check if file is empty
    if (len == 0) {
        delete file;
        return 0;
    }

    // allocate memory
    char* imdata = new char[len];

    if (imdata == 0) {
        delete file;
        return 0;   // allocation failed
    }

    file->read(imdata, len);

    file->close();
    delete file;

    //create devil image
    ILuint ImageName;
    ilGenImages(1, &ImageName);
    ilBindImage(ImageName);

    if (!ilLoadL(IL_TYPE_UNKNOWN, imdata, static_cast<ILuint>(len))) {
        LERROR("Failed to open via ilLoadL " << filename);
        delete[] imdata;
        return 0;
    }
    delete[] imdata;
    imdata = 0;

    //parameters needed for the texture
    GLint format; GLenum dataType;
    size_t bufferSize;
    bool swizzleRG = false; // HACK: map RG texture to illuminace-alpha

    // determine image format
    ILint devilFormat;
    switch (ilGetInteger(IL_IMAGE_FORMAT)) {
        case IL_LUMINANCE:         // intensity channel only
            devilFormat = IL_LUMINANCE;
            format = GL_RED;
            break;
        case IL_LUMINANCE_ALPHA:   // intensity-alpha channels
            devilFormat = IL_LUMINANCE_ALPHA;
            format = GL_RG;
            swizzleRG = true;
            break;
        case IL_RGB:
        case IL_BGR:
            if(extend3To4Channels) {
                devilFormat = IL_RGBA;  // four color channels
                format = GL_RGBA;
            } else {
                devilFormat = IL_RGB;  // three color channels
                format = GL_RGB;
            }
            break;
        case IL_RGBA:
        case IL_BGRA:
            devilFormat = IL_RGBA; // R-G-B-A ordered color channels, convert to RGBA if necessary
            format = GL_RGBA;
            break;
        default:
            LERROR("unsupported format: " << ilGetInteger(IL_IMAGE_FORMAT) << " (" << filename << ")");
            ilDeleteImages(1, &ImageName);
            return 0;
    }

    // determine data type
    ILint devilDataType;
    switch (ilGetInteger(IL_IMAGE_TYPE)) {
    case IL_UNSIGNED_BYTE:
        devilDataType = IL_UNSIGNED_BYTE;
        dataType = GL_UNSIGNED_BYTE;
        break;
    case IL_BYTE:
        devilDataType = IL_BYTE;
        dataType = GL_BYTE;
        break;
    case IL_UNSIGNED_SHORT:
        devilDataType = IL_UNSIGNED_SHORT;
        dataType = GL_UNSIGNED_SHORT;
        break;
    case IL_SHORT:
        devilDataType = IL_SHORT;
        dataType = GL_SHORT;
        break;
    case IL_UNSIGNED_INT:
        devilDataType = IL_UNSIGNED_INT;
        dataType = GL_UNSIGNED_INT;
        break;
    case IL_INT:
        devilDataType = IL_INT;
        dataType = GL_INT;
        break;
    case IL_FLOAT:
        devilDataType = IL_FLOAT;
        dataType = GL_FLOAT;
        break;
    default:
        LERROR("unsupported data type: " << ilGetInteger(IL_IMAGE_TYPE) << " (" << filename << ")");
        ilDeleteImages(1, &ImageName);
        return 0;
    }

    //converts the bgr(a) textures to rgb(a)
    if (!ilConvertImage(devilFormat, devilDataType)) {
        LERROR("failed to convert loaded image: " << filename);
        ilDeleteImages(1, &ImageName);
        return 0;
    }

    int bpt = ilGetInteger(IL_IMAGE_BYTES_PER_PIXEL);

    tgt::ivec3 dims;
    dims.x = ilGetInteger(IL_IMAGE_WIDTH);
    dims.y = ilGetInteger(IL_IMAGE_HEIGHT);
    dims.z = ilGetInteger(IL_IMAGE_DEPTH);
    LDEBUG("Image dimensions: " << dims);
    tgtAssert( dims.z == 1, "depth is not equal 1");


    bufferSize = hmul(dims)*bpt;

    GLubyte* cpuData = new GLubyte[bufferSize];

    memcpy(cpuData, ilGetData(), bufferSize);

    ilDeleteImages(1, &ImageName);

    //create new texture
    Texture* tex = createTgtTexture(cpuData,dims,dataType,format,filter,wrapping,compress,false);
    if (!tex) {
        delete[] cpuData;
        LERROR("Failed to create texture for " << filename);
        return 0;
    }
    tex->setOptionalName(filename);
    if(swizzleRG) { //HACK: map GL_RGto Illuminance-Alpha
        GLint swizzleMask[4] = {GL_RED,GL_RED,GL_RED,GL_GREEN};
        tex->setSwizzle(swizzleMask);
    }

    //upload texture if needed
    if(uploadTexture)
        tex->uploadTexture();

    //delete cpu version
    if (!keepPixels) {
        tex->setCpuTextureData(0,false);
    }

    //return
    return tex;
}

} // namespace tgt

#endif // TGT_HAS_DEVIL
