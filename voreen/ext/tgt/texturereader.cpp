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

#include "tgt/texturereader.h"

#include "tgt/gpucapabilities.h"
#include "tgt/logmanager.h"

namespace tgt {

const std::string TextureReader::loggerCat_("tgt.Texture.Reader");

TextureReader::TextureReader(){
}

TextureReader::~TextureReader(){
}

Texture* TextureReader::createTgtTexture(GLubyte* data, ivec3& dimensions, GLenum dataType, GLint format, Texture::Filter filter,
                                         Texture::Wrapping wrapping, bool compress, bool uploadTexture) {
    GLint internalFormat = format;
    // derived internal format from format and bit depth
    int bpt = Texture::calcFormatBpt(format,dataType);
    switch (format) {
    case GL_RED:
        if (bpt == 1)
            internalFormat = GL_R8;
        else if (bpt == 2)
            internalFormat = GL_R16;
        else {
            LERROR(bpt << " bytes per texel not supported for texture format GL_RED");
            return 0;
        }
        break;
    case GL_RG:
        if (bpt == 2)
            internalFormat = GL_RG;
        else if (bpt == 4)
            internalFormat = GL_RG16;
        else {
            LERROR(bpt << " bytes per texel not supported for texture format GL_RG");
            return 0;
        }
        break;
    case GL_RGB:
        if (bpt == 3)
            compress ? internalFormat = GL_COMPRESSED_RGB : internalFormat = GL_RGB;
        else if (bpt == 6)
            internalFormat = GL_RGB16;
        else {
            LERROR(bpt << " bytes per texel not supported for texture format GL_RGB");
            return 0;
        }
        break;
    case GL_RGBA:
        if (bpt == 4)
            compress ? internalFormat = GL_COMPRESSED_RGBA : internalFormat = GL_RGBA;
        else if (bpt == 8)
            internalFormat = GL_RGBA16;
        else {
            LERROR(bpt << " bytes per texel not supported for texture format GL_RGBA");
            return 0;
        }
        break;

    default:
        // unspecified or unknown format, try to derive from bpp
        LWARNING("Unspecified or unknown texture format, trying to derive format from bbp ...");
        switch (bpt) {
        case 1:
            internalFormat = GL_R8;
            break;
        case 2:
            internalFormat = GL_R16;
            break;
        case 3:
            compress ? internalFormat = GL_COMPRESSED_RGB : internalFormat = GL_RGB;
            break;
        case 4:
            compress ? internalFormat = GL_COMPRESSED_RGBA : internalFormat = GL_RGBA;
            break;
        case 8: // 16-bit-per-channel RGBA
            internalFormat = GL_RGBA16;
            break;
        case 12: //HDR-RGB, cut down to one byte per channel (until proper hdr-handling is implemented)
            internalFormat = GL_RGB;
            break;
        default:
            LERROR(bpt << " bytes per texel not supported!");
            return 0;
        }
    }

    if (!GpuCaps.isAnisotropicFilteringSupported() && filter == Texture::ANISOTROPIC)
        filter = Texture::MIPMAP;

    Texture* tex = new Texture(dimensions,format,internalFormat,dataType,filter,wrapping,data,true);

    if (uploadTexture) {
        tex->uploadTexture();
    }

    return tex;
}

} // namespace tgt
