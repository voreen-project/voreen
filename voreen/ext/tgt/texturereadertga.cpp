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

#include "tgt/texturereadertga.h"

#include "tgt/logmanager.h"
#include "tgt/filesystem.h"

#include <cstring>

namespace tgt {

//------------------------------------------------------------------------------
// TextureReaderTga
//------------------------------------------------------------------------------

const std::string TextureReaderTga::loggerCat_("tgt.TextureReaderTga");

TextureReaderTga::TextureReaderTga()
    : TextureReader()
{
    readerName_ = "TGA Reader";
    supportedEndings_.push_back("tga");
}

Texture* TextureReaderTga::loadTexture(const std::string& filename, Texture::Filter filter, Texture::Wrapping wrapping, bool compress,
                                       bool keepPixels, bool uploadTexture, bool extend3To4Channels)
{

    GLubyte TGAheader[12];
    GLubyte header[6];

    File* file = FileSys.open(filename);

    // Check if file is open
    if (!file) {
        LERROR("Failed to open file " << filename);
        return 0;
    }

    if (!file->isOpen()) {
        LERROR("Failed to open file " << filename);
        delete file;
        return 0;
    }

    size_t len = file->size();

    // check if file is empty
    if (len == 0) {
        delete file;
        return 0;
    }

    if (file->read(reinterpret_cast<char*>(&TGAheader), sizeof(TGAheader)) != sizeof(TGAheader) ||
        file->read(reinterpret_cast<char*>(&header), sizeof(header)) != sizeof(header))
    {
        delete file;
        LERROR("Failed to read header! file: " << filename);
        return 0;
    }

    ivec3 dimensions;
    dimensions.x = header[1] * 256 + header[0];           // determine the TGA width  (highbyte*256+lowbyte)
    dimensions.y = header[3] * 256 + header[2];           // determine the TGA height (highbyte*256+lowbyte)
    dimensions.z = 1;
    LDEBUG("Image dimensions: " << dimensions);

    if (dimensions.x <= 0 || dimensions.y <= 0) {
        delete file;
        LERROR("wrong dimensions: " << dimensions << " file: " << filename);
        return 0;
    }

    if (header[4] != 24 && header[4] != 32) {
        delete file;
        LERROR("Illegal bytes per pixel! file: " << filename);
        return 0;
    }

    int bpp = header[4];
    int bpt = bpp / 8;  // divide by 8 to get the bytes per texel

    GLint format;
    switch (bpt) {
    case 3:
        format = GL_RGB;
        LDEBUG("RGB");
        break;
    case 4:
        format = GL_RGBA;
        LDEBUG("RGBA");
        break;
    default:
        delete file;
        LERROR("unsupported byptes per pixel " << filename);
        return 0;
    }

    //load file into buffer
    size_t bufferSize = hmul(dimensions)*bpt;
    unsigned char* cpuData = new unsigned char[bufferSize];
    unsigned char* currentIndex = cpuData;

    if (TGAheader[2] == 2) {
        // file is not compressed
        LDEBUG("Reading uncompressed TGA file...");
        if (file->read(cpuData, bufferSize) != bufferSize) {
            LERROR("Failed to read uncompressed image! file: " << filename);
            delete file;
            delete[] cpuData;
            return 0;
        }
    } else {
        // file is compressed
        LDEBUG("Reading compressed TGA file " << filename << " ...");

        //TODO: error handling
        unsigned char chunk[4];
        for (unsigned int bytesDone=0; bytesDone < bufferSize; ) {
            unsigned char packetHead;
            file->read(reinterpret_cast<char*>(&packetHead), 1);
            if (packetHead > 128) {
                //RLE
                packetHead -= 127;
                file->read(reinterpret_cast<char*>(&chunk), bpt);
                for (unsigned char b=0; b < packetHead; b++) {
                        std::memcpy(currentIndex, chunk, bpt);
                        bytesDone += bpt;
                        currentIndex += bpt;
                }
            } else {
                //RAW
                packetHead++;
                file->read(currentIndex, bpt * packetHead);
                bytesDone += packetHead * bpt;
                currentIndex += packetHead * bpt;
            }
        }
    }

    file->close();
    delete file;

    //create new texture
    Texture* tex = createTgtTexture(reinterpret_cast<GLubyte*>(cpuData),dimensions,GL_UNSIGNED_BYTE,format,filter,wrapping,compress,false);
    if (!tex) {
        delete[] cpuData;
        LERROR("Failed to create texture for " << filename);
        return 0;
    }
    tex->setOptionalName(filename);

    // switch r & b
    uint8_t temp;
    switch(bpt) {
    case 3:
        for (int i=0; i < (dimensions.x*dimensions.y); ++i) {
            temp = tex->texel<col3>(i).r;
            tex->texel<col3>(i).r = tex->texel<col3>(i).b;
            tex->texel<col3>(i).b = temp;
        }
        break;
    case 4:
        for (int i=0; i < (dimensions.x*dimensions.y); ++i) {
            temp = tex->texel<col4>(i).r;
            tex->texel<col4>(i).r = tex->texel<col4>(i).b;
            tex->texel<col4>(i).b = temp;
        }
        break;
    }

    // flip image horizontally
    //col3 tmp;
    //for (int y=0; y < dimensions.y/2; ++y) {
        //for (int x=0; x < dimensions.x; ++x) {
            //tmp = t->texel<col3>(x,y);
            //t->texel<col3>(x,y) = t->texel<col3>(x,dimensions.y-1-y);
            //t->texel<col3>(x,dimensions.y-1-y) = tmp;
        //}
    //}

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
