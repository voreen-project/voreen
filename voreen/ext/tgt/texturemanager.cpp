/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
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

#include "tgt/texturemanager.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include "tgt/gpucapabilities.h"
#include "tgt/texturereader.h"
#include "tgt/filesystem.h"


namespace {

// cannot use tolower directly, due to some template overloading issues
int lower_case(int c) {
    return tolower(c);
}

} // namespace


namespace tgt {

//------------------------------------------------------------------------------
// TextureManager
//------------------------------------------------------------------------------

const std::string TextureManager::loggerCat_("tgt.TextureManager");

TextureManager::TextureManager()
{}

TextureManager::~TextureManager() {
    // clean up
    for (std::set<TextureReader*>::iterator iter = readerSet_.begin(); iter != readerSet_.end(); ++iter)
        delete *iter;
}

Texture* TextureManager::loadIgnorePath(const std::string& completeFilename, Texture::Filter filter, bool compress,
                              bool keepPixels, bool createOGLTex, bool useCache, bool extend3To4Channels)
{
    size_t found;
    std::string filename;

    found = completeFilename.find_last_of("/\\");

    if (found != std::string::npos){
        filename = completeFilename.substr(found+1);
    }
    else {
        filename = completeFilename;
    }

    if (compress && !GpuCaps.isTextureCompressionSupported())
        compress = false;

    if (useCache && isLoaded(filename)) {
        increaseUsage(filename);
        return get(filename);
    }

    Texture* t = 0;
    std::string ending = FileSys.fileExtension(filename);
    std::transform (ending.begin(), ending.end(), ending.begin(), lower_case);

    if (readers_.size() == 0)
        LWARNING("No TextureReaders are registered while calling TextureManager::loadIgnorePath()!");

    if (readers_.find(ending) != readers_.end()) {
        LDEBUG("Found matching reader: " << readers_[ending]->getReaderName());

        t = readers_[ending]->loadTexture(completeFilename, filter, Texture::REPEAT, compress,
                                            keepPixels, createOGLTex, extend3To4Channels);
        if (!t) {
            t = readers_[ending]->loadTexture(filename, filter, Texture::REPEAT, compress,
                                    keepPixels, createOGLTex, extend3To4Channels);
        }
        if ((t) && (useCache)) {
            reg(t, filename);
            return t;
        }
    }
    return t;
}

Texture* TextureManager::load(const std::string& filename, Texture::Filter filter, bool compress,
                              bool keepPixels, bool createOGLTex, bool useCache, bool extend3To4Channels)
{
    if (compress && !GpuCaps.isTextureCompressionSupported())
        compress = false;

    if (useCache && isLoaded(filename)) {
        increaseUsage(filename);
        return get(filename);
    }

    Texture* t = 0;
    std::string ending = FileSys.fileExtension(filename);
    std::transform (ending.begin(), ending.end(), ending.begin(), lower_case);

    if (readers_.size() == 0)
        LWARNING("No TextureReaders are registered while calling TextureManager::load()!");

    if (readers_.find(ending) != readers_.end()) {
        LDEBUG("Found matching reader: " << readers_[ending]->getReaderName());

        std::string completeFilename = completePath(filename);
        t = readers_[ending]->loadTexture(completeFilename, filter, Texture::REPEAT, compress,
                                          keepPixels, createOGLTex, extend3To4Channels);

        // else try just the filename without path
        if (!t && completeFilename != filename) {
            t = readers_[ending]->loadTexture(filename, filter, Texture::REPEAT, compress,
                                              keepPixels, createOGLTex, extend3To4Channels);
        }

        if (t && useCache)
            reg(t, filename);

        if (!t)
            LERROR("Texture file not found: " + filename);
    }
    else {
        LERROR("No matching reader found for ending '" << ending << "' while loading " << filename);
    }

    return t;
}

void TextureManager::registerReader(TextureReader* r) {
    readerSet_.insert(r);
    LDEBUG("TextureManager: Registering reader: " << r->getReaderName());

    std::string formats = "";
    std::vector<std::string> knownEndings = r->getSupportedEndings();
    std::vector<std::string>::iterator i;
    for ( i = knownEndings.begin(); i != knownEndings.end(); ++i ) {
        readers_[*i] = r;
        formats += *i + " ";
    }
    LDEBUG("Known formats: " << formats);
}

} // namespace tgt
