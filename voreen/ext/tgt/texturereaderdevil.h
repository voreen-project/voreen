/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
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

#ifndef TGT_TEXTUREREADERDEVIL_H
#define TGT_TEXTUREREADERDEVIL_H

#ifdef TGT_HAS_DEVIL
#include "tgt/texturereader.h"

namespace tgt {

    /**
     * Texture reader using the DevIL library.
     * Supports the reading of many file formats.
     */
class TextureReaderDevil : public TextureReader {
public:
    /** Constructor. Sets readerName_ and supportedEndings_. */
    TextureReaderDevil();
    /** Destructor */
    virtual ~TextureReaderDevil();
    /** @override */
    virtual Texture* loadTexture(const std::string& filename, Texture::Filter filter, Texture::Wrapping wrapping, bool compress = false,
                                 bool keepPixels = false, bool uploadTexture = true, bool extend3To4Channels = false) override;
protected:
    static const std::string loggerCat_; ///< used for logging
};


} // namespace tgt

#endif // TGT_HAS_DEVIL
#endif // TGT_TEXTUREREADERDEVIL_H
