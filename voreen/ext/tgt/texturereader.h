/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
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

#ifndef TGT_TEXTUREREADER_H
#define TGT_TEXTUREREADER_H

#include "tgt/texture.h"

#include <vector>

namespace tgt {

    /**
     * This is the base class for a simple texture reading from differnet file formats.
     * @see texturereadertga etc...
     * @note For deriving a new reader, the loadTexture() must be overriden and the readerName_ and
     *       supportedEndings_ must be set in teh constructor.
     */
class TextureReader {
public:
    /** Constructor */
    TextureReader();
    /** Destructor */
    virtual ~TextureReader();

    /** Returns all supported file endings (formats) */
    virtual const std::vector<std::string>& getSupportedEndings() const { return supportedEndings_; }
    /** Returns the name of the reader. For debug purpuse */
    virtual std::string getReaderName() const { return readerName_; }

    /**
     * Main function of each reader. It should load and create the texture. The base
     * implementation is pure virtual and must be overriden.
     * By default, the new texture should be bind automatically.
     */
    virtual Texture* loadTexture(const std::string& filename, Texture::Filter filter, Texture::Wrapping wrapper,
                                 bool compress = false, bool keepPixels = false,
                                 bool uploadTextureAfterCreate = true, bool extend3To4Channels = false) = 0;

protected:
    //------------------
    //  Helper
    //------------------
    /** Called from loadTexture */
    Texture* createTgtTexture(GLubyte* data, ivec3& dimensions, GLenum dataType, GLint format, Texture::Filter filter,
                              Texture::Wrapping wrapping, bool compress, bool uploadTexture = true);

    //------------------
    //  Members
    //------------------
    std::vector<std::string> supportedEndings_; ///< Endings, which are supported by the reader
    std::string readerName_;                    ///< Name of the reader. Used like a class name.

    static const std::string loggerCat_;        ///< Used for logging
};

} // namespace tgt

#endif // TGT_TEXTUREREADER_H
