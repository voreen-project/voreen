/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#ifndef TGT_FONTMANAGER_H
#define TGT_FONTMANAGER_H

#include "tgt/bounds.h"
#include "tgt/tgt_gl.h"
#include "tgt/singleton.h"
#include "tgt/shadermanager.h"

#include <vector>
#include <string>
#include "stb/stb_truetype.h"

namespace tgt {
    class Font;
    /**
     * A FontManager manages the the shader and texture data for the Font class.
     * It includes a LRU Cache for font textures with MAX_FONT_CACHE_SIZE entrys.
     *
     * This class should not be used outside of the Font class.
     */
    class TGT_API FontManager : public Singleton<FontManager> {
    public:
        typedef uint32_t codepoint_t;

        const static int MAX_FONT_CACHE_SIZE = 15;
        //const static int FONT_SCALE = 10;
        const static int FONT_TEXTURE_WIDTH0 = 1024;
        const static int FONT_TEXTURE_HEIGHT_START = 128;
        const static int FONT_TEXTURE_HEIGHT_MAX = 4096;

        FontManager();
        virtual ~FontManager();


    private:
        friend Font;
        struct FontInfo{
            FontInfo(const std::map<codepoint_t,stbtt_packedchar>& chars_)
                :chars(chars_)
            {}

            GLuint texture;
            tgt::vec2 ascentDescent;
            int width;
            int height;
            const std::map<codepoint_t,stbtt_packedchar>& chars;
        };

        FontInfo getFontInfo(std::string fontname, int fontsize);

        /**
         * Get the texture for a font by name (may cache the texture)
         */
        //GLuint getFontTexture(std::string fontname, int fontsize);

        /**
         * Get the char data for a font by name (may cache the texture)
         */
        //const std::map<codepoint_t,stbtt_packedchar>& getChars(std::string fontname, int fontSize);

        tgt::vec2 getAscentDescent(std::string fontname, int fontsize);

        /**
         * Tell fontmanager we want to print all chars in the given text
         */
        void requestCharset(std::string fontname, int fontsize, std::vector<codepoint_t> codepoints);




        /**
         * Gets the shader for fontrendering.Get the texture for a font by name (may cache the texture)
         */
        tgt::Shader* getFontShader();



    private:
        /**
         * A Entry in the cache
         */
        struct CacheEntry{
            void initialize();

            void allocateTexture();

            bool growTexture();

            void deinitialize();

            std::string fontname_;
            int fontsize_;
            float ascentPx_;
            float descentPx_;
            GLuint texture_;
            int texture_height;
            uint64_t lru_;
            std::map<codepoint_t, stbtt_packedchar> chars_;
            stbtt_pack_context packContext_;
            std::vector<uint8_t> textureCPU_;
        };

        /**
         * An id for the lru cache, that increases with each use of a font
         */
        uint64_t currentUseId_;

        /**
         * the lru cache data
         */
        std::vector<CacheEntry> fontcache_;

        /**
         * The shader program for rendering fonts
         */
        tgt::Shader* fontShaderProgram_;


        /**
         * Loads a font and adds it to the cache.
         * This also may delete the least recently used entry from the
         * cache
         */
        CacheEntry* loadFontAndAddToCache(std::string name, int fontSize);

        /**
         * Get a CacheEntry and increses its lru value
         */
        CacheEntry* getCacheEntry(std::string name, int fontSize);




        static const std::string loggerCat_; ///< category used in logging
    };
#define FntMgr tgt::Singleton<tgt::FontManager>::getRef()


} // namespace tgt

#endif
