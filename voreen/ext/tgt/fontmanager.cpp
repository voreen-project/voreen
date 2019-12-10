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

#include "stb/stb_truetype.h"
#include "tgt/fontmanager.h"
#include <limits>

#include "tgt/filesystem.h"


namespace tgt{

const std::string FontManager::loggerCat_("tgt.FontManager");
FontManager::FontManager()
{
    fontShaderProgram_ = 0;
    currentUseId_ = 0;
    fontcache_.reserve(MAX_FONT_CACHE_SIZE);
}

FontManager::~FontManager()
{
    if(fontShaderProgram_)
        ShdrMgr.dispose(fontShaderProgram_);
    for(size_t i = 0; i != fontcache_.size(); i++){
        fontcache_[i].deinitialize();
    }
}

FontManager::CacheEntry* FontManager::loadFontAndAddToCache( std::string name, int fontsize )
{
    size_t newIndex = 0;
    if (fontcache_.size()+1 >= MAX_FONT_CACHE_SIZE){
        uint64_t minId = std::numeric_limits<uint64_t>::max();
        for(size_t i = 0; i != fontcache_.size(); i++){
            if (fontcache_[i].lru_ < minId){
                minId = fontcache_[i].lru_;
                newIndex = i;
            }
        }
        fontcache_[newIndex].deinitialize();

    }else{
        newIndex = static_cast<int>(fontcache_.size());
        fontcache_.resize(newIndex+1);
    }

    CacheEntry &ce = fontcache_[newIndex];

    ce.initialize();

    ce.lru_ = currentUseId_;
    ce.fontname_ = name;
    ce.fontsize_ = fontsize;


    RegularFile file(name);
    file.seek(0, File::SeekDir::END);
    std::vector<unsigned char> fontBuffer;
    fontBuffer.resize(file.tell());
    file.seek(0, File::SeekDir::BEGIN);
    file.read(fontBuffer.data(), fontBuffer.size());

    stbtt_fontinfo fontinfo;
    stbtt_InitFont(&fontinfo, fontBuffer.data(), stbtt_GetFontOffsetForIndex(fontBuffer.data(), 0));
    int ascent, descent, lineGap;
    stbtt_GetFontVMetrics(&fontinfo, &ascent, &descent, &lineGap);
    ce.ascentPx_  = static_cast<float>(ascent)*stbtt_ScaleForPixelHeight(&fontinfo, fontsize);
    ce.descentPx_ = static_cast<float>(descent)*stbtt_ScaleForPixelHeight(&fontinfo, fontsize);

    return &fontcache_[newIndex];
}

Shader* FontManager::getFontShader()
{
    if (fontShaderProgram_ == 0){
        fontShaderProgram_ =  ShdrMgr.load("fontrendering", "", false);
    }
    return fontShaderProgram_;
}

/*GLuint FontManager::getFontTexture( std::string fontname, int fontsize )
{
    fontname = FileSys.cleanupPath(fontname, false);
    return getCacheEntry(fontname, fontsize)->texture_;
}

const std::map<FontManager::codepoint_t,stbtt_packedchar>& FontManager::getChars( std::string fontname, int fontSize )
{
    fontname = FileSys.cleanupPath(fontname, false);
    return getCacheEntry(fontname, fontSize)->chars_;
}

tgt::vec2 FontManager::getAscentDescent(std::string fontname, int fontSize)
{
    fontname = FileSys.cleanupPath(fontname, false);
    CacheEntry* ce = getCacheEntry(fontname, fontSize);
    return tgt::vec2(ce->ascentPx_, ce->descentPx_);
}*/

FontManager::FontInfo FontManager::getFontInfo(std::string fontname, int fontSize){
    fontname = FileSys.cleanupPath(fontname, false);
    CacheEntry* ce = getCacheEntry(fontname, fontSize);


    FontInfo info = FontInfo(ce->chars_);
    info.texture = ce->texture_;
    info.ascentDescent = tgt::vec2(ce->ascentPx_, ce->descentPx_);
    info.width = FONT_TEXTURE_WIDTH0;
    info.height = ce->texture_height;
    return info;
}

void FontManager::requestCharset(std::string fontname, int fontsize, std::vector<codepoint_t> codepoints)
{
    fontname = FileSys.cleanupPath(fontname, false);
    CacheEntry* ce =  getCacheEntry(fontname, fontsize);

    // find all different codepoints in the text
    std::sort(codepoints.begin(), codepoints.end());
    std::vector<codepoint_t>::iterator end = std::unique(codepoints.begin(), codepoints.end());

    // find codepoints that aren't loaded yet
    std::vector<int> codepointsRequested;
    for (std::vector<codepoint_t>::iterator begin = codepoints.begin(); begin != end; begin++){
        if (ce->chars_.count(*begin) == 0 && *begin != 10 && *begin != 13){
            // 10, 13 are CR and LF and we handle linebreaks in font.cpp

            // we need to add this char
            codepointsRequested.push_back(static_cast<int>(*begin));
        }
    }

    if (codepointsRequested.size() != 0){
        // load font from disk
        // (cache this?)
        RegularFile file(fontname);
        file.seek(0, File::SeekDir::END);
        std::vector<unsigned char> fontBuffer;
        fontBuffer.resize(file.tell());
        file.seek(0, File::SeekDir::BEGIN);
        file.read(fontBuffer.data(), fontBuffer.size());

        // rasterize missing codepoints
        std::vector<stbtt_packedchar> packedchars(codepointsRequested.size());
        stbtt_pack_range r;
        r.first_unicode_codepoint_in_range = 0;
        r.array_of_unicode_codepoints = codepointsRequested.data();
        r.font_size = static_cast<float>(fontsize);
        r.num_chars = static_cast<int>(codepointsRequested.size());
        r.chardata_for_range = packedchars.data();
        bool success = stbtt_PackFontRanges(&ce->packContext_, fontBuffer.data(), 0, &r, 1);
        if (!success){


            if (ce->growTexture()){
                LDEBUG("Texture growing");
                // add already loaded codepoints
                for(auto i = ce->chars_.begin(); i != ce->chars_.end(); i++){
                    codepoints.push_back(i->first);
                }

                ce->chars_.clear();

                requestCharset(fontname, fontsize, codepoints);

                return;
            }
            else{
                LERROR("Texture full. Some chars might not be shown.");
            }


        }


        for (size_t i = 0; i != codepointsRequested.size(); i++){
            ce->chars_.insert(std::pair<codepoint_t, stbtt_packedchar>(codepointsRequested[i], packedchars[i]));
        }

        // update texture on the gpu
        glBindTexture(GL_TEXTURE_2D, ce->texture_);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, FONT_TEXTURE_WIDTH0, ce->texture_height, 0, GL_RED, GL_UNSIGNED_BYTE, ce->textureCPU_.data());
        glGenerateMipmap(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

}

FontManager::CacheEntry* FontManager::getCacheEntry(std::string fontname, int fontsize)
{
    CacheEntry *c = nullptr;
    for(size_t i = 0; i != fontcache_.size(); i++){
        if (fontcache_[i].fontname_ == fontname
            && fontcache_[i].fontsize_ == fontsize){
            c = &fontcache_[i];
        }
    }
    if (!c){
        c = loadFontAndAddToCache(fontname, fontsize);
    }

    // update last used for lru cache
    currentUseId_++;
    c->lru_ = currentUseId_;
    return c;
}


void FontManager::CacheEntry::initialize()
{
    texture_height = FONT_TEXTURE_HEIGHT_START;
    allocateTexture();

}

void FontManager::CacheEntry::deinitialize()
{
    glDeleteTextures(1, &texture_);
    stbtt_PackEnd(&packContext_);
    chars_.clear();
}

void FontManager::CacheEntry::allocateTexture()
{
    textureCPU_.resize(FONT_TEXTURE_WIDTH0*texture_height);
    stbtt_PackBegin(&packContext_, textureCPU_.data(), FONT_TEXTURE_WIDTH0, texture_height, 0, 0, 0);
    stbtt_PackSetOversampling(&packContext_, 4, 4);
    chars_.clear();

    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, FONT_TEXTURE_WIDTH0, texture_height, 0, GL_RED, GL_UNSIGNED_BYTE, textureCPU_.data());
    glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
}

bool FontManager::CacheEntry::growTexture()
{

    if (texture_height*2 <= FONT_TEXTURE_HEIGHT_MAX){
        texture_height *= 2;
    }else{
        LERROR("Can't grow texture. Maximal size reached!");
        return false;
    }
    glDeleteTextures(1, &texture_);
    allocateTexture();
    return true;
}

}
