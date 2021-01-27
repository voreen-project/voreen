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

#define STB_RECT_PACK_IMPLEMENTATION
#include "stb/stb_rect_pack.h"
#undef STB_RECT_PACK_IMPLEMENTATION

#define STB_TRUETYPE_IMPLEMENTATION  // force following include to generate implementation
#include "stb/stb_truetype.h"
#undef STB_TRUETYPE_IMPLEMENTATION

#include "tgt/font.h"

#include "tgt/logmanager.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/tgt_gl.h"
#include "tgt/shadermanager.h"
#include "tgt/matrixstack.h"
#include "tgt/fontmanager.h"
#include "tgt/framebufferobject.h"

#include <algorithm>
#include <sstream>

namespace tgt {
static std::vector<FontManager::codepoint_t> decodeUTF8(const unsigned char* string){

const unsigned char* input = string;
bool err = false; //< set when the text is not correct utf8
std::vector<FontManager::codepoint_t> ucs32;
while (*input && !err){
    FontManager::codepoint_t codepoint = 0;
    unsigned char first = *input++;
    if (first > 0xF0){ // read 3 bytes
        codepoint = first & 0x03;
        if ((*input & 0xC0) != 0x80){ err = true; break; }
        codepoint = (codepoint << 6) | ((*input++) & 0x3F);
        if ((*input & 0xC0) != 0x80){ err = true; break; }
        codepoint = (codepoint << 6) | ((*input++) & 0x3F);
        if ((*input & 0xC0) != 0x80){ err = true; break; }
        codepoint = (codepoint << 6) | ((*input++) & 0x3F);
    }
    else if (first > 0xD0){ // read 2 bytes
        codepoint = first & 0x07;
        if ((*input & 0xC0) != 0x80){ err = true; break; }
        codepoint = (codepoint << 6) | ((*input++) & 0x3F);
        if ((*input & 0xC0) != 0x80){ err = true; break; }
        codepoint = (codepoint << 6) | ((*input++) & 0x3F);
    }
    else if (first > 0xB0){
        codepoint = first & 0x0F;
        if ((*input & 0xC0) != 0x80){ err = true; break; }
        codepoint = (codepoint << 6)|((*input++) & 0x3F);
    }
    else if (first < 0x80){
        codepoint = first;
    }
    else{
        err = true;
        break;
    }
    ucs32.push_back(codepoint);
}
if (!err){
    return ucs32;
}
else{
    // fall back to latin1 encoding
    ucs32.clear();

    for (input = string; *input; input++)
        ucs32.push_back(*input);
    return ucs32;
}
}

static const float POINT_TO_PX = 1.15; // choosen to be close to old font rendering

Font::Font(const std::string& fontName, int size, float lineWidth, TextAlignment textAlignment)
    : fontName_(fontName)
    , fontSize_(size)
    , lineWidth_(lineWidth)
    , alignment_(textAlignment)
    , padding_(tgt::vec2::zero)
    , vao_(0)
    , vbo_(0)
    , color_(vec4(1.0f))
{
}

Font::~Font() {
    if (vao_)
        glDeleteVertexArrays(1, &vao_);
    if (vbo_)
        glDeleteBuffers(1, &vbo_);
}



float Font::getLineWidth() {
    return lineWidth_;
}

void Font::setLineWidth(float lineWidth) {
    lineWidth_ = lineWidth;
}

void Font::setFontSize(int size) {
    fontSize_ = size;
}

int Font::getFontSize() {
    return fontSize_;
}

void Font::setFontName(const std::string& fontName) {
    fontName_ = fontName;
}

std::string Font::getFontName() {
    return fontName_;
}

void Font::setTextAlignment(TextAlignment textAlignment) {
    alignment_ = textAlignment;
}
Font::TextAlignment Font::getTextAlignment(){
    return alignment_;
}

void Font::setFontColor(vec4 color){
    color_ = color;
}

vec4 Font::getFontColor(){
    return color_;
}

void Font::setPadding(tgt::vec2 padding){
    padding_ = padding;
}
void Font::setPadding(float padding){
    setPadding(tgt::vec2(padding));
}
tgt::vec2 Font::getPadding() const{
    return padding_;
}

float Font::getLineHeight(){
    return (float)static_cast<int>(getEffektiveFontSize(1.0f*fontSize_)*getScaleFactor(1.0f*fontSize_)*POINT_TO_PX);
}

void Font::render(const tgt::vec3 pos, const std::string& text, const tgt::ivec2 screensize) {
    if (text.empty())
        return;

    if (screensize.x <= 0 || screensize.y <= 0)
        return;

    // store blend status to reset it afterwards
    tgt::vec3 realpos = tgt::round(pos);
    GLboolean blend;
    GLint blend_src_rgb;
    GLint blend_src_alpha;
    GLint blend_dst_rgb;
    GLint blend_dst_alpha;
    glGetBooleanv(GL_BLEND, &blend);
    glGetIntegerv(GL_BLEND_DST_ALPHA, &blend_dst_alpha);
    glGetIntegerv(GL_BLEND_DST_RGB, &blend_dst_rgb);
    glGetIntegerv(GL_BLEND_SRC_ALPHA, &blend_src_alpha);
    glGetIntegerv(GL_BLEND_SRC_RGB, &blend_src_rgb);


    const bool left  = alignment_ == TopLeft || alignment_ == MiddleLeft || alignment_ == BottomLeft;
    const bool right = alignment_ == TopRight || alignment_ == MiddleRight || alignment_ == BottomRight;
    const bool center = !left && !right;

    const bool top = alignment_ == TopLeft || alignment_ == TopCentered || alignment_ == TopRight;
    const bool bottom = alignment_ == BottomLeft || alignment_ == BottomCentered || alignment_ == BottomRight;
    const bool middle = !top && ! bottom;

    float lineWidth = lineWidth_;
    if (lineWidth == 0){
        if (left){
            lineWidth = screensize.x-realpos.x;
        }else if (right){
            lineWidth = realpos.x;
        }else if(center){
            lineWidth = std::min(screensize.x-realpos.x, realpos.x);
        }
    }

    lineWidth -= padding_.x;

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    tgt::Shader* shader = FntMgr.getFontShader();

    float scaleFactor = getScaleFactor(1.0f*fontSize_);
    float fontSizePx = getEffektiveFontSize(1.0f*fontSize_)*POINT_TO_PX;

    std::vector<FontManager::codepoint_t> textUCS32 = decodeUTF8(reinterpret_cast<const unsigned char*>(text.c_str()));
    FntMgr.requestCharset(fontName_, static_cast<int>(fontSizePx), textUCS32);

    FontManager::FontInfo info = FntMgr.getFontInfo(fontName_, static_cast<int>(fontSizePx));
    info.ascentDescent *= scaleFactor;

    std::map<FontManager::codepoint_t, stbtt_packedchar>& chars = const_cast<std::map<FontManager::codepoint_t, stbtt_packedchar>&>(info.chars);
    std::vector<CharacterVertex> characterVerticies;

    std::vector<LineBuffer> layout = calculateLayoutGeometry(textUCS32, lineWidth, scaleFactor, chars, info);
    float offsety = 0;

    if (top){
        offsety = -getLineHeight()*(layout.size())-padding_.y+info.ascentDescent.x;
    }else if(middle){
        offsety = +(info.ascentDescent.x)/2-getLineHeight()*(layout.size()-1)/2.0f;
        //offsety = 0;
    }else if(bottom){
        offsety = +info.ascentDescent.x+padding_.y;
    }

    //vec2 size = getSizeFromLayout(layout)*scaleFactor;

    for(size_t i = 0; i != layout.size(); i++){
        float offsetx = 0;
        if (left){
            offsetx = padding_.x;
        }else if(center){
            offsetx = (lineWidth-layout[i].width_)/2.0f;
        }else if(right){
            offsetx = -layout[i].width_-padding_.x;
        }
        addVerticesFromLineBuffer(offsetx, (offsety+getLineHeight()*i), scaleFactor,
                                    layout[i], chars, characterVerticies, info);
    }

    // render
    mat4 M = mat4::createScale(vec3(1.0f, -1.0f , 1.0));
    M = mat4::createTranslation(realpos)*M;
    M = mat4::createOrtho(0, (float)screensize.x, 0, (float)screensize.y, -1, 1)*M;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, info.texture);

    shader->activate();
    shader->setUniform("MVP", M);
    shader->setUniform("sampler_", 0);
    shader->setUniform("color", color_);

    if (!vbo_){
        glGenBuffers(1, &vbo_);
    }
    if (!vao_){
        glGenVertexArrays(1, &vao_);
    }

    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(CharacterVertex)*characterVerticies.size(),
                    characterVerticies.data(), GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(CharacterVertex), (GLvoid*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(CharacterVertex),
                            (GLvoid*)offsetof(CharacterVertex, coord));
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)characterVerticies.size());

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    shader->deactivate();
    glBindTexture(GL_TEXTURE_2D, 0);

    // reset opengl state to previous values
    if(!blend){
        glDisable(GL_BLEND);
    }
    glBlendFuncSeparate(blend_src_rgb, blend_dst_rgb, blend_src_alpha, blend_dst_alpha);
}

tgt::vec2 Font::getSize(const tgt::vec3& pos, const std::string& text, const tgt::ivec2 screensize) {
    float lineWidth = lineWidth_;
    if (lineWidth == 0){
        lineWidth = screensize.x-pos.x;
    }

    std::vector<FontManager::codepoint_t> textUCS32 = decodeUTF8(reinterpret_cast<const unsigned char*>(text.c_str()));
    FntMgr.requestCharset(fontName_, static_cast<int>(getEffektiveFontSize(1.0f*fontSize_)*POINT_TO_PX), textUCS32);
    FontManager::FontInfo info = FntMgr.getFontInfo(fontName_, static_cast<int>(getEffektiveFontSize(1.0f*fontSize_)*POINT_TO_PX));
    std::map<FontManager::codepoint_t, stbtt_packedchar>& chars = const_cast<std::map<FontManager::codepoint_t, stbtt_packedchar>&>(info.chars);
    std::vector<LineBuffer> layout = calculateLayoutGeometry(textUCS32, lineWidth, getScaleFactor(1.0f*fontSize_), chars, info);
    return getSizeFromLayout(layout);
}

std::vector<Font::LineBuffer> Font::calculateLayoutGeometry(std::vector<FontManager::codepoint_t> text, float maxLineWidth, float scale, std::map<FontManager::codepoint_t,stbtt_packedchar>& chars, FontManager::FontInfo info){
    std::vector<LineBuffer> layout;
    LineBuffer line;
    float x = 0;
    float y = 0;
    for(size_t i = 0; i != text.size(); i++){
        FontManager::codepoint_t c = text[i];

        if (c == '\n'){
            layout.push_back(line);
            line = LineBuffer();
            x = 0;
            y = 0;
            continue;
        }

        if (chars.count(c) == 0){
            // should not happen
            tgtAssert(false, "Unexpected char.")
            continue;
        }

        stbtt_aligned_quad q;
        stbtt_GetPackedQuad(&chars[c],
            info.width, info.height,
            0, &x, &y, &q, 1);
        if (x*scale < maxLineWidth){
            line.width_ = x*scale;
            line.line_.push_back(c);
        }else{
            layout.push_back(line);
            line = LineBuffer();
            x = 0;
            y = 0;
            stbtt_GetPackedQuad(&chars[c], info.width, info.height,
                                0, &x,&y,&q,1);
            line.width_ = x*scale;
            line.line_.push_back(c);
        }
    }
    layout.push_back(line);
    return layout;
}


void Font::addVerticesFromLineBuffer(float xoffset, float yoffset, float scale, LineBuffer line, std::map<FontManager::codepoint_t,stbtt_packedchar>& chars, std::vector<CharacterVertex> &verts, FontManager::FontInfo info){
    float x = 0;
    float y = 0;

    for(size_t i = 0; i != line.line_.size(); i++){
        FontManager::codepoint_t c = line.line_[i];

        //if (c < 32 || c >= 128) continue; // TODO: necessary?

        stbtt_aligned_quad q;
        stbtt_GetPackedQuad(&chars[c], info.width, info.height, 0, &x, &y, &q, 1);

        CharacterVertex vertex;

        // First Triangle
        vertex.pos.x = q.x0;
        vertex.pos.y = q.y1;
        vertex.coord.x = q.s0;
        vertex.coord.y = q.t1;
        vertex.pos = scale*vertex.pos+vec2(xoffset, yoffset);
        verts.push_back(vertex);

        vertex.pos.x = q.x1;
        vertex.pos.y = q.y0;
        vertex.coord.x = q.s1;
        vertex.coord.y = q.t0;
        vertex.pos = scale*vertex.pos+vec2(xoffset, yoffset);
        verts.push_back(vertex);

        vertex.pos.x = q.x0;
        vertex.pos.y = q.y0;
        vertex.coord.x = q.s0;
        vertex.coord.y = q.t0;
        vertex.pos = scale*vertex.pos+vec2(xoffset, yoffset);
        verts.push_back(vertex);

        // Second Triangle
        vertex.pos.x = q.x0;
        vertex.pos.y = q.y1;
        vertex.coord.x = q.s0;
        vertex.coord.y = q.t1;
        vertex.pos = scale*vertex.pos+vec2(xoffset, yoffset);
        verts.push_back(vertex);

        vertex.pos.x = q.x1;
        vertex.pos.y = q.y0;
        vertex.coord.x = q.s1;
        vertex.coord.y = q.t0;
        vertex.pos = scale*vertex.pos+vec2(xoffset, yoffset);
        verts.push_back(vertex);

        vertex.pos.x = q.x1;
        vertex.pos.y = q.y1;
        vertex.coord.x = q.s1;
        vertex.coord.y = q.t1;
        vertex.pos = scale*vertex.pos+vec2(xoffset, yoffset);
        verts.push_back(vertex);
    }
}

tgt::vec2 Font::getSizeFromLayout(std::vector<LineBuffer>& layout){
    tgt::vec2 size(0.0f);
    size.y = static_cast<float>(layout.size())*getLineHeight();
    for(size_t i = 0; i != layout.size(); i++){
        size.x = std::max(size.x, layout[i].width_);
    }
    return size;
}

Texture* Font::renderToTexture(const std::string& text)
{
    // Restore state afterwards.
    GLFramebufferObjectGuard guard;
    GLint oldViewport[4];
    glGetIntegerv( GL_VIEWPORT, oldViewport );
    float oldLinewidth = getLineWidth();
    
    TextAlignment alignment = getTextAlignment();
    //tgt::vec4 color = getFontColor();
    setTextAlignment(TopLeft);
    //setFontColor(tgt::vec4::one);
    setLineWidth(0);

    tgt::ivec2 size = tgt::ceil(getSize(tgt::vec3::zero, text, tgt::ivec2(4096, 4096)))+tgt::vec2(1, 1);
    Texture *tex = new Texture(tgt::ivec3(size, 1), GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE);

    // Create fbo.
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glBindTexture(GL_TEXTURE_2D, tex->getId());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA, tex->getDimensions().x, tex->getDimensions().y);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size.x, size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex->getId(), 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    //-------------------------
    //Attach depth buffer to FBO.
    GLuint depth_rb;
    glGenRenderbuffers(1, &depth_rb);
    glBindRenderbuffer(GL_RENDERBUFFER, depth_rb);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, size.x, size.y);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_rb);
    GLuint err = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    GLenum buffers[1] = {GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, buffers);
    tgtAssert(FramebufferObject::isComplete(), "tgt::Font::renderToTexture: Framebuffer not complete!");
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0, 0, size.x, size.y);
    render(tgt::vec3::zero, text, size);
    glDeleteFramebuffers(1, &fbo);
    glDeleteRenderbuffers(1, &depth_rb);

    // Restore state.
    glViewport(oldViewport[0], oldViewport[1], oldViewport[2], oldViewport[3]);
    setTextAlignment(alignment);
    //setFontColor(color);
    setLineWidth(oldLinewidth);
    return tex;
}

float Font::getEffektiveFontSize(float fontSize)
{
    if (fontSize > MAX_RENDER_FONT_SIZE)
        return MAX_RENDER_FONT_SIZE;
    else
        return fontSize;
}

float Font::getScaleFactor(float fontSize)
{
    return fontSize/getEffektiveFontSize(fontSize);
}

} // namespace tgt
