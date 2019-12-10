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

#ifndef TGT_FONT_H
#define TGT_FONT_H

#include "tgt/bounds.h"
#include <vector>
#include <string>
#include <map>
#include "tgt/tgt_gl.h"
#include "tgt/texture.h"
#include "stb/stb_truetype.h"
#include "fontmanager.h"


namespace tgt {

/**
 * The maximal font size where Text is rendered at native resolution
 * Bigger fontsizes are possible, but the text is upscaled from
 * this limit font size.
 */
const static float MAX_RENDER_FONT_SIZE = 60;

/**
 * This class represents a font to be rendered by OpenGL. A FTTextureFont will be created and
 * it is possible to render a text to a location or get the bounding box for the text beforehand.
 * By default a font size of 72 is set, but can be changed via the \sa setSize method.
 */
class TGT_API Font {
public:
    enum TextAlignment {
        TopLeft,        ///< Text is under and right to it's position
        TopCentered,    ///< Text is under and centered to it's position
        TopRight,       ///< Text is under and left to it's position

        MiddleLeft,     ///< Text is vertical centered and right to it's position
        Centered,       ///< Text is centered to it's position
        MiddleRight,    ///< Text is vertical centered and left to it's position

        BottomLeft,     ///< Text is above and right to it's position
        BottomCentered, ///< Text is above and centered to it's position
        BottomRight     ///< Text is above and left to it's position
    };

    /**
     * Creates a Font object from the given font file. A font size of 72 will be used by default
     * but that can be changed later on.
     * \param fontName The path to the font file, which should be used for this font object.
     * \param size The font size.
     */
    Font(const std::string& fontName
        , int size = 72
        , float lineWidth = 0.0f // this makes sure there are no unexpected line breaks
                                    // where exactly one line is expected
        , TextAlignment textAlignment = TopLeft);

    /**
     * Destructor - deletes the font object.
     */
    virtual ~Font();

    /**
     * Get the width of one line.
     */
    float getLineWidth();

    /**
     * Sets the line width to use.
     */
    void setLineWidth(float lineWidth);

    /**
     * Sets the font size which should be used from now on.
     * \param size The new font size
     */
    void setFontSize(int size);

    /**
     * Gets the font size.
     */
    int getFontSize();

    /**
     * Set font name.
     */
    void setFontName(const std::string& fontName);

    /**
     * Get font name.
     */
    std::string getFontName();

    /**
     * Get text alignment.
     */
    TextAlignment getTextAlignment();

    /**
     * Sets the font text alignment of the font.
     */
    void setTextAlignment(TextAlignment textAlignment);

    /**
     * Sets the color of the font
     */
    void setFontColor(vec4 color);

    /**
     * Get font color
     */
    vec4 getFontColor();

    /**
     * Gets the hight of a single line of text in pixel
     */
    float getLineHeight();

    /**
     * Sets horizontal and vertical padding in pixels of the font, so it always has this
     * distance to the position (except for middle/centered directions)
     */
    void setPadding(tgt::vec2 padding);

    /**
     * Sets horizontal and vertical padding in pixels of the font to the
     * same value, so it always has this
     * distance to the position (except for middle/centered directions)
     */
    void setPadding(float padding);

    /**
     * Gets the padding in pixels
     */
    tgt::vec2 getPadding() const;

    /**
     * Renders the text 'text' to the position 'pos'.
     * \sa pos The pen position of the first character
     * \sa text The text to be rendered
     * \sa screensize size of the current render target
     */
    void render(const tgt::vec3 pos, const std::string& text, const tgt::ivec2 screensize);


    /**
     * Renders the text to a new texture.
     * \sa text The text to be rendered
     */
    Texture* renderToTexture(const std::string& text);

    /**
     * Computes the bounding box for the the text 'text' beginning at position 'pos'.
     * \sa pos The pen position of the first character
     * \sa text The text for which the bounding box should be computed
     */
    tgt::vec2 getSize(const tgt::vec3& pos, const std::string& text, const tgt::ivec2 screensize);

private:
    struct CharacterVertex{
        vec2 pos;
        vec2 coord;
    };

    struct LineBuffer{
        LineBuffer() :width_(0){}
        float width_;
        std::vector<FontManager::codepoint_t> line_;
    };

    int           fontSize_;
    std::string   fontName_;
    TextAlignment alignment_;
    tgt::vec2     padding_;
    float         lineWidth_;
    tgt::vec4     color_;

    GLuint vao_;
    GLuint vbo_;

    float getEffektiveFontSize(float fontSize);
    float getScaleFactor(float fontSize);


    std::vector<LineBuffer> calculateLayoutGeometry(std::vector<FontManager::codepoint_t> text, float maxLineWidth, float sclae,
                                                    std::map<FontManager::codepoint_t,stbtt_packedchar>& chars, FontManager::FontInfo info);
    void addVerticesFromLineBuffer(float xoffset, float yoffset, float scale, LineBuffer line,
                                   std::map<FontManager::codepoint_t,stbtt_packedchar>& chars,
                                   std::vector<CharacterVertex> &verts, FontManager::FontInfo info);
    tgt::vec2 getSizeFromLayout(std::vector<LineBuffer>& layout);

};

} // namespace tgt

#endif
