/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "touchtablepreviewicon.h"


namespace voreen {
    TouchTablePreviewIcon::TouchTablePreviewIcon(tgt::ivec2 posLL, int width, int height, tgt::Texture* symbolTex, std::string filename, std::string caption, bool isSelected)
        : bounds_(tgt::vec3(posLL, 1.f), tgt::vec3(posLL.x + (float) width, posLL.y + (float) height, 1.f))
        , width_(width)
        , height_(height)
        , symbolTex_(symbolTex)
        , caption_(caption)
        , filename_(filename)
        , isSelected_(isSelected)
    {

    }

    tgt::ivec2 TouchTablePreviewIcon::getLL() const {
        return tgt::ivec2(static_cast<int>(bounds_.getLLF().x),
            static_cast<int>(bounds_.getLLF().y));
    }

    void TouchTablePreviewIcon::setLL(tgt::ivec2 ll){
        tgt::vec3 ll3(tgt::vec2(ll),1.f);
        bounds_ = tgt::Bounds(ll3, tgt::vec3(ll.x + (float) width_, ll.y + (float) height_, 1.f));
    }

    int TouchTablePreviewIcon::getWidth() const { return width_; }

    void TouchTablePreviewIcon::setWidth(int width) {
        width_ = width;
        bounds_ = tgt::Bounds(bounds_.getLLF(), bounds_.getLLF() + tgt::vec3((float) width_, (float) height_, 0.f));
    }

    int TouchTablePreviewIcon::getHeight() const { return height_; }

    void TouchTablePreviewIcon::setHeight(int height){
        height_ = height;
        bounds_ = tgt::Bounds(bounds_.getLLF(), bounds_.getLLF() + tgt::vec3((float) width_, (float) height_, 0.f));
    }

    tgt::Texture* TouchTablePreviewIcon::getTexture() const { return symbolTex_; }

    void TouchTablePreviewIcon::setTexture(tgt::Texture* symbolTex) {
        symbolTex_ = symbolTex;
    }

    std::string TouchTablePreviewIcon::getCaption() const { return caption_; }

    void TouchTablePreviewIcon::setCaption(std::string caption){ caption_ = caption;}

    std::string TouchTablePreviewIcon::getFilename() const { return filename_; }

    void TouchTablePreviewIcon::setFilename(std::string filename) { filename_ = filename;}

    bool TouchTablePreviewIcon::isSelected() const { return isSelected_;}

    void TouchTablePreviewIcon::setIsSelected(bool isSelected){isSelected_ = isSelected;}

    bool TouchTablePreviewIcon::contains(tgt::ivec2 tp) {
        return bounds_.containsPoint(tgt::vec3(tgt::vec2(tp), 1.f));
    }

    bool TouchTablePreviewIcon::intersects(tgt::ivec2 ll, tgt::ivec2 ur) const {
        return bounds_.intersects(tgt::Bounds(tgt::vec3(tgt::vec2(ll), 1.f),
            tgt::vec3(tgt::vec2(ur), 1.f)));
    }
}
