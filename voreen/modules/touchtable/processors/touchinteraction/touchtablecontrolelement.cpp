/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "touchtablecontrolelement.h"
#include "tgt/texturemanager.h"
#include "voreen/core/voreenapplication.h"



namespace voreen {

    TouchTableControlElement::TouchTableControlElement()
    {

    }

    TouchTableControlElement::TouchTableControlElement(int radius, tgt::ivec2 pos, tgt::Texture* symbolTexture,
        float opacity, tgt::vec4 colorMod, tgt::vec4 propColor, int id, bool usePropColor, bool isActivatable)
        : radius_(radius)
        , pos_(pos)
        , symbolTexture_(symbolTexture)
        , opacity_(opacity)
        , colorMod_(colorMod)
        , propColor_(propColor)
        , usePropColor_(usePropColor)
        , id_(id)
        , active_(false)
        , isActivatable_(isActivatable)
    {
        //TexMgr.addPath(VoreenApplication::app()->getModulePath("touchtable") + "/textures");
        //backGround_ = TexMgr.load("newPuckOn2.png");
    }

    int TouchTableControlElement::getId() const {
        return id_;
    }

    void TouchTableControlElement::setColorMod(tgt::vec4 colorMod){
        colorMod_ = colorMod;
    }

    tgt::vec4 TouchTableControlElement::getColorMod() const {
        return colorMod_;
    }

    bool TouchTableControlElement::isActive() const {
        return active_;
    }

    void TouchTableControlElement::setActive(bool active){
        active_ = active;
    }

    tgt::Texture* TouchTableControlElement::getSymbolTexture() const {
        return symbolTexture_;
    }

    void TouchTableControlElement::setSymbolTexture(tgt::Texture* symbolTexture){
        symbolTexture_ = symbolTexture;
    }

    float TouchTableControlElement::getOpacity() const {
        return opacity_;
    }

    void TouchTableControlElement::setOpacity(float opacity) {
        opacity_ = opacity;
    }

    tgt::ivec2 TouchTableControlElement::getPos() const {
        return pos_;
    }

    int TouchTableControlElement::getRadius() const {
        return radius_;
    }

    void TouchTableControlElement::setPosition(tgt::ivec2 pos) {
        pos_ = pos;
    }

    void TouchTableControlElement::setRadius(int radius) {
        radius_ = radius;
    }

    tgt::ivec2 TouchTableControlElement::getPosition() const {
        return pos_;
    }

    bool TouchTableControlElement::isActivatable() const {
        return isActivatable_;
    }

} // namespace
