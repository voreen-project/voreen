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

#include "touchtablemenuframe.h"

using tgt::ivec2;

namespace voreen {

/******************************** TouchTableMenuFrame*************************************/

     TouchTableMenuFrame:: TouchTableMenuFrame(ivec2 anchor, tgt::vec4 color)
        : anchor_(anchor)
        , color_(color)
    { }

    ivec2  TouchTableMenuFrame::getAnchor() const {
        return anchor_;
    }

    void  TouchTableMenuFrame::setAnchor(ivec2 anchor){
        anchor_ = anchor;
    }

    ivec2  TouchTableMenuFrame::getUR() const {
        return ivec2(static_cast<int>(bounds_.getURB().x),
            static_cast<int>(bounds_.getURB().y));
    }

    void  TouchTableMenuFrame::setUR(ivec2 ur){
        tgt::vec3 ur3(tgt::vec3(tgt::vec2(ur),1.f));
        bounds_ = tgt::Bounds(tgt::min(bounds_.getLLF(), ur3) , ur3);
    }

    ivec2  TouchTableMenuFrame::getLL() const {
        return ivec2(static_cast<int>(bounds_.getLLF().x),
            static_cast<int>(bounds_.getLLF().y));
    }

    void  TouchTableMenuFrame::setLL(ivec2 ll){
        tgt::vec3 ll3(tgt::vec2(ll),1.f);
        bounds_ = tgt::Bounds(ll3, tgt::max(bounds_.getURB(), ll3));
    }

    tgt::vec4  TouchTableMenuFrame::getColor() const {
        return color_;
    }
    void  TouchTableMenuFrame::setColor(tgt::vec4 color){
        color_ = color;
    }

    bool  TouchTableMenuFrame::contains(ivec2 tp) {
        return bounds_.containsPoint(tgt::vec3(tgt::vec2(tp), 1.f));
    }

    bool  TouchTableMenuFrame::intersects(ivec2 ll, ivec2 ur) const {
        return bounds_.intersects(tgt::Bounds(tgt::vec3(tgt::vec2(ll), 1.f),
            tgt::vec3(tgt::vec2(ur), 1.f)));
    }

/*************************************************************************************/

/***************************TouchTableClippingMenu************************************/
    TouchTableClippingMenu::TouchTableClippingMenu(ivec2 anchor, tgt::vec4 color)
        :  TouchTableMenuFrame(anchor, color)
    {
    }
/*************************************************************************************/

/***************************TouchTableOverlayMenu************************************/
    TouchTableOverlayMenu::TouchTableOverlayMenu(tgt::ivec2 anchor, tgt::vec4 color) :  TouchTableMenuFrame(anchor, color){

    }

    int TouchTableOverlayMenu::getID(){
        return id_;
    }

    void TouchTableOverlayMenu::setID(int id){
        id_ = id;
    }

    tgt::ivec2 TouchTableOverlayMenu::getButton(){
        return button_;
    }
    void TouchTableOverlayMenu::setButton(tgt::ivec2 button){
        button_ = button;
    }

    std::vector<tgt::ivec2>* TouchTableOverlayMenu::getButtons(){
        return &buttons_;
    }

    void TouchTableOverlayMenu::addButton(tgt::ivec2 button){
        buttons_.push_back(button);
    }

    void TouchTableOverlayMenu::clearButtons(){
        buttons_.clear();
    }


/*************************************************************************************/
}
