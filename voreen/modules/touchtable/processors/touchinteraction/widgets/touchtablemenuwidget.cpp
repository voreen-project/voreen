/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "touchtablemenuwidget.h"

namespace voreen {

TouchTableMenuWidget::TouchTableMenuWidget(std::string textureFile)
    : TouchTableWidget()
    , menu_(tgt::ivec2(100,100), tgt::vec4(0.4f, 0.4f, 0.4f, 0.5f))
    , rotationFrame_(tgt::ivec2(100,100), tgt::vec4(0.4f, 0.4f, 0.4f, 0.5f))
    , menuDimensions_("menudim", "Menu Dimensions", tgt::ivec2(650, 190), tgt::ivec2(20, 20), tgt::ivec2(1200, 850))
    , controlElementRadius_("elemradius", "Control Element Radius", 25, 10, 30)
    , invertMenuOrientation_("invertmenu", "Invert Menu Orientation", false)
    , lastOutportSize_(0,0)
    , symbolFileString_(textureFile)
    , rotate_(25, tgt::ivec2(50,50))
    , close_(25, tgt::ivec2(50,50))
    , rotateSymbol_(0)
    , closeSymbol_(0)
    , rotateButtonID_(-1)
    , closeButtonID_(-1)
{
    //some properties update menu coordinates
    addProperty(menuDimensions_);
    menuDimensions_.onChange(MemberFunctionCallback<TouchTableMenuWidget>(this, &TouchTableMenuWidget::updateMenuCoordinates));
    position_.onChange(MemberFunctionCallback<TouchTableMenuWidget>(this, &TouchTableMenuWidget::updateMenuCoordinates));
    isOnTable_.onChange(MemberFunctionCallback<TouchTableMenuWidget>(this, &TouchTableMenuWidget::updateMenuCoordinates));

    //update control elements if radius changes
    addProperty(controlElementRadius_);
    controlElementRadius_.onChange(MemberFunctionCallback<TouchTableMenuWidget>(this, &TouchTableMenuWidget::updateComponentAttributes));

    addProperty(invertMenuOrientation_);
    invertMenuOrientation_.onChange(MemberFunctionCallback<TouchTableMenuWidget>(this, &TouchTableMenuWidget::updateComponentAttributes));

    //symbol should not be changed: set property invisible
    symbolFile_.setVisibleFlag(false);
}

void TouchTableMenuWidget::process() {

    if (!overlay_)
        return;

    if (overlay_->getOutportSize() != lastOutportSize_) {
        lastOutportSize_ = overlay_->getOutportSize();
        //since outport size changed: compute menu position if widget is not on the table
        if (!isOnTable_.get())
            updateMenuCoordinates();
    }

    updateComponentAttributes();
}

void TouchTableMenuWidget::render() {

    if (!overlay_)
        return;

    //render menu frame
    overlay_->renderMenuFrame(menu_, false);

    //render rotation frame and rotate button
    overlay_->renderMenuFrame(rotationFrame_, false);
    overlay_->renderControlElement(&rotate_, rotationFrame_.getLL(), 0.f, false);
    overlay_->renderControlElement(&close_, rotationFrame_.getLL(), 0.f, false);

    renderComponents();
}

void TouchTableMenuWidget::pressed() {
    if(isActive()){
        if (overlay_)
            overlay_->setMenuWidget(0);
        setActive(false);
    }else{
        if (overlay_)
            overlay_->setMenuWidget(this);
        updateMenuCoordinates();
        setActive();
    }
}

void TouchTableMenuWidget::initialize() {
    TouchTableWidget::initialize();

    //dispose texture loaded by subclass
    TexMgr.dispose(symbolTex_); symbolTex_ = 0;

    //set texture
    TexMgr.addPath(VoreenApplication::app()->getModulePath("touchtable") + "/textures");
    symbolTex_ = TexMgr.load(symbolFileString_);
    symbolFile_.set(VoreenApplication::app()->getModulePath("touchtable") + "/textures" + "/" + symbolFileString_);

    rotateSymbol_ = TexMgr.load("rotate.png");
    rotate_.setSymbolTexture(rotateSymbol_);

    closeSymbol_ = TexMgr.load("exit.png");
    close_.setSymbolTexture(closeSymbol_);
}

void TouchTableMenuWidget::deinitialize() {
    if (overlay_ && isActive_.get()) {
        overlay_->setMenuWidget(0);
    }

    rotate_.setSymbolTexture(0);
    TexMgr.dispose(rotateSymbol_);
    rotateSymbol_ = 0;

    close_.setSymbolTexture(0);
    TexMgr.dispose(closeSymbol_);
    closeSymbol_ = 0;

    TouchTableWidget::deinitialize();
}

void TouchTableMenuWidget::setOverlay(TouchTableOverlay* overlay) {
    TouchTableWidget::setOverlay(overlay);
    if (overlay && isActive_.get())
        overlay->setMenuWidget(this);
    else
        isActive_.set(false);

    if (overlay)
        lastOutportSize_ = overlay->getOutportSize();

    updateMenuCoordinates();
}

void TouchTableMenuWidget::updateMenuCoordinates() {

    if (!overlay_)
        return;

    int widgetRadius = overlay_->getWidgetRadius();
    tgt::ivec2 outportSize = overlay_->getOutportSize();

    tgt::ivec2 frameDimensions = menuDimensions_.get() + tgt::ivec2(rotationFrame_.getUR().x - rotationFrame_.getLL().x, 0);

    if (!isOnTable_.get()) {
        //if the widget is located within the overlay menu: render the menu frame at the bottom of the screen
        menu_.setAnchor(getPosition());
        menu_.setLL(tgt::ivec2((outportSize.x - frameDimensions.x) / 2, 0));
        menu_.setUR(menu_.getLL() + menuDimensions_.get());

        //set rotate frame to the right
        rotationFrame_.setLL(tgt::ivec2(menu_.getUR().x, menu_.getLL().y));
        rotationFrame_.setUR(tgt::ivec2(rotationFrame_.getLL().x + controlElementRadius_.get() * 2 + 20, menu_.getUR().y));
    }
    else {

        menu_.setAnchor(getPosition());

        //compute new position
        int llx, urx, lly, ury;

        if (menuPlacement_.horizontal_ == MenuPlacement::RIGHT)
            llx = static_cast<int>(getPosition().x) + widgetRadius;
        else    //LEFT
            llx = static_cast<int>(getPosition().x) - frameDimensions.x - widgetRadius;

        urx = llx + menuDimensions_.get().x;

        if (menuPlacement_.vertical_ == MenuPlacement::UP)
            lly = static_cast<int>(getPosition().y) + widgetRadius;
        else    //DOWN
            lly = static_cast<int>(getPosition().y) - menuDimensions_.get().y - widgetRadius;

        ury = lly + menuDimensions_.get().y;

        //set coordinates
        menu_.setLL(tgt::ivec2(llx, lly));
        menu_.setUR(tgt::ivec2(urx, ury));
        rotationFrame_.setLL(tgt::ivec2(menu_.getUR().x, menu_.getLL().y));
        rotationFrame_.setUR(tgt::ivec2(rotationFrame_.getLL().x + controlElementRadius_.get() * 2 + 20, menu_.getUR().y));

        //check placement
        bool placementChanged = false;
        if ((menuPlacement_.horizontal_ == MenuPlacement::RIGHT) && (rotationFrame_.getUR().x > outportSize.x)) {
            menuPlacement_.horizontal_ = MenuPlacement::LEFT;
            placementChanged = true;
        }

        if ((menuPlacement_.horizontal_ == MenuPlacement::LEFT) && (menu_.getLL().x < 0)) {
            menuPlacement_.horizontal_ = MenuPlacement::RIGHT;
            placementChanged = true;
        }

        if ((menuPlacement_.vertical_ == MenuPlacement::UP) && (rotationFrame_.getUR().y > outportSize.y)) {
            menuPlacement_.vertical_ = MenuPlacement::DOWN;
            placementChanged = true;
        }

        if ((menuPlacement_.vertical_ == MenuPlacement::DOWN) && (menu_.getLL().y < 0)) {
            menuPlacement_.vertical_ = MenuPlacement::UP;
            placementChanged = true;
        }

        //compute new position if necessary
        if (placementChanged) {

            if (menuPlacement_.horizontal_ == MenuPlacement::RIGHT)
                llx = static_cast<int>(getPosition().x) + widgetRadius;
            else    //LEFT
                llx = static_cast<int>(getPosition().x) - frameDimensions.x - widgetRadius;

            urx = llx + menuDimensions_.get().x;

            if (menuPlacement_.vertical_ == MenuPlacement::UP)
                lly = static_cast<int>(getPosition().y) + widgetRadius;
            else    //DOWN
                lly = static_cast<int>(getPosition().y) - menuDimensions_.get().y - widgetRadius;

            ury = lly + menuDimensions_.get().y;

            //set coordinates
            menu_.setLL(tgt::ivec2(llx, lly));
            menu_.setUR(tgt::ivec2(urx, ury));
            rotationFrame_.setLL(tgt::ivec2(menu_.getUR().x, menu_.getLL().y));
            rotationFrame_.setUR(tgt::ivec2(rotationFrame_.getLL().x + controlElementRadius_.get() * 2 + 20, menu_.getUR().y));
        }
    }

    //also update position and status of control elements and other components
    updateComponentAttributes();
}

void TouchTableMenuWidget::updateComponentAttributes() {
    rotate_.setRadius(controlElementRadius_.get());
    rotate_.setActive(invertMenuOrientation_.get());

    close_.setRadius(controlElementRadius_.get());
    close_.setActive(false);

    int posX = (rotationFrame_.getUR().x - rotationFrame_.getLL().x) / 2;
    int posYUp = rotationFrame_.getUR().y - rotationFrame_.getLL().y - 10 - controlElementRadius_.get();
    int posYDown = 10 + controlElementRadius_.get();

    /*if (invertMenuOrientation_.get()) {
        close_.setPosition(tgt::ivec2(posX,posYDown));
        rotate_.setPosition(tgt::ivec2(posX, posYUp));
    }
    else {
        close_.setPosition(tgt::ivec2(posX,posYUp));
        rotate_.setPosition(tgt::ivec2(posX, posYDown));
    }*/

    close_.setPosition(tgt::ivec2(posX,posYUp));
    rotate_.setPosition(tgt::ivec2(posX, posYDown));

    updateComponents();
}

bool TouchTableMenuWidget::isInMenu(tgt::ivec2 tp) {
    return (isActive() && (menu_.contains(tp) || rotationFrame_.contains(tp)));
}

void TouchTableMenuWidget::handleWidgetTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

    //check if a touch point presses the rotate button, remember the ID and do not foward it to the subclass method
    //when released: remove the remembered ID

    std::deque<tgt::TouchPoint> pointsForSubclass;

    for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {

        //convert position according to viewport
        tgt::ivec2 viewport(0,0);
        if (overlay_)
            viewport = overlay_->getOutportSize();
        tgt::ivec2 pos = tgt::ivec2(static_cast<int>(currentTp->pos().x), viewport.y - static_cast<int>(currentTp->pos().y));

        if (currentTp->state() == tgt::TouchPoint::TouchPointPressed) {
            //check if the touch point is within the rotation frame and hits the rotate button
            if (rotationFrame_.contains(pos)) {
                if (hitsControlElement(pos - rotationFrame_.getLL(), rotate_))
                    rotateButtonID_ = currentTp->id();
                else if (hitsControlElement(pos - rotationFrame_.getLL(), close_))
                    closeButtonID_ = currentTp->id();
            }
            else
                pointsForSubclass.push_back(*currentTp);

            continue;
        }

        if (currentTp->state() == tgt::TouchPoint::TouchPointMoved) {
            if (((rotateButtonID_ == -1) || (currentTp->id() != rotateButtonID_)) && ((closeButtonID_ == -1) || (currentTp->id() != closeButtonID_)))
                pointsForSubclass.push_back(*currentTp);

            continue;
        }

        if (currentTp->state() == tgt::TouchPoint::TouchPointReleased) {
            if (currentTp->id() == rotateButtonID_) {
                rotateButtonID_ = -1;
                if (hitsControlElement(pos - rotationFrame_.getLL(), rotate_))
                    invertMenuOrientation_.set(!invertMenuOrientation_.get());
            }
            else if (currentTp->id() == closeButtonID_) {
                closeButtonID_ = -1;
                if (hitsControlElement(pos - rotationFrame_.getLL(), close_))
                    pressed();
            }
            else
                pointsForSubclass.push_back(*currentTp);

            continue;
        }
    }

    handleTouchPoints(pointsForSubclass);
}

bool TouchTableMenuWidget::isMenuInverted() const {
    return invertMenuOrientation_.get();
}

tgt::ivec2 TouchTableMenuWidget::getMenuDimensions() const {
    return menuDimensions_.get();
}

tgt::ivec2 TouchTableMenuWidget::convertToMenuCoordinates(tgt::ivec2 p) const {
    tgt::ivec2 viewport(0,0);
    if (overlay_)
        viewport = overlay_->getOutportSize();

    //transform y-coordinate
    tgt::ivec2 pos = tgt::ivec2(p.x, viewport.y - p.y);

    //transform coordinates into menu coordinate system
    pos = pos - menu_.getLL();

    //if menu orientation is inverted: invert coordinates
    if (invertMenuOrientation_.get())
        pos = menuDimensions_.get() - pos;

    return pos;
}

std::deque<tgt::TouchPoint> TouchTableMenuWidget::convertTouchPointsToMenuCoordinates(const std::deque<tgt::TouchPoint>& points) const {
    std::deque<tgt::TouchPoint> convertedPoints;

    for (std::deque<tgt::TouchPoint>::const_iterator i = points.begin(); i != points.end(); ++i) {
        tgt::TouchPoint p = *i;
        p.setPos(tgt::vec2(convertToMenuCoordinates(tgt::ivec2(p.pos()))));
        convertedPoints.push_back(p);
    }

    return convertedPoints;
}

tgt::ivec2 TouchTableMenuWidget::getMenuLL() const {
    return menu_.getLL();
}

bool TouchTableMenuWidget::hitsControlElement(tgt::ivec2 tp, const TouchTableControlElement& element) const {
    return (tgt::length(tgt::vec2(tp) - tgt::vec2(element.getPosition())) <= static_cast<float>(controlElementRadius_.get()));
}

} // namespace
