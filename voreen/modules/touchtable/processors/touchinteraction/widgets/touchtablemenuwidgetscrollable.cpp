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

#include "touchtablemenuwidgetscrollable.h"
#include "tgt/immediatemode/immediatemode.h"
namespace voreen {

TouchTableMenuWidgetScrollable::TouchTableMenuWidgetScrollable(std::string textureFile)
    : TouchTableMenuWidget(textureFile)
    , scrollableMenuIsOpen_("scrollablemenuopen", "Scrollable Menu Open", true)
    , fontProp_("font", "Font")
    , scrollableMenu_(tgt::ivec2(100,100), tgt::vec4(1.f))
{
    addProperty(scrollableMenuIsOpen_);
    scrollableMenuIsOpen_.setVisibleFlag(false);

    addProperty(fontProp_);

    //set default value
    fontProp_.get()->setFontSize(22);
    fontProp_.get()->setFontColor(tgt::vec4(0.65f,0.65f,0.65f,1.f));
}

void TouchTableMenuWidgetScrollable::render() {

    if (!overlay_)
        return;

    //call parent's class render method
    TouchTableMenuWidget::render();

    //render scrollable menu if open
    if (scrollableMenuIsOpen_.get())
        overlay_->renderScrollableMenu(&scrollableMenu_, getMenuLL());
}

void TouchTableMenuWidgetScrollable::initialize() {
    TouchTableMenuWidget::initialize();

    scrollableMenu_.setSelectionHandler(std::bind1st(std::mem_fun(&TouchTableMenuWidgetScrollable::handleScrollableMenuSelection), this));
}

void TouchTableMenuWidgetScrollable::deinitialize() {
    TouchTableMenuWidget::deinitialize();
}

void TouchTableMenuWidgetScrollable::updateComponentAttributes() {
    TouchTableMenuWidget::updateComponentAttributes();

    scrollableMenu_.setSelectionHandler(std::bind1st(std::mem_fun(&TouchTableMenuWidgetScrollable::handleScrollableMenuSelection), this));
    updateScrollableMenuPosition();

    scrollableMenu_.placeContentMenu();
    scrollableMenu_.placeSlider();

    fontProp_.get()->setLineWidth(static_cast<float>(scrollableMenu_.getTextColumnWidth()));
    fontProp_.invalidate();

    scrollableMenu_.setFont(fontProp_.get());
}

void TouchTableMenuWidgetScrollable::handleWidgetTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

    if (scrollableMenuIsOpen_.get()) {
        std::deque<tgt::TouchPoint> pointsForMenuWidget;
        std::deque<tgt::TouchPoint> pointsForScrollableMenu;

        for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {

            if (currentTp->state() == tgt::TouchPoint::TouchPointPressed) {
                //check if preset menu contains touch point
                if (scrollableMenu_.contains(convertToMenuCoordinates(tgt::ivec2(currentTp->pos())))) {
                    pointsForScrollableMenu.push_back(*currentTp);
                    scrollableMenuIDs_.push_back(currentTp->id());
                }
                else {
                    pointsForMenuWidget.push_back(*currentTp);
                }
            }
            else if (currentTp->state() == tgt::TouchPoint::TouchPointMoved) {
                //check if touch point is associated with the scrollable menu
                if (std::find(scrollableMenuIDs_.begin(), scrollableMenuIDs_.end(), currentTp->id()) != scrollableMenuIDs_.end())
                    pointsForScrollableMenu.push_back(*currentTp);
                else
                    pointsForMenuWidget.push_back(*currentTp);

            }
            else if (currentTp->state() == tgt::TouchPoint::TouchPointReleased) {
                //check if tp is associated with the scrollable menu
                std::vector<int>::iterator i = std::find(scrollableMenuIDs_.begin(), scrollableMenuIDs_.end(), currentTp->id());
                if (i != scrollableMenuIDs_.end()) {
                    pointsForScrollableMenu.push_back(*currentTp);
                    scrollableMenuIDs_.erase(i);
                }
                else
                    pointsForMenuWidget.push_back(*currentTp);
            }
            else
                pointsForMenuWidget.push_back(*currentTp);
        }

        if (!pointsForScrollableMenu.empty()) {
            scrollableMenu_.handleTouchPoints(convertTouchPointsToMenuCoordinates(pointsForScrollableMenu));
            invalidate();
        }

        TouchTableMenuWidget::handleWidgetTouchPoints(pointsForMenuWidget);
    }
    else
        TouchTableMenuWidget::handleWidgetTouchPoints(tp);
}

} // namespace
