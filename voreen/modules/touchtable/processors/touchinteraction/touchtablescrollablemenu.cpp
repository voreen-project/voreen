/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "touchtablescrollablemenu.h"

namespace voreen{

    TouchTableScrollableMenu::TouchTableScrollableMenu(tgt::ivec2 anchor, tgt::vec4 color, tgt::Font* font)
        :  TouchTableMenuFrame(anchor, color)
        , menuSlider_(tgt::ivec2(100,100), 10, 10, VERTICAL, 25, 25, 0)
        , sliderMargin_(tgt::ivec2(5,5))
        , menuMargin_(tgt::ivec2(10,10))
        , contentMenu_(anchor, tgt::vec4(1.f))
        , sliderID_(-1)
        , handle_(0)
        , font_(font)
        , maxElementsInMenu_(4)
        , contentOffset_(0)
        , rowHeight_(25)
        , isVisualizingSelection_(true)
    {
        if(font_){
            font_->setFontSize(20);
            font_->setLineWidth(400);
            font_->setTextAlignment(tgt::Font::Centered);
        }

        placeSlider();

    }

    void TouchTableScrollableMenu::placeSlider(){

        tgt::ivec2 tmp = tgt::ivec2(menuSlider_.getIndicatorWidth()/2 + sliderMargin_.x, menuSlider_.getIndicatorHeight()/2 + sliderMargin_.y);
        menuSlider_.setBarOrigin(tgt::ivec2(static_cast<int>(bounds_.getURB().x) - tmp.x, static_cast<int>(bounds_.getURB().y) - tmp.y));
        menuSlider_.setBarLength(static_cast<int>(bounds_.getURB().y - bounds_.getLLF().y) - 2*tmp.y);
    }

    void TouchTableScrollableMenu::placeContentMenu(){

        tgt::ivec2 tmp = tgt::ivec2(menuSlider_.getIndicatorWidth()/2 + sliderMargin_.x, menuSlider_.getIndicatorHeight()/2 + sliderMargin_.y);
        contentMenu_.setLL(tgt::ivec2(static_cast<int>(bounds_.getLLF().x) + menuMargin_.x, static_cast<int>(bounds_.getLLF().y) + menuMargin_.y));
        contentMenu_.setUR(tgt::ivec2(static_cast<int>(bounds_.getURB().x) - (tmp.x * 2), static_cast<int>(bounds_.getURB().y) - menuMargin_.y));
    }

    void TouchTableScrollableMenu::setLL(tgt::ivec2 ll){

        tgt::vec3 ll3(tgt::vec2(ll),1.f);
        bounds_ = tgt::Bounds(ll3, tgt::max(bounds_.getURB(), ll3));
        placeSlider();
        placeContentMenu();
        computeContentMenuProperties();
    }

    void TouchTableScrollableMenu::setUR(tgt::ivec2 ur){

        tgt::vec3 ur3(tgt::vec3(tgt::vec2(ur),1.f));
        bounds_ = tgt::Bounds(tgt::min(bounds_.getLLF(), ur3) , ur3);
        placeSlider();
        placeContentMenu();
        computeContentMenuProperties();
    }

    void TouchTableScrollableMenu::setContent(std::vector<std::pair<std::string, tgt::Texture*> > content){
        content_ = content;
        menuSlider_.setIndicatorPos(0.f);
        contentOffset_ = 0;
    }

    std::vector<std::pair<std::string, tgt::Texture*> > TouchTableScrollableMenu::getContent() const {
        return content_;
    }

    void TouchTableScrollableMenu::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp){

        //iterate Touchpoints
        for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {

            tgt::vec2 pos = currentTp->pos();

            //check if TP state is pressed and hits indicator
            if(currentTp->state() == tgt::TouchPoint::TouchPointPressed){
                if(menuSlider_.checkIndicatorHit(pos))
                    //Indicator got hit -> save TPs ID
                    sliderID_ = currentTp->id();
                continue;
            }

            //if TP corresponding with indicator is moved update indicator position
            if((currentTp->state() == tgt::TouchPoint::TouchPointMoved) && (currentTp->id() == sliderID_)){
                if(content_.size() > maxElementsInMenu_){
                    menuSlider_.updateIndicatorPosition(pos);
                    int tmp = static_cast<int>(content_.size() - maxElementsInMenu_) * menuSlider_.getIndicatorPos();
                    contentOffset_ = tgt::clamp(tmp, 0, static_cast<int>(content_.size()));
                }
                continue;
                //TODO: write function that scrolls contentMenu
            }


            if(currentTp->state() ==  tgt::TouchPoint::TouchPointReleased){
                if(currentTp->id() == sliderID_){
                    sliderID_ = -1;
                    continue;
                }else if(contentMenu_.contains(pos)){
                    handleContentMenu(pos);
                    continue;
                }
            }
        }
    }

    void TouchTableScrollableMenu::setSelectionHandler(boost::function<void (std::string)> handle) {
        handle_ = handle;
    }

    void TouchTableScrollableMenu::handleContentMenu(tgt::vec2 pos){


        int index = tgt::clamp(static_cast<int>((contentMenu_.getUR().y - pos.y) / static_cast<float>(rowHeight_)),
            0, static_cast<int>(content_.size()-1));
        index = index + contentOffset_;
        selectedRow_ = index;
        if (handle_)
            handle_(content_.at(index).first);
    }

    void TouchTableScrollableMenu::computeContentMenuProperties(){

        int textHeight = 1;
        if (font_ && (font_->getFontSize() > 0))
            textHeight = font_->getFontSize();

        if(textHeight > rowHeight_)
            font_->setFontSize(rowHeight_);

        maxElementsInMenu_ = static_cast<int>((contentMenu_.getUR().y - contentMenu_.getLL().y) / rowHeight_);

    }

     TouchTableMenuFrame& TouchTableScrollableMenu::getContentMenu() {
        return contentMenu_;
    }

    int TouchTableScrollableMenu::getContentOffset() const {
        return contentOffset_;
    }

    void TouchTableScrollableMenu::setContentOffset(int offset){
        contentOffset_ = offset;
    }

    int TouchTableScrollableMenu::getMaxElementsInMenu() const {
        return maxElementsInMenu_;
    }

    TouchTableSlider& TouchTableScrollableMenu::getSlider() {
        return menuSlider_;
    }

    tgt::Font* TouchTableScrollableMenu::getFont() const {
        return font_;
    }

    void TouchTableScrollableMenu::setFont(tgt::Font* font){
        font_ = font;
        computeContentMenuProperties();
    }

    int TouchTableScrollableMenu::getTextColumnWidth() const {
        return textColumnWidth_;
    }

    void TouchTableScrollableMenu::setTextColumnWidth(int width){
        textColumnWidth_ = width;
    }

    int TouchTableScrollableMenu::getRowHeight() const {
        return rowHeight_;
    }

    void TouchTableScrollableMenu::setRowHeight(int height){
        rowHeight_ = height;
    }

    void TouchTableScrollableMenu::setVisualizingSelection(bool selection){
        isVisualizingSelection_ = selection;
    }
    bool TouchTableScrollableMenu::isVisualizingSelection(){
        return isVisualizingSelection_;
    }

    int TouchTableScrollableMenu::getSelectedRow(){
        return selectedRow_;
    }

    void TouchTableScrollableMenu::setSelectedRow(int row){
        selectedRow_ = row;
    }

    void TouchTableScrollableMenu::clearSelection() {
        setSelectedRow(-1);
    }
}
