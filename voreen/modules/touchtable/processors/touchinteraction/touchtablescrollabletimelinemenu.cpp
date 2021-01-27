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

#include "touchtablescrollabletimelinemenu.h"

namespace voreen{

    TouchTableScrollableTimelineMenu::TouchTableScrollableTimelineMenu(tgt::ivec2 anchor, tgt::vec4 color, int  duration, int stepSize, tgt::Font* font)
        :  TouchTableMenuFrame(anchor, color)
        , menuSlider_(tgt::ivec2(100,100), 10, 10, HORIZONTAL, 25, 25, 0)
        , sliderMargin_(tgt::ivec2(5,5))
        , menuMargin_(tgt::ivec2(10,10))
        , contentMenu_(anchor, tgt::vec4(1.f))
        , sliderID_(-1)
        , handle_(0)
        , font_(font)
        , contentOffset_(0)
        , timeline_(duration, stepSize)
    {
        //updateSliderMargin();
        placeSlider();
        placeContentMenu();
    }

    void TouchTableScrollableTimelineMenu::placeSlider(){

        tgt::ivec2 tmp = tgt::ivec2(menuSlider_.getIndicatorWidth()/2 + sliderMargin_.x, menuSlider_.getIndicatorHeight()/2 + sliderMargin_.y);
        menuSlider_.setBarOrigin(tgt::ivec2(static_cast<int>(bounds_.getLLF().x) + tmp.x, static_cast<int>(bounds_.getLLF().y) + tmp.y));
        menuSlider_.setBarLength(static_cast<int>(bounds_.getURB().x - bounds_.getLLF().x) - 2*tmp.x);
    }

    void TouchTableScrollableTimelineMenu::placeContentMenu(){

        tgt::ivec2 tmp = tgt::ivec2(menuSlider_.getIndicatorWidth()/2 + sliderMargin_.x, menuSlider_.getIndicatorHeight()/2 + sliderMargin_.y);
        contentMenu_.setLL(tgt::ivec2(static_cast<int>(bounds_.getLLF().x) + menuMargin_.x, static_cast<int>(bounds_.getLLF().y) + tmp.y *2));
        contentMenu_.setUR(tgt::ivec2(static_cast<int>(bounds_.getURB().x) - menuMargin_.x, static_cast<int>(bounds_.getURB().y) - menuMargin_.y));
        timeline_.setMainLine((contentMenu_.getUR().y - contentMenu_.getLL().y)/2);
    }

    void TouchTableScrollableTimelineMenu::setLL(tgt::ivec2 ll){

        tgt::vec3 ll3(tgt::vec2(ll),1.f);
        bounds_ = tgt::Bounds(ll3, tgt::max(bounds_.getURB(), ll3));
        placeSlider();
        placeContentMenu();
        computeContentMenuProperties();
    }

    void TouchTableScrollableTimelineMenu::setUR(tgt::ivec2 ur){

        tgt::vec3 ur3(tgt::vec3(tgt::vec2(ur),1.f));
        bounds_ = tgt::Bounds(tgt::min(bounds_.getLLF(), ur3) , ur3);
        placeSlider();
        placeContentMenu();
        computeContentMenuProperties();
    }


    void TouchTableScrollableTimelineMenu::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp){

        //iterate Touchpoints
        for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {

            tgt::vec2 tpPos = currentTp->pos();

            //check if TP state is pressed and hits indicator or a KeyValue
            if(currentTp->state() == tgt::TouchPoint::TouchPointPressed){
                int contentMenuXPos = static_cast<int>(tpPos.x) - contentMenu_.getLL().x;
                int contentMenuYPos = static_cast<int>(tpPos.y) - contentMenu_.getLL().y;
                //check if keyvalue has been hit
                KeyValue* keyValue = timeline_.keyValueHit(contentMenuXPos, contentMenuYPos);
                if(keyValue){
                    keyValueTouchPointMappping_.insert(std::make_pair(currentTp->id(), std::make_pair(keyValue->time_, keyValue)));
                    continue; //touchpoint has hit keyValue, dont do anything else
                }

                if(menuSlider_.checkIndicatorHit(tpPos))
                    //Indicator got hit -> save TPs ID
                    sliderID_ = currentTp->id();
                continue;
            }

            //if TP corresponding with indicator is moved update indicator position
            if((currentTp->state() == tgt::TouchPoint::TouchPointMoved) ){

                //check if moving touchpoint is associated with a keyValue
                std::map<int, std::pair<float, KeyValue*> >::iterator mappingIter = keyValueTouchPointMappping_.find(currentTp->id());

                if(mappingIter != keyValueTouchPointMappping_.end()
                    && tpPos.x >= contentMenu_.getLL().x + mappingIter->second.second->width_/2 && tpPos.x <= contentMenu_.getUR().x - mappingIter->second.second->width_/2){
                    mappingIter->second.second->pos_ = static_cast<int>(tpPos.x) - contentMenu_.getLL().x + timeline_.getTimeOffset();
                    //check if new position of keyValue is valid
                    if(timeline_.validateKeyValuePos(mappingIter->second.second->pos_, mappingIter->second.second->lastValidPos_)){
                        mappingIter->second.second->lastValidPos_ = mappingIter->second.second->pos_;
                    }
                    continue;
                }

                //check if moving touchpoint is associated with slider
                if( (currentTp->id() == sliderID_)){
                    menuSlider_.updateIndicatorPosition(tpPos);
                    //set timeline timeoffset: time user sees within the content menu depending on slider position
                    int contentMenuLength = contentMenu_.getUR().x - contentMenu_.getLL().x;
                    int lengthDiff = (timeline_.getDuration() > contentMenuLength)? timeline_.getDuration() - contentMenuLength : 0;
                    timeline_.setTimeOffset(static_cast<int>(menuSlider_.getIndicatorPos() * lengthDiff));
                }
            }


            if(currentTp->state() ==  tgt::TouchPoint::TouchPointReleased){
                if(currentTp->id() == sliderID_){
                    sliderID_ = -1;
                    continue;
                }else if(contentMenu_.contains(tpPos)){
                    tgt::TouchPoint tp = *currentTp;
                    handleContentMenu(tp);
                    continue;
                }
            }
        }
    }

    void TouchTableScrollableTimelineMenu::setSelectionHandler(boost::function<void (boost::tuple<float,bool,float>)> handle) {
        handle_ = handle;
    }

    void TouchTableScrollableTimelineMenu::handleContentMenu(tgt::TouchPoint& tp){
        int contentMenuXPos = static_cast<int>(tp.pos().x) - contentMenu_.getLL().x;
        int contentMenuYPos = static_cast<int>(tp.pos().y) - contentMenu_.getLL().y;

        bool keyValueHit = false; //true if key value has been hit -> different handling through "handle_"

        float oldKeyValueTime = -1;
        float newTime = -1;
        std::map<int, std::pair<float, KeyValue*> >::iterator mappingIter = keyValueTouchPointMappping_.find(tp.id());

        if(mappingIter != keyValueTouchPointMappping_.end()){
            mappingIter->second.second->pos_ = mappingIter->second.second->lastValidPos_;
            mappingIter->second.second->time_ = timeline_.getTimeForPos(mappingIter->second.second->pos_);
            keyValueHit = true;
            oldKeyValueTime = mappingIter->second.first;
            newTime = mappingIter->second.second->time_;
            keyValueTouchPointMappping_.erase(mappingIter);

        }else{
            //set current position in timeline
            timeline_.setCurrentPos(contentMenuXPos);
            timeline_.updateCurrentTime();
            newTime = timeline_.getCurrentTime();
        }
        if(handle_)
            handle_(boost::make_tuple(newTime, keyValueHit, oldKeyValueTime));

    }

    void TouchTableScrollableTimelineMenu::computeContentMenuProperties(){

        int textHeight = 1;
        if (font_ && (font_->getFontSize() > 0))
            textHeight = font_->getFontSize();

        //tgt::Bounds tmp = font_->getBounds(tgt::vec3(0,0,0), "KSP");
        //int textHeight = static_cast<int>(tmp.getURB().y - tmp.getLLF().y);
        //maxElementsInMenu_ = static_cast<int>((contentMenu_.getUR().y - contentMenu_.getLL().y) / textHeight);

    }

    void TouchTableScrollableTimelineMenu::setIndicatorPosition(float pos){
        menuSlider_.setIndicatorPos(pos);
        //set timeline timeoffset: time user sees within the content menu depending on slider position
        int contentMenuLength = contentMenu_.getUR().x - contentMenu_.getLL().x;
        int lengthDiff = (timeline_.getDuration() > contentMenuLength)? timeline_.getDuration() - contentMenuLength : 0;
        timeline_.setTimeOffset(static_cast<int>(menuSlider_.getIndicatorPos() * lengthDiff));
    }

    bool TouchTableScrollableTimelineMenu::addKeyValue(){
        int newKeyValue = timeline_.getCurrentPos();

        return timeline_.addKeyValueAt(newKeyValue);
    }

     TouchTableMenuFrame& TouchTableScrollableTimelineMenu::getContentMenu() {
        return contentMenu_;
    }

    int TouchTableScrollableTimelineMenu::getContentOffset() const {
        return timeline_.getTimeOffset();
    }


    TouchTableSlider& TouchTableScrollableTimelineMenu::getSlider() {
        return menuSlider_;
    }

    tgt::Font* TouchTableScrollableTimelineMenu::getFont() const {
        return font_;
    }

    void TouchTableScrollableTimelineMenu::setFont(tgt::Font* font){
        font_ = font;
        computeContentMenuProperties();
    }

    TouchTableTimeline& TouchTableScrollableTimelineMenu::getTimeline(){
        return timeline_;
    }

}
