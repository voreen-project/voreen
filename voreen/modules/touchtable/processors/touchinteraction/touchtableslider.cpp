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

#include "touchtableslider.h"

namespace voreen {

    const std::string TouchTableSlider::loggerCat_("voreen.touchtable.TouchTableSlider");
    const std::string TouchTableTwoIndicatorSlider::loggerCat_("voreen.touchtable.TouchTableTwoIndicatorSlider");

/********************************TouchTableSlider**************************************/

    TouchTableSlider::TouchTableSlider(tgt::ivec2 origin, int length, int width, Orientation orientation, int indHeight, int indWidth, float indPos)
        /*: barOrigin_(origin)
        , barLength_(length)
        , barWidth_(width)
        , orientation_(orientation)
        , indicatorPos_(indPos)
        , indicatorHeight_(indHeight)
        , indicatorWidth_(indWidth)*/
    {
        indicator_.height_= indHeight;
        indicator_.width_ = indWidth;
        indicator_.position_ = indPos;
        bar_.length_ = length;
        bar_.orientation_ = orientation;
        bar_.origin_ = origin;
        bar_.width_ = width;
        updateIndicatorBounds();
    }

    void TouchTableSlider::setIndicatorPos(float pos){

        indicator_.position_ = pos;
        updateIndicatorBounds();
    }

    float TouchTableSlider::getIndicatorPos() const{
        return indicator_.position_;
    }

    tgt::ivec2 TouchTableSlider::getBarOrigin() const{
        return bar_.origin_;
    }

    void TouchTableSlider::setBarOrigin(tgt::ivec2 origin){
        bar_.origin_ = origin;
        updateIndicatorBounds();
    }

    void TouchTableSlider::setBarLength(int length){
        bar_.length_ = length;
        updateIndicatorBounds();
    }

    int TouchTableSlider::getBarLength() const{
        return bar_.length_;
    }

    void TouchTableSlider::setOrientation(Orientation orientation){
        bar_.orientation_ = orientation;
    }

    Orientation TouchTableSlider::getOrientation() const {
        return bar_.orientation_;
    }

    void TouchTableSlider::updateIndicatorPosition(tgt::ivec2 newPos){

        if(bar_.orientation_ == HORIZONTAL)
            newPos.y = bar_.origin_.y;
        else
            newPos.x = bar_.origin_.x;
        tgt::vec2 tmp = (tgt::vec2(newPos - bar_.origin_))/static_cast<float>(bar_.length_);
        indicator_.position_ = (float)(tmp.x + tmp.y);
        if (bar_.orientation_ == VERTICAL)
            indicator_.position_ = indicator_.position_ * -1;
        indicator_.position_ = tgt::clamp(indicator_.position_, 0.f, 1.f);
        updateIndicatorBounds();
    }

    void TouchTableSlider::setIndicatorHeight(int height){

        indicator_.height_ = height;
        updateIndicatorBounds();
    }

    void TouchTableSlider::setIndicatorWidth(int width){

        indicator_.width_ = width;
        updateIndicatorBounds();
    }

    void TouchTableSlider::updateIndicatorBounds(){

        tgt::ivec2 indPos = interpolate();
        indicator_.bounds_ = tgt::Bounds(tgt::vec3(indPos.x - indicator_.width_/2, indPos.y - indicator_.height_/2, 1.f), tgt::vec3(indPos.x + indicator_.width_/2, indPos.y + indicator_.height_/2, 1.f));

    }

    tgt::ivec2 TouchTableSlider::getIndicatorLL() const {

        return tgt::ivec2(static_cast<int>(indicator_.bounds_.getLLF().x),
            static_cast<int>(indicator_.bounds_.getLLF().y));
    }

    tgt::ivec2 TouchTableSlider::getIndicatorUR() const {
        return tgt::ivec2(static_cast<int>(indicator_.bounds_.getURB().x),
            static_cast<int>(indicator_.bounds_.getURB().y));
    }

    tgt::ivec2 TouchTableSlider::getBarLL() const {

        if (bar_.orientation_ == HORIZONTAL)
            return tgt::ivec2(bar_.origin_.x, bar_.origin_.y-bar_.width_/2);
        else
            return tgt::ivec2(bar_.origin_.x -bar_.width_/2, bar_.origin_.y - bar_.length_);
    }

    tgt::ivec2 TouchTableSlider::getBarUR() const {

        if (bar_.orientation_ == HORIZONTAL)
            return tgt::ivec2(bar_.origin_.x + bar_.length_, bar_.origin_.y + bar_.width_/2);
        else
            return tgt::ivec2(bar_.origin_.x + bar_.width_/2, bar_.origin_.y);
    }

    tgt::ivec2 TouchTableSlider::interpolate(){

        tgt::vec2 tmp;
        if (bar_.orientation_ == HORIZONTAL){
            tmp = tgt::vec2(indicator_.position_ * bar_.length_,0);
            return bar_.origin_ + tgt::ivec2(static_cast<int>(tmp.x), static_cast<int>(tmp.y));
        }else{
            tmp = tgt::vec2(0, indicator_.position_ * bar_.length_);
            return bar_.origin_ - tgt::ivec2(static_cast<int>(tmp.x), static_cast<int>(tmp.y));
        }


    }

    bool TouchTableSlider::checkIndicatorHit(tgt::ivec2 tp){

        if(indicator_.bounds_.containsPoint(tgt::vec3(tp, 1.f))){
            //isIndicatorHit_ = true;
            return true;
        }
        else
            return false;
    }

    int TouchTableSlider::getIndicatorHeight() const {
        return indicator_.height_;
    }

    int TouchTableSlider::getIndicatorWidth() const {
        return indicator_.width_;
    }

    /*void TouchTableSlider::unHitIndicator(){
        isIndicatorHit_ = false;
    }

    int TouchTableSlider::getIndicatorHeight();
    int TouchTableSlider::getIndicatorWidth();
    */

    int TouchTableSlider::getBarWidth() const {
        return bar_.width_;
    }

    void TouchTableSlider::setBarWidth(int width) {
        bar_.width_ = width;
    }

/*************************************************************************************/

/**************************TouchTableTwoIndicatorSlider******************************/

    TouchTableTwoIndicatorSlider::TouchTableTwoIndicatorSlider(tgt::ivec2 origin, int length, int width, Orientation orientation, int indHeight, int indWidth)
        : TouchTableSlider(origin, length, width, orientation, 0, 0, 0.f)
    {
        leftIndicator_.height_= indHeight;
        leftIndicator_.width_ = indWidth;
        leftIndicator_.position_ = 0;
        leftIndicator_.isLeft_ = true;
        updateIndicatorBounds(&leftIndicator_);

        rightIndicator_.height_= indHeight;
        rightIndicator_.width_ = indWidth;
        rightIndicator_.position_ = 1;
        rightIndicator_.isLeft_ = false;
        updateIndicatorBounds(&rightIndicator_);

        bar_.length_ = length;
        bar_.orientation_ = orientation;
        bar_.origin_ = origin;
        bar_.width_ = width;

    }

    void TouchTableTwoIndicatorSlider::setIndicatorPositions(tgt::vec2 pos){
        leftIndicator_.position_= std::min(pos.x, pos.y);
        rightIndicator_.position_= std::max(pos.x, pos.y);
        updateIndicatorBounds(&leftIndicator_);
        updateIndicatorBounds(&rightIndicator_);
    }

    void TouchTableTwoIndicatorSlider::updateIndicatorBounds(Indicator* indicator){

        tgt::ivec2 indPos = interpolate(indicator);
        if(indicator->isLeft_)
            indicator->bounds_ = tgt::Bounds(tgt::vec3(indPos.x - indicator->width_, indPos.y - indicator->height_/2, 1.f), tgt::vec3(indPos.x, indPos.y + indicator->height_/2, 1.f));
        else
            indicator->bounds_ = tgt::Bounds(tgt::vec3(indPos.x, indPos.y - indicator->height_/2, 1.f), tgt::vec3(indPos.x + indicator->width_, indPos.y + indicator->height_/2, 1.f));
    }

    bool TouchTableTwoIndicatorSlider::checkIndicatorHit(tgt::vec2 pos, int id){

        if(leftIndicator_.bounds_.containsPoint(tgt::vec3(pos, 1.f))){
            leftIndicator_.associatedTouchPointID_ = id;
            return true;
        }
        else if(rightIndicator_.bounds_.containsPoint(tgt::vec3(pos, 1.f))){
            rightIndicator_.associatedTouchPointID_ = id;

            return true;
        }else
            return false;
    }

    tgt::ivec2 TouchTableTwoIndicatorSlider::interpolate(Indicator* indicator){

        tgt::vec2 tmp;
        if (bar_.orientation_ == HORIZONTAL)
            tmp = tgt::vec2(indicator->position_ * bar_.length_,0);
        else
            tmp = tgt::vec2(0,indicator->position_ * bar_.length_);
        return bar_.origin_ + tgt::ivec2(static_cast<int>(tmp.x), static_cast<int>(tmp.y));

    }

    tgt::ivec2 TouchTableTwoIndicatorSlider::getIndicatorLL(Indicator* indicator) const {

        return tgt::ivec2(static_cast<int>(indicator->bounds_.getLLF().x),
            static_cast<int>(indicator->bounds_.getLLF().y));
    }

    tgt::ivec2 TouchTableTwoIndicatorSlider::getIndicatorUR(Indicator* indicator) const {
        return tgt::ivec2(static_cast<int>(indicator->bounds_.getURB().x),
            static_cast<int>(indicator->bounds_.getURB().y));
    }

    void TouchTableTwoIndicatorSlider::releaseIndicator(int id){

        if(leftIndicator_.associatedTouchPointID_== id){
            leftIndicator_.associatedTouchPointID_ = -1;
        }else if(rightIndicator_.associatedTouchPointID_== id){
            rightIndicator_.associatedTouchPointID_ = -1;
        }else{
            //passed TP's ID does not match with that of an indicator
        }

    }

    tgt::vec2 TouchTableTwoIndicatorSlider::getMinMax() const {
       return tgt::vec2(leftIndicator_.position_, rightIndicator_.position_);
    }

    void TouchTableTwoIndicatorSlider::updateIndicatorPosition(tgt::vec2 pos, int id){

        Indicator* activeIndicator;
        Indicator* passivIndicator;
        tgt::ivec2 newPos = tgt::ivec2(pos);

        if(leftIndicator_.associatedTouchPointID_ == id){
            activeIndicator = &leftIndicator_;
            passivIndicator = &rightIndicator_;
        }else if (rightIndicator_.associatedTouchPointID_ == id){
            activeIndicator = &rightIndicator_;
            passivIndicator = &leftIndicator_;
        }else
            //passed TP's ID does not match with that of an indicator
            return;

        if(bar_.orientation_ == HORIZONTAL)
            newPos.y = bar_.origin_.y;
        else
            newPos.x = bar_.origin_.x;


        tgt::vec2 tmp = (tgt::vec2(newPos - bar_.origin_))/static_cast<float>(bar_.length_);
        activeIndicator->position_ = (float)(tmp.x + tmp.y);
        activeIndicator->position_ = tgt::clamp(activeIndicator->position_, 0.f, 1.f);

        if(activeIndicator->isLeft_ && (activeIndicator->position_ >= passivIndicator->position_)){
            passivIndicator->position_ = activeIndicator->position_;
        }else if (!activeIndicator->isLeft_ && (activeIndicator->position_ <= passivIndicator->position_)){
            passivIndicator->position_ = activeIndicator->position_;
        }

        updateIndicatorBounds(activeIndicator);
        updateIndicatorBounds(passivIndicator);

    }

    tgt::vec2 TouchTableTwoIndicatorSlider::getIndicatorDimensions() const {

        return tgt::ivec2(leftIndicator_.position_, rightIndicator_.position_);

    }

    std::deque<tgt::ivec2> TouchTableTwoIndicatorSlider::getIndicatorPositions() {

        std::deque<tgt::ivec2> tmp;
        tmp.push_back(getIndicatorLL(&leftIndicator_));
        tmp.push_back(getIndicatorUR(&leftIndicator_));
        tmp.push_back(getIndicatorLL(&rightIndicator_));
        tmp.push_back(getIndicatorUR(&rightIndicator_));
        return tmp;
    }

    void TouchTableTwoIndicatorSlider::setBarOrigin(tgt::ivec2 origin){
        bar_.origin_ = origin;
        updateIndicatorBounds(&leftIndicator_);
        updateIndicatorBounds(&rightIndicator_);
    }

    void TouchTableTwoIndicatorSlider::setBarLength(int length){
        bar_.length_ = length;
        updateIndicatorBounds(&leftIndicator_);
        updateIndicatorBounds(&rightIndicator_);
    }

    void TouchTableTwoIndicatorSlider::setIndicatorHeight(int height){

        leftIndicator_.height_ = height;
        rightIndicator_.height_ = height;
        updateIndicatorBounds(&leftIndicator_);
        updateIndicatorBounds(&rightIndicator_);
    }

    void TouchTableTwoIndicatorSlider::setIndicatorWidth(int width){

        leftIndicator_.width_ = width;
        rightIndicator_.width_ = width;
        updateIndicatorBounds(&leftIndicator_);
        updateIndicatorBounds(&rightIndicator_);
    }

    int TouchTableTwoIndicatorSlider::getIndicatorHeight() const {
        return leftIndicator_.height_;
    }

    int TouchTableTwoIndicatorSlider::getIndicatorWidth() const {
        return leftIndicator_.width_;
    }

/*************************************************************************************/
}
