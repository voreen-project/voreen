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

#include "touchtabletimeline.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen {

    TouchTableTimeline::TouchTableTimeline(int duration, int stepSize)
        : duration_(duration)
        , stepSize_(stepSize)
        , timePerStep_(1)
        , timeOffset_(0)
        , currentPosition_(1)
        , currentTime_(0.f)
    { }

    int TouchTableTimeline::getStepSize(){
        return stepSize_;
    }

    void TouchTableTimeline::setStepSize(int stepSize){
        stepSize_ = stepSize;
    }

    int TouchTableTimeline::getDuration(){
        return duration_;
    }

    void TouchTableTimeline::setDuration(int duration){
        duration_  = duration;
    }

    int TouchTableTimeline::getTimeOffset() const{
        return timeOffset_;
    }

    void TouchTableTimeline::setTimeOffset(int offset){
        timeOffset_ = offset;
    }

    int TouchTableTimeline::getCurrentPos(){
        return currentPosition_;
    }

    void TouchTableTimeline::setCurrentPos(int pos){
        currentPosition_= timeOffset_ + pos;
        updateCurrentTime();
    }

    int TouchTableTimeline::getMainLine(){
        return mainLinePos_;
    }

    void TouchTableTimeline::setMainLine(int pos){
        mainLinePos_ = pos;
    }

    float TouchTableTimeline::getCurrentTime(){
        return currentTime_;
    }

    void TouchTableTimeline::setCurrentTime(float time){
        currentTime_ = time;
        currentPosition_ = static_cast<int>(time * stepSize_);
    }

    std::vector<KeyValue>& TouchTableTimeline::getKeyValues() {
        return keyValues_;
    }

    int TouchTableTimeline::getTimePerStep(){
        return timePerStep_;
    }

    void TouchTableTimeline::updateCurrentTime(){

        currentTime_ = static_cast<float>(currentPosition_) / static_cast<float>(stepSize_) * timePerStep_;

    }

    std::string TouchTableTimeline::getCurrentTimeAsString(){

        std::pair<int,int> timePair = std::make_pair((int) std::floor( static_cast<float>(currentTime_ / 60)), static_cast<int>(currentTime_) % 60);

        std::string sec =( timePair.second >=10) ? itos(timePair.second) : "0" + itos(timePair.second);
        std::string min =( timePair.first >=10) ? itos(timePair.first) : "0" + itos(timePair.first);
        std::string time=  min + ":" + sec;

        return time;
    }

    std::string TouchTableTimeline::getTimeForPosAsString(int pos){
        int currentSeconds= (int ) std::floor( static_cast<float>((timeOffset_ + pos) / stepSize_)) * timePerStep_;
        std::pair<int,int> time = std::make_pair((int) std::floor( static_cast<float>(currentSeconds / 60)), currentSeconds % 60);

        std::string min =( time.first >=10) ? itos( time.first) : "0" + itos( time.first);
        std::string sec =( time.second >=10) ? itos(time.second) : "0" + itos(time.second);
        return (min + ":" + sec);


    }
    float TouchTableTimeline::getTimeForPos(int pos){

        float time= static_cast<float>(pos) / static_cast<float>(stepSize_) * timePerStep_;

        return time;


    }


    bool TouchTableTimeline::addKeyValueAt(int pos, int height, int width){
        for(std::vector<KeyValue>::iterator keyValueIter = keyValues_.begin(); keyValueIter != keyValues_.end(); ++keyValueIter){
            if(std::abs(pos- keyValueIter->pos_) <= keyValueIter->width_){
                return false;
            }
        }
        KeyValue newKeyValue = KeyValue(pos, width, height, currentTime_);
        keyValues_.push_back(newKeyValue);

        return true;
    }

    KeyValue* TouchTableTimeline::keyValueHit(int xPos, int yPos){
        KeyValue* rtnKeyValue = 0;

        for(std::vector<KeyValue>::iterator keyValueIter = keyValues_.begin(); keyValueIter != keyValues_.end(); ++keyValueIter){
            if( (xPos+timeOffset_) >= keyValueIter->pos_- keyValueIter->width_/2 && (xPos+timeOffset_) <= keyValueIter->pos_+ keyValueIter->width_/2
                && yPos >= mainLinePos_ - keyValueIter->height_/2 && yPos <= mainLinePos_+ keyValueIter->height_/2){
                keyValueIter->isSelected_ = true;
                rtnKeyValue =  &(*keyValueIter);
            }else{
                keyValueIter->isSelected_ = false;
            }
        }
        return rtnKeyValue;
    }

    KeyValue* TouchTableTimeline::getSelectedKeyValue(){
        for(std::vector<KeyValue>::iterator keyValueIter = keyValues_.begin(); keyValueIter != keyValues_.end(); ++keyValueIter){
            if(keyValueIter->isSelected_)
                return &(*keyValueIter);
        }

        return 0;
    }

    void TouchTableTimeline::eraseSelectedKeyValue(){
        for(std::vector<KeyValue>::iterator keyValueIter = keyValues_.begin(); keyValueIter != keyValues_.end(); ++keyValueIter){
            if(keyValueIter->isSelected_){
                keyValues_.erase(keyValueIter);
                return;
            }
        }
    }

    bool TouchTableTimeline::anyKeyValueSelected(){
        for(std::vector<KeyValue>::iterator keyValueIter = keyValues_.begin(); keyValueIter != keyValues_.end(); ++keyValueIter){
            if(keyValueIter->isSelected_)
                return true;
        }
        return false;
    }

    bool TouchTableTimeline::validateKeyValuePos(int pos, int lastValidPos){
        for(std::vector<KeyValue>::iterator keyValueIter = keyValues_.begin(); keyValueIter != keyValues_.end(); ++keyValueIter){
            if(lastValidPos == keyValueIter->lastValidPos_) continue;

            if(std::abs(pos-keyValueIter->pos_) <= keyValueIter->width_)
                return false;
        }

        return true;
    }


}
