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

#include "voreen/qt/widgets/property/transfunc/utils/transfuncpropertywidgetpainterslider.h"

#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

//**********************************************************************************
// slider base functions
//**********************************************************************************
TransFuncPropertyWidgetPainterBaseSlider::TransFuncPropertyWidgetPainterBaseSlider(size_t associatedDimension, SliderOrientation orientation, int sliderSize)
    : canvasSizeInPixel_(100,40)
    , subCanvasOffsetInPixel_(0,0)
    , subCanvasSizeInPixel_(100,40)
    , valueRange_(0.f,1.f)
    , orientation_(orientation)
    , associatedDimension_(associatedDimension)
    , sliderSize_(sliderSize)
{
}

void TransFuncPropertyWidgetPainterBaseSlider::setCanvasSize(tgt::ivec2 size) {
    if(canvasSizeInPixel_ != size) {
        canvasSizeInPixel_ = size;
        updateSliderPosition();
    }
}

void TransFuncPropertyWidgetPainterBaseSlider::setSubCanvasOffset(tgt::ivec2 offset) {
    if(subCanvasOffsetInPixel_ != offset) {
        subCanvasOffsetInPixel_ = offset;
        updateSliderPosition();
    }
}

void TransFuncPropertyWidgetPainterBaseSlider::setSubCanvasSize(tgt::ivec2 size) {
    if(subCanvasSizeInPixel_ != size) {
        subCanvasSizeInPixel_ = size;
        updateSliderPosition();
    }
}

void TransFuncPropertyWidgetPainterBaseSlider::setValueRange(tgt::vec2 range) {
    if(valueRange_ != range) {
        valueRange_ = range;
        updateSliderPosition();
    }
}

//**********************************************************************************
// slider gamma functions
//**********************************************************************************
TransFuncPropertyWidgetPainterGammaSlider::TransFuncPropertyWidgetPainterGammaSlider(size_t associatedDimension, SliderOrientation orientation, int sliderSize)
    : TransFuncPropertyWidgetPainterBaseSlider(associatedDimension, orientation, sliderSize)
    , position_(0)
    , gamma_(1.f)
{
    valueRange_ = tgt::vec2(0.1f,5.f);
}

void TransFuncPropertyWidgetPainterGammaSlider::paint() {
    //define variables
    tgt::vec3 point1, point2, point3;
    //switch orientation
    switch(orientation_) {
    case SO_TOP:
        point1 = tgt::vec3(position_-sliderSize_,canvasSizeInPixel_.y-1,-0.5f);
        point2 = tgt::vec3(position_+sliderSize_,canvasSizeInPixel_.y-1,-0.5f);
        point3 = tgt::vec3(position_,canvasSizeInPixel_.y-1 - sliderSize_,-0.5f);
        break;
    case SO_BOTTOM:
        point1 = tgt::vec3(position_-sliderSize_,0.0f,-0.5f);
        point2 = tgt::vec3(position_+sliderSize_,0.0f,-0.5f);
        point3 = tgt::vec3(position_,sliderSize_,-0.5f);
        break;
    case SO_LEFT:
        point1 = tgt::vec3(0.0f,position_-sliderSize_,-0.5f);
        point2 = tgt::vec3(0.0f,position_+sliderSize_,-0.5f);
        point3 = tgt::vec3(sliderSize_,position_,-0.5f);
        break;
    case SO_RIGHT:
        point1 = tgt::vec3(canvasSizeInPixel_.x-1,position_-sliderSize_,-0.5f);
        point2 = tgt::vec3(canvasSizeInPixel_.x-1,position_+sliderSize_,-0.5f);
        point3 = tgt::vec3(canvasSizeInPixel_.y-1-sliderSize_,position_,-0.5f);
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
    //draw black arrow
    IMode.begin(tgt::ImmediateMode::TRIANGLES);
        IMode.color(tgt::vec4(0.f,0.f,0.f,1.f));
        IMode.vertex(point1);
        IMode.vertex(point2);
        IMode.vertex(point3);
    IMode.end();
    //draw white border
    glLineWidth(2.f);
    IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.color(tgt::vec4(1.f,1.f,1.f,1.f));
        IMode.vertex(point1);
        IMode.vertex(point2);
        IMode.vertex(point3);
    IMode.end();
    glLineWidth(1.f);
}

bool TransFuncPropertyWidgetPainterGammaSlider::isHit(tgt::ivec2 pos) {
    int toleranceInPixel = 5;

    switch(orientation_) {
    case SO_TOP:
        if((std::abs(pos.x-position_) <= (sliderSize_+toleranceInPixel)) && (std::max(0,pos.y) <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_BOTTOM:
        if((std::abs(pos.x-position_) <= (sliderSize_+toleranceInPixel)) && (std::abs(canvasSizeInPixel_.y-1 - pos.y) <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_LEFT:
        if((std::abs((canvasSizeInPixel_.y-1 -pos.y) - position_) <= (sliderSize_+toleranceInPixel)) && (std::max(0,pos.x) <= (sliderSize_+toleranceInPixel)))
            return true;
    case SO_RIGHT:
        if((std::abs((canvasSizeInPixel_.y-1 -pos.y) - position_) <= (sliderSize_+toleranceInPixel)) && (std::abs(canvasSizeInPixel_.x-1 - pos.x) <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
    //return default false
    return false;
}

void TransFuncPropertyWidgetPainterGammaSlider::setPosition(tgt::ivec2 pos) {
    //invert y and clamp position
    pos.y = canvasSizeInPixel_.y - 1 - pos.y;
    pos = tgt::clamp(pos,subCanvasOffsetInPixel_,subCanvasOffsetInPixel_+subCanvasSizeInPixel_-1);

    switch(orientation_) {
    case SO_TOP:
    case SO_BOTTOM:
        if(position_ != pos.x) {
            position_ = pos.x;
            updateStoredValue();
        }
        break;
    case SO_LEFT:
    case SO_RIGHT:
        if(position_ != pos.y) {
            position_ = pos.y;
            updateStoredValue();
        }
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
}

int TransFuncPropertyWidgetPainterGammaSlider::getPosition() {
    return position_;
}

void TransFuncPropertyWidgetPainterGammaSlider::setGammaValue(float gamma) {
    gamma = tgt::clamp(gamma,valueRange_.x,valueRange_.y);
    if(gamma_ != gamma) {
        gamma_ = gamma;
        updateSliderPosition();
    }
}

float TransFuncPropertyWidgetPainterGammaSlider::getGammaValue() {
    return gamma_;
}

void TransFuncPropertyWidgetPainterGammaSlider::updateSliderPosition() {
    float position = 1.f;
    int firstCanvasHalf, secondCanvasHalf;

    switch(orientation_) {
    case SO_TOP:
    case SO_BOTTOM:
        firstCanvasHalf = subCanvasSizeInPixel_.x/2;
        secondCanvasHalf = subCanvasSizeInPixel_.x-firstCanvasHalf;
        if(gamma_ < 1.f) {
            //gamma between min and 1
            position = static_cast<int>(((gamma_-valueRange_.x)*static_cast<float>(firstCanvasHalf))/(1.f-valueRange_.x))+subCanvasOffsetInPixel_.x;
        } else {
           //gamma between 1 and max
            position = static_cast<int>(((gamma_-1.f)*static_cast<float>(secondCanvasHalf-1))/(valueRange_.y-1.f))+subCanvasOffsetInPixel_.x+firstCanvasHalf;
        }
        break;
    case SO_LEFT:
    case SO_RIGHT:
        firstCanvasHalf = subCanvasSizeInPixel_.y/2;
        secondCanvasHalf = subCanvasSizeInPixel_.y-firstCanvasHalf;
        if(gamma_ < 1.f) {
            //gamma between min and 1
            position = static_cast<int>(((gamma_-valueRange_.x)*static_cast<float>(firstCanvasHalf))/(1.f-valueRange_.x))+subCanvasOffsetInPixel_.y;
        } else {
           //gamma between 1 and max
            position = static_cast<int>(((gamma_-1.f)*static_cast<float>(secondCanvasHalf-1))/(valueRange_.y-1.f))+subCanvasOffsetInPixel_.y+firstCanvasHalf;
        }
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
    //update position
    position_ = position;
}

void TransFuncPropertyWidgetPainterGammaSlider::updateStoredValue() {
    float gamma = 1.f;
    int firstCanvasHalf, secondCanvasHalf;

    switch(orientation_) {
    case SO_TOP:
    case SO_BOTTOM:
        firstCanvasHalf = subCanvasSizeInPixel_.x/2;
        secondCanvasHalf = subCanvasSizeInPixel_.x-firstCanvasHalf;
        if(position_ < subCanvasOffsetInPixel_.x+firstCanvasHalf) {
            //gamma between min and 1
            gamma = ((1.f - valueRange_.x)*static_cast<float>(position_-subCanvasOffsetInPixel_.x))/static_cast<float>(firstCanvasHalf)+valueRange_.x;
        } else {
           //gamma between 1 and max
            gamma = ((valueRange_.y - 1.f)*static_cast<float>(position_-firstCanvasHalf-subCanvasOffsetInPixel_.x))/static_cast<float>(secondCanvasHalf-1)+1.f;
        }
        break;
    case SO_LEFT:
    case SO_RIGHT:
        firstCanvasHalf = subCanvasSizeInPixel_.y/2;
        secondCanvasHalf = subCanvasSizeInPixel_.y-firstCanvasHalf;
        if(position_ < subCanvasOffsetInPixel_.y+firstCanvasHalf) {
            //gamma between min and 1
            gamma = ((1.f - valueRange_.x)*static_cast<float>(position_-subCanvasOffsetInPixel_.y))/static_cast<float>(firstCanvasHalf)+valueRange_.x;
        } else {
           //gamma between 1 and max
            gamma = ((valueRange_.y - 1.f)*static_cast<float>(position_-firstCanvasHalf-subCanvasOffsetInPixel_.y))/static_cast<float>(secondCanvasHalf-1)+1.f;
        }
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }

    //update gamma, if it has been changed
    if(gamma_ != gamma) {
        gamma_ = gamma;
        emit gammaChanged(gamma, associatedDimension_);
    }
}

//***********************************************************************************
// slider domain functions
//***********************************************************************************
TransFuncPropertyWidgetPainterDomainSlider::TransFuncPropertyWidgetPainterDomainSlider(size_t associatedDimension, SliderOrientation orientation, int sliderSize)
    : TransFuncPropertyWidgetPainterBaseSlider(associatedDimension, orientation, sliderSize)
    , leftPosition_(0)
    , rightPosition_(100)
    , domain_(0.f,1.f)
{
}

void TransFuncPropertyWidgetPainterDomainSlider::paint() {
   tgt::vec3 pointLeft1, pointLeft2, pointLeft3;
   tgt::vec3 pointRight1, pointRight2, pointRight3;

    switch(orientation_) {
    case SO_TOP:
        pointLeft1 = tgt::vec3(leftPosition_,canvasSizeInPixel_.y-1,-0.5f);
        pointLeft2 = tgt::vec3(leftPosition_,canvasSizeInPixel_.y-1-sliderSize_,-0.5f);
        pointLeft3 = tgt::vec3(leftPosition_+sliderSize_,canvasSizeInPixel_.y-1,-0.5f);
        pointRight1 = tgt::vec3(rightPosition_,canvasSizeInPixel_.y-1,-0.5f);
        pointRight2 = tgt::vec3(rightPosition_,canvasSizeInPixel_.y-1-sliderSize_,-0.5f);
        pointRight3 = tgt::vec3(rightPosition_-sliderSize_,canvasSizeInPixel_.y-1,-0.5f);
        break;
    case SO_BOTTOM:
        pointLeft1 = tgt::vec3(leftPosition_,0.0f,-0.5f);
        pointLeft2 = tgt::vec3(leftPosition_,sliderSize_,-0.5f);
        pointLeft3 = tgt::vec3(leftPosition_+sliderSize_,0.0f,-0.5f);
        pointRight1 = tgt::vec3(rightPosition_,0.0f,-0.5f);
        pointRight2 = tgt::vec3(rightPosition_,sliderSize_,-0.5f);
        pointRight3 = tgt::vec3(rightPosition_-sliderSize_,0.0f,-0.5f);
        break;
    case SO_LEFT:
        pointLeft1 = tgt::vec3(0.0f,leftPosition_,-0.5f);
        pointLeft2 = tgt::vec3(sliderSize_,leftPosition_,-0.5f);
        pointLeft3 = tgt::vec3(0.0f,leftPosition_+sliderSize_,-0.5f);
        pointRight1 = tgt::vec3(0.0f,rightPosition_,-0.5f);
        pointRight2 = tgt::vec3(sliderSize_,rightPosition_,-0.5f);
        pointRight3 = tgt::vec3(0.0f,rightPosition_-sliderSize_,-0.5f);
        break;
    case SO_RIGHT:
        pointLeft1 = tgt::vec3(canvasSizeInPixel_.x-1,leftPosition_,-0.5f);
        pointLeft2 = tgt::vec3(canvasSizeInPixel_.x-1-sliderSize_,leftPosition_,-0.5f);
        pointLeft3 = tgt::vec3(canvasSizeInPixel_.x-1,leftPosition_+sliderSize_,-0.5f);
        pointRight1 = tgt::vec3(canvasSizeInPixel_.x-1,rightPosition_,-0.5f);
        pointRight2 = tgt::vec3(canvasSizeInPixel_.x-1-sliderSize_,rightPosition_,-0.5f);
        pointRight3 = tgt::vec3(canvasSizeInPixel_.x-1,rightPosition_-sliderSize_,-0.5f);
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
    //LEFT
    //draw black arrow
    IMode.begin(tgt::ImmediateMode::TRIANGLES);
        IMode.color(tgt::vec4(0.f,0.f,0.f,1.f));
        IMode.vertex(pointLeft1); IMode.vertex(pointLeft2); IMode.vertex(pointLeft3);
    IMode.end();
    //draw white border
    glLineWidth(2.f);
    IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.color(tgt::vec4(1.f,1.f,1.f,1.f));
        IMode.vertex(pointLeft1);
        IMode.vertex(pointLeft2);
        IMode.vertex(pointLeft3);
    IMode.end();
    glLineWidth(1.f);

    //RIGHT
    //draw black arrow
    IMode.begin(tgt::ImmediateMode::TRIANGLES);
        IMode.color(tgt::vec4(0.f,0.f,0.f,1.f));
        IMode.vertex(pointRight1); IMode.vertex(pointRight2); IMode.vertex(pointRight3);
    IMode.end();
    //draw white border
    glLineWidth(2.f);
    IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.color(tgt::vec4(1.f,1.f,1.f,1.f));
        IMode.vertex(pointRight1);
        IMode.vertex(pointRight2);
        IMode.vertex(pointRight3);
    IMode.end();
    glLineWidth(1.f);
}

bool TransFuncPropertyWidgetPainterDomainSlider::isHit(tgt::ivec2 pos) {
    return isLeftHit(pos) || isRightHit(pos);
}

bool TransFuncPropertyWidgetPainterDomainSlider::isLeftHit(tgt::ivec2 pos) {
    int toleranceInPixel = 5;

    switch(orientation_) {
    case SO_TOP:
        if( (pos.x <= leftPosition_+sliderSize_+toleranceInPixel) && (pos.x >= (leftPosition_-toleranceInPixel)) &&
            (pos.y <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_BOTTOM:
        if( (pos.x <= leftPosition_+sliderSize_+toleranceInPixel) && (pos.x >= (leftPosition_-toleranceInPixel)) &&
            (std::abs(canvasSizeInPixel_.y - 1 -pos.y) <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_LEFT:
        if( (std::abs(canvasSizeInPixel_.y - 1 - pos.y) <= leftPosition_+sliderSize_+toleranceInPixel) &&
            (std::abs(canvasSizeInPixel_.y - 1 - pos.y) >= (leftPosition_-toleranceInPixel)) &&
            (pos.x <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_RIGHT:
        if( (std::abs(canvasSizeInPixel_.y - 1 - pos.y) <= leftPosition_+sliderSize_+toleranceInPixel) &&
            (std::abs(canvasSizeInPixel_.y - 1 - pos.y) >= (leftPosition_-toleranceInPixel)) &&
            (std::abs(canvasSizeInPixel_.x - 1 - pos.x) <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }

    //default return false
    return false;
}

bool TransFuncPropertyWidgetPainterDomainSlider::isRightHit(tgt::ivec2 pos) {
    int toleranceInPixel = 5;

    switch(orientation_) {
    case SO_TOP:
        if( (pos.x <= rightPosition_+toleranceInPixel) && (pos.x >= (rightPosition_-sliderSize_-toleranceInPixel)) &&
            (pos.y <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_BOTTOM:
        if( (pos.x <= rightPosition_+toleranceInPixel) && (pos.x >= (rightPosition_-sliderSize_-toleranceInPixel)) &&
            (std::abs(canvasSizeInPixel_.y - 1 -pos.y) <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_LEFT:
        if( (std::abs(canvasSizeInPixel_.y - 1 - pos.y) <= rightPosition_+toleranceInPixel) && (std::abs(canvasSizeInPixel_.y - 1 - pos.y) >= (rightPosition_-sliderSize_-toleranceInPixel)) &&
            (pos.x <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    case SO_RIGHT:
        if( (std::abs(canvasSizeInPixel_.y - 1 - pos.y) <= rightPosition_+toleranceInPixel) && (std::abs(canvasSizeInPixel_.y - 1 - pos.y) >= (rightPosition_-sliderSize_-toleranceInPixel)) &&
            (std::abs(canvasSizeInPixel_.x - 1 - pos.x) <= (sliderSize_+toleranceInPixel)))
            return true;
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }

    //default return false
    return false;
}

bool TransFuncPropertyWidgetPainterDomainSlider::isIntervalHit(tgt::ivec2 pos) {
    int toleranceInPixel = 5;

    switch(orientation_) {
    case SO_TOP:
        if( (pos.x <= rightPosition_) && (pos.x >= leftPosition_) &&
            (pos.y <= (sliderSize_+toleranceInPixel)) )
            return true;
        break;
    case SO_BOTTOM:
        if( (pos.x <= rightPosition_) && (pos.x >= leftPosition_) &&
            (std::abs(canvasSizeInPixel_.y - 1 - pos.y) <= (sliderSize_+toleranceInPixel)) )
            return true;
        break;
    case SO_LEFT:
        if( (canvasSizeInPixel_.y - 1 - pos.y <= rightPosition_) && (canvasSizeInPixel_.y - 1 - pos.y >= leftPosition_) &&
            (pos.x <= (sliderSize_+toleranceInPixel)) )
            return true;
        break;
    case SO_RIGHT:
        if( (canvasSizeInPixel_.y - 1 - pos.y <= rightPosition_) && (canvasSizeInPixel_.y - 1 - pos.y >= leftPosition_) &&
            (std::abs(canvasSizeInPixel_.x - 1 - pos.x) <= (sliderSize_+toleranceInPixel)) )
            return true;
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
    //default return false
    return false;
}

void TransFuncPropertyWidgetPainterDomainSlider::setLeftPosition(tgt::ivec2 pos) {
    //invert y and clamp
    pos.y = canvasSizeInPixel_.y - 1 - pos.y;
    pos = tgt::clamp(pos,subCanvasOffsetInPixel_,subCanvasOffsetInPixel_+subCanvasSizeInPixel_-1);

    switch(orientation_) {
    case SO_TOP:
    case SO_BOTTOM:
        if(leftPosition_ != pos.x) {
            leftPosition_ = pos.x;
            if(leftPosition_ > rightPosition_)
                rightPosition_ = leftPosition_;
            updateStoredValue();
        }
        break;
    case SO_LEFT:
    case SO_RIGHT:
        if(leftPosition_ != pos.y) {
            leftPosition_ = pos.y;
            if(leftPosition_ > rightPosition_)
                rightPosition_ = leftPosition_;
            updateStoredValue();
        }
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
}

void TransFuncPropertyWidgetPainterDomainSlider::setRightPosition(tgt::ivec2 pos) {
    //invert y and clamp
    pos.y = canvasSizeInPixel_.y - 1 - pos.y;
    pos = tgt::clamp(pos,subCanvasOffsetInPixel_,subCanvasOffsetInPixel_+subCanvasSizeInPixel_-1);

    switch(orientation_) {
    case SO_TOP:
    case SO_BOTTOM:
        if(rightPosition_ != pos.x) {
            rightPosition_ = pos.x;
            if(rightPosition_ < leftPosition_)
                leftPosition_ = rightPosition_;
            updateStoredValue();
        }
        break;
    case SO_LEFT:
    case SO_RIGHT:
        if(rightPosition_ != pos.y) {
            rightPosition_ = pos.y;
            if(rightPosition_ < leftPosition_)
                leftPosition_ = rightPosition_;
            updateStoredValue();
        }
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
}

tgt::ivec2 TransFuncPropertyWidgetPainterDomainSlider::getPosition() {
    return tgt::ivec2(leftPosition_,rightPosition_);
}

void TransFuncPropertyWidgetPainterDomainSlider::setDomainValue(tgt::vec2 domain) {
    domain.x = tgt::clamp(domain.x,valueRange_.x,valueRange_.y);
    domain.y = tgt::clamp(domain.y,valueRange_.x,valueRange_.y);
    if(domain_ != domain) {
        domain_ = domain;
        updateSliderPosition();
    }
}

tgt::vec2 TransFuncPropertyWidgetPainterDomainSlider::getDomainValue() {
    return domain_;
}

void TransFuncPropertyWidgetPainterDomainSlider::updateSliderPosition() {
    //initialize variables
    float range = valueRange_.y-valueRange_.x;

    //determine dimension
    int subCanvasOffset = 0, subCanvasSize = 0;
    switch(orientation_) {
    case SO_TOP:
    case SO_BOTTOM:
        subCanvasOffset = subCanvasOffsetInPixel_.x;
        subCanvasSize = subCanvasSizeInPixel_.x;
        break;
    case SO_LEFT:
    case SO_RIGHT:
        subCanvasOffset = subCanvasOffsetInPixel_.y;
        subCanvasSize = subCanvasSizeInPixel_.y;
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }
    //set new positions
    leftPosition_ = static_cast<int>((domain_.x-valueRange_.x)*static_cast<float>(subCanvasSize-1)/range)+subCanvasOffset;
    rightPosition_ = static_cast<int>((domain_.y-valueRange_.x)*static_cast<float>(subCanvasSize-1)/range)+subCanvasOffset;
}

void TransFuncPropertyWidgetPainterDomainSlider::updateStoredValue() {
    //initialize variables
    tgt::vec2 domain(0.f,1.f);
    float range = valueRange_.y-valueRange_.x;

    //determine dimension
    int subCanvasOffset, subCanvasSize;
    switch(orientation_) {
    case SO_TOP:
    case SO_BOTTOM:
        subCanvasOffset = subCanvasOffsetInPixel_.x;
        subCanvasSize = subCanvasSizeInPixel_.x;
        break;
    case SO_LEFT:
    case SO_RIGHT:
        subCanvasOffset = subCanvasOffsetInPixel_.y;
        subCanvasSize = subCanvasSizeInPixel_.y;
        break;
    default:
        tgtAssert(false, "unknwon orientation");
    }

    //calculate new domain
    domain.x = static_cast<float>(leftPosition_-subCanvasOffset)*range/static_cast<float>(subCanvasSize-1)+valueRange_.x;
    domain.y = static_cast<float>(rightPosition_-subCanvasOffset)*range/static_cast<float>(subCanvasSize-1)+valueRange_.x;

    //test, if domain is 0
    if(domain.x == domain.y) {
        if(domain.x - 0.0001f < valueRange_.x)
            domain.y += 0.0001f;
        else
            domain.x -= 0.0001f;
    }

    //update domain, if is was changed
    if(domain_ != domain) {
        domain_ = domain;
        emit domainChanged(domain, associatedDimension_);
    }
}


} // namespace voreen
