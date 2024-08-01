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

#include "voreen/qt/widgets/property/transfunc/utils/quadrupleslider.h"

#include "tgt/tgt_math.h"
#include "voreen/core/utils/stringutils.h"
#include "tgt/logmanager.h"

#include <math.h>

#include <QColor>
#include <QMouseEvent>
#include <QPainter>
#include <QApplication>
#include <QToolTip>

namespace voreen {

QuadrupleSlider::QuadrupleSlider(QWidget* parent)
    : QWidget(parent)
{
    leftOffset_ = 0;
    rightOffset_ = 0;

    minValueMain_ = 0.f;
    maxValueMain_ = 1.f;
    minValueExtra_ = 0.f;
    maxValueExtra_ = 1.f;

    leftMainSliderActive_ = false;
    rightMainSliderActive_ = false;
    leftExtraSliderActive_ = false;
    rightExtraSliderActive_ = false;

    useExtraSlider_ = true;

    minimalAllowedSliderDistance_ = 0.0001f;
    sliderMainWidth_ = 9; //6
    sliderExtraWidth_ = 6;
    minRWValue_ = 0.0f;
    maxRWValue_ = 1.0f;
    unit_ = "";
    showToolTip_ = false;
    setFixedHeight(16); //16

    setToolTip(QString::fromStdString(generateToolTipText(false, false, false, false)));
}

void QuadrupleSlider::setOffsets(int left, int right) {
    leftOffset_ = left;
    rightOffset_ = right;
}

void QuadrupleSlider::paintEvent(QPaintEvent* event) {
    event->accept();
    int leftMainMarker =  tgt::iround(minValueMain_ * (width()-leftOffset_-rightOffset_) + leftOffset_);
    int rightMainMarker = tgt::iround(maxValueMain_ * (width()-leftOffset_-rightOffset_) + leftOffset_);
    int leftExtraMarker =  tgt::iround(minValueExtra_ * (width()-leftOffset_-rightOffset_) + leftOffset_);
    int rightExtraMarker = tgt::iround(maxValueExtra_ * (width()-leftOffset_-rightOffset_) + leftOffset_);
    QPoint leftMainSlider[4] = {
        QPoint(leftMainMarker, height()),
        QPoint(leftMainMarker + sliderMainWidth_, height()),
        QPoint(leftMainMarker + sliderMainWidth_, static_cast<int>(0.25f * height())),
        QPoint(leftMainMarker, 0)
    };
    QPoint leftExtraSlider[4] = {
        QPoint(leftExtraMarker, height()),
        QPoint(leftExtraMarker + sliderExtraWidth_, height()),
        QPoint(leftExtraMarker + sliderExtraWidth_, static_cast<int>(0.6f * height())),
        QPoint(leftExtraMarker, static_cast<int>(0.4f * height()))
    };
    QPoint rightMainSlider[4] = {
        QPoint(rightMainMarker - sliderMainWidth_, static_cast<int>(0.25f * height())),
        QPoint(rightMainMarker - sliderMainWidth_, height()),
        QPoint(rightMainMarker, height()),
        QPoint(rightMainMarker, 0)
    };
    QPoint rightExtraSlider[4] = {
        QPoint(rightExtraMarker, height()),
        QPoint(rightExtraMarker - sliderExtraWidth_, height()),
        QPoint(rightExtraMarker - sliderExtraWidth_, static_cast<int>(0.6f * height())),
        QPoint(rightExtraMarker, static_cast<int>(0.4f * height()))
    };
    QPoint centerQuad[4] = {
        QPoint(rightMainMarker - sliderMainWidth_, static_cast<int>(0.25f * height())),
        QPoint(rightMainMarker - sliderMainWidth_, height()),
        QPoint(leftMainMarker + sliderMainWidth_, height()),
        QPoint(leftMainMarker + sliderMainWidth_, static_cast<int>(0.25f * height())),
    };

    QPalette pal = QApplication::palette();

    QColor sliderMainColor = pal.color(QPalette::Base);
    QColor sliderExtraColor = pal.color(QPalette::Button);
    QColor sliderDarkColor = pal.color(QPalette::Mid);
    QColor lineColor = pal.color(QPalette::Dark);

    QPainter paint(this);
    paint.setRenderHint(QPainter::Antialiasing);

    //draw horizontal line
    paint.setPen(lineColor);
    paint.drawLine(leftOffset_, height()/2, width()-rightOffset_, height()/2);

    //draw center
    paint.setBrush(sliderDarkColor);
    paint.setPen(lineColor);
    paint.save();
    paint.drawConvexPolygon(centerQuad, 4);
    paint.restore();

    paint.setBrush(sliderMainColor);
    paint.setPen(sliderDarkColor);

    //draw left main marker
    paint.save();
    paint.drawConvexPolygon(leftMainSlider, 4);
    paint.restore();
    //draw right main marker
    paint.save();
    paint.drawConvexPolygon(rightMainSlider, 4);
    paint.restore();

    paint.setBrush(sliderExtraColor);
    paint.setPen(sliderDarkColor);

    if(useExtraSlider_) {
        //draw left extra marker
        paint.save();
        paint.drawConvexPolygon(leftExtraSlider, 4);
        paint.restore();
        //draw right extra marker
        paint.save();
        paint.drawConvexPolygon(rightExtraSlider, 4);
        paint.restore();
    }
}

void QuadrupleSlider::mousePressEvent(QMouseEvent* e) {
    e->accept();
    globalMousePos_ = e->globalPos();
    //calculate which marker is nearest to mouse position
    normalizedMousePos_ = static_cast<float>(e->pos().x()-leftOffset_) / static_cast<float>(width()-leftOffset_-rightOffset_);
    int leftMainMarker =  tgt::iround(minValueMain_ * (width()-leftOffset_-rightOffset_) + leftOffset_);
    int rightMainMarker = tgt::iround(maxValueMain_ * (width()-leftOffset_-rightOffset_) + leftOffset_);
    int leftExtraMarker =  tgt::iround(minValueExtra_ * (width()-leftOffset_-rightOffset_) + leftOffset_);
    int rightExtraMarker = tgt::iround(maxValueExtra_ * (width()-leftOffset_-rightOffset_) + leftOffset_);

    mV1_ = minValueMain_;
    mV2_ = maxValueMain_;

    if (e->button() == Qt::LeftButton) {
        if( useExtraSlider_ && ( (sliderExtraWidth_)*(e->pos().y()-(0.4f*height())) - (0.2f*height())*(e->pos().x() - leftExtraMarker) >= 0 ) &&
          (e->pos().x() <= (leftExtraMarker + sliderExtraWidth_)) && (e->pos().x() >= leftExtraMarker) && (e->pos().y() <= height())) {
            leftMainSliderActive_ = false;
            rightMainSliderActive_ = false;
            leftExtraSliderActive_ = true;
            rightExtraSliderActive_ = false;
        }
        else if( useExtraSlider_ && ( (sliderExtraWidth_)*(e->pos().y()-(0.6f*height())) - (-0.2f*height())*(e->pos().x() - (rightExtraMarker - sliderExtraWidth_)) >= 0 ) &&
          (e->pos().x() <= rightExtraMarker) && (e->pos().x() >= (rightExtraMarker -  + sliderExtraWidth_)) && (e->pos().y() <= height())) {
            leftMainSliderActive_ = false;
            rightMainSliderActive_ = false;
            leftExtraSliderActive_ = false;
            rightExtraSliderActive_ = true;
        }
        else if (e->pos().x() < (leftMainMarker + sliderMainWidth_)) {
            leftMainSliderActive_ = true;
            rightMainSliderActive_ = false;
            leftExtraSliderActive_ = false;
            rightExtraSliderActive_ = false;
        }
        else if(e->pos().x() > (rightMainMarker - sliderMainWidth_)) {
            leftMainSliderActive_ = false;
            rightMainSliderActive_ = true;
            leftExtraSliderActive_ = false;
            rightExtraSliderActive_ = false;
        }
        else if((leftMainMarker < e->pos().x()) && (rightMainMarker > e->pos().x())){
            leftMainSliderActive_ = true;
            rightMainSliderActive_ = true;
            leftExtraSliderActive_ = false;
            rightExtraSliderActive_ = false;
        }
    }
    else if (e->button() == Qt::RightButton) {
        leftMainSliderActive_ = true;
        rightMainSliderActive_ = true;
    }

    moveSlider(normalizedMousePos_);
    emit toggleInteractionMode(true);
}

void QuadrupleSlider::mouseMoveEvent(QMouseEvent* e){
    e->accept();
    globalMousePos_ = e->globalPos();
    float normalizedMousePosTmp = static_cast<float>((e->pos()).x()-leftOffset_) / static_cast<float>(width()-leftOffset_-rightOffset_);
    if (normalizedMousePosTmp > 1.f)
        normalizedMousePosTmp = 1.f;
    else if (normalizedMousePosTmp < 0.f)
        normalizedMousePosTmp = 0.f;
    moveSlider(normalizedMousePosTmp);

    emit mainValuesChanged(minValueMain_, maxValueMain_);
    emit extraValuesChanged(minValueExtra_, maxValueExtra_);
}

void QuadrupleSlider::showToolTip(std::string text) {
    if(!showToolTip_)
        return;

    //QToolTip::showText(globalMousePos_, "");
    QToolTip::showText(globalMousePos_, QString::fromStdString(text));
}

void QuadrupleSlider::mouseReleaseEvent(QMouseEvent* e) {
    leftMainSliderActive_ = false;
    rightMainSliderActive_ = false;
    leftExtraSliderActive_ = false;
    rightExtraSliderActive_ = false;
    e->accept();
    emit toggleInteractionMode(false);
}

std::string QuadrupleSlider::generateToolTipText(bool minMainBold, bool maxMainBold, bool minExtraBold, bool maxExtraBold) {
    std::string text = "<p style='white-space:pre'>";

    if(minMainBold)
        text += "<b>";
    text += ftos(getMappedValue(getMainMinValue()));
    if(minMainBold)
        text += "</b>";
    text += " -> ";
    if(maxMainBold)
        text += "<b>";
    text += ftos(getMappedValue(getMainMaxValue()));
    if(maxMainBold)
        text += "</b>";
    if(unit_ != "")
        text += " [" + unit_ + "]";

    if(useExtraSlider_) {
        text += "<br>";

        if(minExtraBold)
            text += "<b>";
        text += ftos(getMappedValue(getExtraMinValue()));
        if(minExtraBold)
            text += "</b>";
        text += " -> ";
        if(maxExtraBold)
            text += "<b>";
        text += ftos(getMappedValue(getExtraMaxValue()));
        if(maxExtraBold)
            text += "</b>";
        if(unit_ != "")
            text += " [" + unit_ + "]";
    }

    text += "<br>w: " + ftos(getMappedValue(getMainMaxValue()) - getMappedValue(getMainMinValue()));
    text += "<br>l: " + ftos( (getMappedValue(getMainMaxValue()) + getMappedValue(getMainMinValue())) * 0.5f);

    return text;
}

void QuadrupleSlider::moveSlider(float mousePos) {
    if(leftExtraSliderActive_) {
        setExtraMinValue(mousePos);
        showToolTip(generateToolTipText(false, false, true, false));
    }
    else if(rightExtraSliderActive_) {
        setExtraMaxValue(mousePos);
        showToolTip(generateToolTipText(false, false, false, true));
    }
    else if (leftMainSliderActive_ && !rightMainSliderActive_) {
        setMainMinValue(mousePos);
        showToolTip(generateToolTipText(true, false, false, false));
    }
    else if (rightMainSliderActive_ && !leftMainSliderActive_) {
        setMainMaxValue(mousePos);
        showToolTip(generateToolTipText(false, true, false, false));
    }
    else if (rightMainSliderActive_ && leftMainSliderActive_) {
        float mouseDiff = normalizedMousePos_ - mousePos;
        setMainMinValue(mV1_ - mouseDiff);
        setMainMaxValue(mV2_ - mouseDiff);
        showToolTip(generateToolTipText(true, true, false, false));
    }
}

void QuadrupleSlider::setMainMinValue(float val) {
    if (val == minValueMain_)
        return;

    if (val < 0.f)
        val = 0.f;
    if (val + minimalAllowedSliderDistance_ < maxValueMain_)
        minValueMain_ = val;
    else {
        maxValueMain_ = val + minimalAllowedSliderDistance_;
        if (maxValueMain_ > 1.f) {
            maxValueMain_ = 1.f;
            minValueMain_ = 1.f - minimalAllowedSliderDistance_;
        }
        else
            minValueMain_ = val;
    }
    emit mainValuesChanged(minValueMain_, maxValueMain_);
    setExtraMinValue(std::max(minValueMain_,minValueExtra_));
    setToolTip(QString::fromStdString(generateToolTipText(false, false, false, false)));
    update();
}

void QuadrupleSlider::setMainMaxValue(float val) {
    if (val == maxValueMain_)
        return;

    if (val > 1.f)
        val = 1.f;
    if (minValueMain_ + minimalAllowedSliderDistance_ < val)
        maxValueMain_ = val;
    else {
        minValueMain_ = val - minimalAllowedSliderDistance_;
        if (minValueMain_ < 0.f) {
            minValueMain_ = 0.f;
            maxValueMain_ = minimalAllowedSliderDistance_;
        }
        else
            maxValueMain_ = val;
    }
    emit mainValuesChanged(minValueMain_, maxValueMain_);
    setExtraMaxValue(std::min(maxValueMain_,maxValueExtra_));
    setToolTip(QString::fromStdString(generateToolTipText(false, false, false, false)));
    update();
}

void QuadrupleSlider::setExtraMinValue(float val) {
    if (val == minValueExtra_)
        return;

    if (val < minValueMain_)
        val = minValueMain_;
    if (val + minimalAllowedSliderDistance_ < maxValueExtra_)
        minValueExtra_ = val;
    else {
        maxValueExtra_ = val + minimalAllowedSliderDistance_;
        if (maxValueExtra_ > maxValueMain_) {
            maxValueExtra_ = maxValueMain_;
            minValueExtra_ = maxValueMain_ - minimalAllowedSliderDistance_;
        }
        else
            minValueExtra_ = val;
    }
    setToolTip(QString::fromStdString(generateToolTipText(false, false, false, false)));
    update();
    emit extraValuesChanged(minValueExtra_, maxValueExtra_);
}

void QuadrupleSlider::setExtraMaxValue(float val) {
    if (val == maxValueExtra_)
        return;

    if (val > maxValueMain_)
        val = maxValueMain_;
    if (minValueExtra_ + minimalAllowedSliderDistance_ < val)
        maxValueExtra_ = val;
    else {
        minValueExtra_ = val - minimalAllowedSliderDistance_;
        if (minValueExtra_ < minValueMain_) {
            minValueExtra_ = minValueMain_;
            maxValueExtra_ = minimalAllowedSliderDistance_;
        }
        else
            maxValueExtra_ = val;
    }
    setToolTip(QString::fromStdString(generateToolTipText(false, false, false, false)));
    update();
    emit extraValuesChanged(minValueExtra_, maxValueExtra_);
}

void QuadrupleSlider::setSliderMainWidth(int sliderW) {
    sliderMainWidth_ = sliderW;
}

void QuadrupleSlider::setSliderExtraWidth(int sliderW) {
    sliderExtraWidth_ = sliderW;
}

void QuadrupleSlider::setMainValues(float val1, float val2) {
    if (val1 < val2) {
        setMainMinValue(val1);
        setMainMaxValue(val2);
    }
    else {
        setMainMinValue(val2);
        setMainMaxValue(val1);
    }
}

void QuadrupleSlider::setExtraValues(float val1, float val2) {
    if (val1 < val2) {
        setExtraMinValue(val1);
        setExtraMaxValue(val2);
    }
    else {
        setExtraMinValue(val2);
        setExtraMaxValue(val1);
    }
}

float QuadrupleSlider::getMainMinValue() {
    return minValueMain_;
}

float QuadrupleSlider::getMainMaxValue() {
    return maxValueMain_;
}

float QuadrupleSlider::getExtraMinValue() {
    return minValueExtra_;
}

float QuadrupleSlider::getExtraMaxValue() {
    return maxValueExtra_;
}

void QuadrupleSlider::setUseExtraSlider(bool b) {
    useExtraSlider_ = b;
    update();
}

void QuadrupleSlider::setMapping(float min, float max) {
    minRWValue_ = min;
    maxRWValue_ = max;
    setToolTip(QString::fromStdString(generateToolTipText(false, false, false, false)));
}

void QuadrupleSlider::setUnit(std::string unit) {
    unit_ = unit;
    setToolTip(QString::fromStdString(generateToolTipText(false, false, false, false)));
}

void QuadrupleSlider::showToolTip(bool stt) {
    showToolTip_ = stt;
}

float QuadrupleSlider::getMappedValue(float norm) {
    return minRWValue_ + (norm * (maxRWValue_ - minRWValue_));
}

void QuadrupleSlider::setMinimalAllowedSliderDistance(float dist) {
    minimalAllowedSliderDistance_ = dist;
}

} // namespace voreen
