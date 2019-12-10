/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/transfunc/utils/doubleslider.h"

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

DoubleSlider::DoubleSlider(QWidget* parent, bool vertical)
    : QWidget(parent)
    , minValue_(0.f), maxValue_(1.f)
    , minimalAllowedSliderDistance_(0.0001f), sliderWidthHeight_(6)
    , minRWValue_(0.0f), maxRWValue_(1.0f)
    , unit_(""), showToolTip_(false)
    , leftTopOffset_(0), rightBottomOffset_(0)
    , leftTopSliderActive_(false), rightBottomSliderActive_(false)
    , vertical_(vertical)
{
    if(vertical_)
        setFixedWidth(16);
    else
        setFixedHeight(16);

    setToolTip(QString::fromStdString(generateToolTipText(false, false)));
}

void DoubleSlider::setOffsets(int leftTop, int rightBottom) {
    leftTopOffset_ = leftTop;
    rightBottomOffset_ = rightBottom;
}

void DoubleSlider::paintEvent(QPaintEvent* event) {
    event->accept();
    int oldWidth = (vertical_ ? height() : width());
    int leftTopMarker =  tgt::iround(minValue_ * (oldWidth-leftTopOffset_-rightBottomOffset_) + leftTopOffset_);
    int rightBottomMarker = tgt::iround(maxValue_ * (oldWidth-leftTopOffset_-rightBottomOffset_) + leftTopOffset_);
    QPoint leftTopSlider[4], rightBottomSlider[4], centerQuad[4];
    if(vertical_) {
        leftTopSlider[0] = QPoint(0,leftTopMarker);
        leftTopSlider[1] = QPoint(width(),leftTopMarker);
        leftTopSlider[2] = QPoint(static_cast<int>(0.6f * width()) ,leftTopMarker + sliderWidthHeight_);
        leftTopSlider[3] = QPoint(0,leftTopMarker + sliderWidthHeight_);
        rightBottomSlider[0]   = QPoint(0, rightBottomMarker);
        rightBottomSlider[1]   = QPoint(width(), rightBottomMarker);
        rightBottomSlider[2]   = QPoint(static_cast<int>(0.6f * width()), rightBottomMarker - sliderWidthHeight_);
        rightBottomSlider[3]   = QPoint(0,rightBottomMarker - sliderWidthHeight_);
        centerQuad[0]       = QPoint(0,rightBottomMarker - sliderWidthHeight_);
        centerQuad[1]       = QPoint(static_cast<int>(0.6f * width()),rightBottomMarker - sliderWidthHeight_);
        centerQuad[2]       = QPoint(static_cast<int>(0.6f * width()),leftTopMarker + sliderWidthHeight_);
        centerQuad[3]       = QPoint(0,leftTopMarker + sliderWidthHeight_);
    } else {
        leftTopSlider[0] = QPoint(leftTopMarker, height());
        leftTopSlider[1] = QPoint(leftTopMarker + sliderWidthHeight_, height());
        leftTopSlider[2] = QPoint(leftTopMarker + sliderWidthHeight_, static_cast<int>(0.4f * height()));
        leftTopSlider[3] = QPoint(leftTopMarker, 0);
        rightBottomSlider[0]   = QPoint(rightBottomMarker - sliderWidthHeight_, static_cast<int>(0.4f * height()));
        rightBottomSlider[1]   = QPoint(rightBottomMarker - sliderWidthHeight_, height());
        rightBottomSlider[2]   = QPoint(rightBottomMarker, height());
        rightBottomSlider[3]   = QPoint(rightBottomMarker, 0);
        centerQuad[0]       = QPoint(rightBottomMarker - sliderWidthHeight_, static_cast<int>(0.4f * height()));
        centerQuad[1]       = QPoint(rightBottomMarker - sliderWidthHeight_, height());
        centerQuad[2]       = QPoint(leftTopMarker + sliderWidthHeight_, height());
        centerQuad[3]       = QPoint(leftTopMarker + sliderWidthHeight_, static_cast<int>(0.4f * height()));
    }

    QPalette pal = QApplication::palette();

    QColor sliderColor = pal.color(QPalette::Base);
    QColor sliderDarkColor = pal.color(QPalette::Mid);
    QColor lineColor = pal.color(QPalette::Dark);

    QPainter paint(this);
    paint.setRenderHint(QPainter::Antialiasing);

    //draw horizontal line
    paint.setPen(lineColor);
    if(vertical_)
        paint.drawLine(width()/2, leftTopOffset_, width()/2, height()-rightBottomOffset_);
    else
        paint.drawLine(leftTopOffset_, height()/2, width()-rightBottomOffset_, height()/2);

    //draw center
    paint.setBrush(sliderDarkColor);
    paint.setPen(lineColor);
    paint.save();
    paint.drawConvexPolygon(centerQuad, 4);
    paint.restore();

    paint.setBrush(sliderColor);
    paint.setPen(sliderDarkColor);

    //draw left marker
    paint.save();
    paint.drawConvexPolygon(leftTopSlider, 4);
    paint.restore();

    //draw right marker
    paint.save();
    paint.drawConvexPolygon(rightBottomSlider, 4);
    paint.restore();
}

void DoubleSlider::mousePressEvent(QMouseEvent* e) {
    e->accept();
    globalMousePos_ = e->globalPos();
    //calculate which marker is nearest to mouse position
    int leftTopMarker, rightBottomMarker;
    if(vertical_) {
        normalizedMousePos_ = static_cast<float>(e->pos().y()-leftTopOffset_) / static_cast<float>(height()-leftTopOffset_-rightBottomOffset_);
        leftTopMarker =  tgt::iround(minValue_ * (height()-leftTopOffset_-rightBottomOffset_) + leftTopOffset_);
        rightBottomMarker = tgt::iround(maxValue_ * (height()-leftTopOffset_-rightBottomOffset_) + leftTopOffset_);
    } else {
        normalizedMousePos_ = static_cast<float>(e->pos().x()-leftTopOffset_) / static_cast<float>(width()-leftTopOffset_-rightBottomOffset_);
        leftTopMarker =  tgt::iround(minValue_ * (width()-leftTopOffset_-rightBottomOffset_) + leftTopOffset_);
        rightBottomMarker = tgt::iround(maxValue_ * (width()-leftTopOffset_-rightBottomOffset_) + leftTopOffset_);
    }

    mV1_ = minValue_; mV2_ = maxValue_;
    if (e->button() == Qt::LeftButton) {
        int pos = (vertical_ ? e->pos().y() : e->pos().x());
        if (pos < (leftTopMarker + sliderWidthHeight_)) {
            leftTopSliderActive_ = true;
            rightBottomSliderActive_ = false;
        }
        else if(pos > (rightBottomMarker - sliderWidthHeight_)) {
            leftTopSliderActive_ = false;
            rightBottomSliderActive_ = true;
        }
        else if((leftTopMarker < pos) && (rightBottomMarker > pos)){
            leftTopSliderActive_ = true;
            rightBottomSliderActive_ = true;
        }
    }
    else if (e->button() == Qt::RightButton) {
        leftTopSliderActive_ = true;
        rightBottomSliderActive_ = true;
    }
    moveSlider(normalizedMousePos_);
    emit toggleInteractionMode(true);
}

void DoubleSlider::mouseMoveEvent(QMouseEvent* e){
    e->accept();
    globalMousePos_ = e->globalPos();
    float normalizedMousePosTmp;
    if(vertical_)
        normalizedMousePosTmp = static_cast<float>((e->pos()).y()-leftTopOffset_) / static_cast<float>(height()-leftTopOffset_-rightBottomOffset_);
    else
        normalizedMousePosTmp = static_cast<float>((e->pos()).x()-leftTopOffset_) / static_cast<float>(width()-leftTopOffset_-rightBottomOffset_);

    if (normalizedMousePosTmp > 1.f)
        normalizedMousePosTmp = 1.f;
    else if (normalizedMousePosTmp < 0.f)
        normalizedMousePosTmp = 0.f;
    moveSlider(normalizedMousePosTmp);

    emit valuesChanged(minValue_, maxValue_);
}

void DoubleSlider::showToolTip(std::string text) {
    if(!showToolTip_)
        return;

    //QToolTip::showText(globalMousePos_, "");
    QToolTip::showText(globalMousePos_, QString::fromStdString(text));
}

void DoubleSlider::mouseReleaseEvent(QMouseEvent* e) {
    leftTopSliderActive_ = false;
    rightBottomSliderActive_ = false;
    e->accept();
    emit toggleInteractionMode(false);
}

std::string DoubleSlider::generateToolTipText(bool minBold, bool maxBold) {
    std::string text = "<p style='white-space:pre'>";

    if(minBold)
        text += "<b>";
    text += ftos(getMappedValue(getMinValue()));
    if(minBold)
        text += "</b>";

    text += " -> ";

    if(maxBold)
        text += "<b>";
    text += ftos(getMappedValue(getMaxValue()));
    if(maxBold)
        text += "</b>";

    if(unit_ != "")
        text += " [" + unit_ + "]";

    text += "<br>w: " + ftos(getMappedValue(getMaxValue()) - getMappedValue(getMinValue()));
    text += "<br>l: " + ftos( (getMappedValue(getMaxValue()) + getMappedValue(getMinValue())) * 0.5f);

    return text;
}

void DoubleSlider::moveSlider(float mousePos) {
    if (leftTopSliderActive_ && !rightBottomSliderActive_) {
        setMinValue(mousePos);
        showToolTip(generateToolTipText(true, false));
    }
    if (rightBottomSliderActive_ && !leftTopSliderActive_) {
        setMaxValue(mousePos);
        showToolTip(generateToolTipText(false, true));
    }
    if (rightBottomSliderActive_ && leftTopSliderActive_) {
        float mouseDiff = normalizedMousePos_ - mousePos;
        setMinValue(mV1_ - mouseDiff);
        setMaxValue(mV2_ - mouseDiff);
        showToolTip(generateToolTipText(true, true));
    }
}

void DoubleSlider::setMinValue(float val) {
    if (val == minValue_)
        return;

    if (val < 0.f)
        val = 0.f;
    if (val + minimalAllowedSliderDistance_ < maxValue_)
        minValue_ = val;
    else {
        maxValue_ = val + minimalAllowedSliderDistance_;
        if (maxValue_ > 1.f) {
            maxValue_ = 1.f;
            minValue_ = 1.f - minimalAllowedSliderDistance_;
        }
        else
            minValue_ = val;
    }
    setToolTip(QString::fromStdString(generateToolTipText(false, false)));
    update();
    emit valuesChanged(minValue_, maxValue_);
}

void DoubleSlider::setMaxValue(float val) {
    if (val == maxValue_)
        return;

    if (val > 1.f)
        val = 1.f;
    if (minValue_ + minimalAllowedSliderDistance_ < val)
        maxValue_ = val;
    else {
        minValue_ = val - minimalAllowedSliderDistance_;
        if (minValue_ < 0.f) {
            minValue_ = 0.f;
            maxValue_ = minimalAllowedSliderDistance_;
        }
        else
            maxValue_ = val;
    }
    setToolTip(QString::fromStdString(generateToolTipText(false, false)));
    update();
    emit valuesChanged(minValue_, maxValue_);
}

void DoubleSlider::setSliderWidthOrHeight(int sliderWH) {
    sliderWidthHeight_ = sliderWH;
}

void DoubleSlider::setValues(float val1, float val2) {
    if (val1 < val2) {
        setMinValue(val1);
        setMaxValue(val2);
    }
    else {
        setMinValue(val2);
        setMaxValue(val1);
    }
}

float DoubleSlider::getMinValue() {
    return minValue_;
}

float DoubleSlider::getMaxValue() {
    return maxValue_;
}

void DoubleSlider::setMapping(float min, float max) {
    minRWValue_ = min;
    maxRWValue_ = max;
    setToolTip(QString::fromStdString(generateToolTipText(false, false)));
}

void DoubleSlider::setUnit(std::string unit) {
    unit_ = unit;
    setToolTip(QString::fromStdString(generateToolTipText(false, false)));
}

void DoubleSlider::showToolTip(bool stt) {
    showToolTip_ = stt;
}

float DoubleSlider::getMappedValue(float norm) {
    return minRWValue_ + (norm * (maxRWValue_ - minRWValue_));
}

void DoubleSlider::setMinimalAllowedSliderDistance(float dist) {
    minimalAllowedSliderDistance_ = dist;
}

} // namespace voreen
