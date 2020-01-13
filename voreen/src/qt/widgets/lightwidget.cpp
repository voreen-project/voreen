/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/floatvec4propertywidget.h"
#include "voreen/qt/widgets/lightwidget.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QStyle>
#include <QWidget>
#include <QPushButton>
#include <QSlider>

#include "math.h"

namespace voreen {

LightWidget::LightWidget(QWidget* parent)
    : QWidget(parent)
    , radius_(50)
    , hemisphere_(false)
    , distance_(4.f)
    , offset_(tgt::vec2(0.f))
{
    distanceSlider_ = new DoubleSliderSpinBoxWidget(this);
    distanceSlider_->setMaxValue(10.0f);
    distanceSlider_->setMinValue(1.0f);
    distanceSlider_->setValue(distance_);
    connect(distanceSlider_, SIGNAL(valueChanged(double)), this, SLOT(updateDistance(double)));
    distanceSlider_->move(QPoint(0,110));
    distanceSlider_->setMaximumWidth(110);
    distanceSlider_->setSingleStep(0.1f);
    distanceSlider_->show();

    hemisphereButton_ = new QPushButton(this);
    hemisphereButton_->setCheckable(true);
    hemisphereButton_->setChecked(hemisphere_);
    hemisphereButton_->setText("Z");
    hemisphereButton_->setGeometry(QRect(0, 0, 20, 20));
    connect(hemisphereButton_, SIGNAL(toggled(bool)), this, SLOT(setHemisphere(bool)));
    hemisphereButton_->move(QPoint(120, 110));
    hemisphereButton_->show();

    lightPosition_ = QPointF(50,50);

    setMinimumSize(75, 140);
}

void LightWidget::paintEvent(QPaintEvent* /*event*/) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    QRadialGradient radialGradient(QPointF(50, 50), radius_ * 2  , QPointF(lightPosition_.x(), lightPosition_.y()));

    float colDist = (distance_ / distanceSlider_->getMaxValue()) * 10.f;

    QColor highlight(static_cast<int>(colDist*3)+100,
                     static_cast<int>(colDist*3)+100,
                     static_cast<int>(colDist*3)+100);
    if (!hemisphere_) {
        radialGradient.setColorAt(0.0, QColor(255,255,255));
        radialGradient.setColorAt((colDist /30)+ 0.5, highlight);
        radialGradient.setColorAt(1.0, QColor(0,0,0));
    }
    else {
        radialGradient.setColorAt(0.0, QColor(0,0,0));
        radialGradient.setColorAt((colDist /30) + 0.5, highlight);
        radialGradient.setColorAt(1.0, QColor(255,255,255));
    }

    painter.setBrush(QBrush(radialGradient));
    painter.setPen(highlight);
    painter.drawEllipse(QRectF(2, 2, radius_ * 2, radius_ * 2));
}

void LightWidget::setHemisphere(bool hemibool) {
    hemisphere_ = hemibool;
    int hem = 1;
    if(hemisphere_)
        hem = -1;
    emit(lightWidgetChanged(distance_ * tgt::vec4(offset_, hem * sqrt(std::max(0.f, 1.f - lengthSq(offset_))), 0.f)));
    update();
}

void LightWidget::setMinDist(double d) {
    distanceSlider_->setMinValue(d);
}

void LightWidget::setMaxDist(double d) {
    distanceSlider_->setMaxValue(d);
}

float LightWidget::getMaxDist() const {
    return distanceSlider_->getMaxValue();
}

void LightWidget::mousePressEvent(QMouseEvent *event) {
    mouseMoveEvent(event);
    update();
}

void LightWidget::mouseMoveEvent(QMouseEvent *event) {
    tgt::vec2 tmp = tgt::vec2(event->pos().x() - 60, event->pos().y() - 60) / 60.f;
    if(lengthSq(tmp) > 1.f)
        tmp = normalize(tmp);

    int hem = 1;
    if(hemisphere_) {
        hem = -1;
    }

    lightPosition_ = QPointF((tmp.x + 1.f)*60.f, (tmp.y + 1.f)*60.f);
    offset_ = tgt::vec2(tmp.x, -1.f * tmp.y);

    emit(lightWidgetChanged(distance_ * tgt::vec4(offset_, hem * sqrt(std::max(0.f, 1.f - lengthSq(offset_))), 0.f)));
    update();
}

void LightWidget::updateDistance(double distance) {
    distance_ = distance;
    update();
    int hem = 1;
    if(hemisphere_) {
        hem = -1;
    }
    emit(lightWidgetChanged(distance_ * tgt::vec4(offset_, hem * sqrt(std::max(0.f, 1.f - lengthSq(offset_))), 0.f)));
}

void LightWidget::setLightPosition(const tgt::vec4& position) {
    hemisphere_ = position.z < 0.f ? true : false;
    distance_ = std::max((float)(distanceSlider_->getMinValue()), std::min((float)(distanceSlider_->getMaxValue()), length(position.xyz())));

    tgt::vec2 pos(0.f);
    if(length(position.xyz()) > 0.0001f) {
        offset_ = position.xy() / length(position.xyz());
        if(lengthSq(offset_) > 1.f)
            offset_ = normalize(offset_);
    }
    else
        offset_ = tgt::vec2(0.0f);

    lightPosition_ = QPointF((offset_.x + 1.f)*60.f, (1.f - offset_.y)*60.f);

    distanceSlider_->blockSignals(true);
    distanceSlider_->setValue(distance_);
    distanceSlider_->blockSignals(false);
    hemisphereButton_->blockSignals(true);
    hemisphereButton_->setChecked(hemisphere_);
    hemisphereButton_->blockSignals(false);
    update();
}

tgt::vec4 LightWidget::getLightPosition() const {
    return distance_ * tgt::vec4(offset_, (hemisphere_ ? -1.f : 1.f) * sqrt(std::max(0.f, 1.f - lengthSq(offset_))), 0.f);
}

bool LightWidget::getHemisphereStatus() const {
    return hemisphereButton_->isChecked();
}

} // namespace voreen
