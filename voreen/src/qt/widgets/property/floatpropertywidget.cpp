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

#include "voreen/qt/widgets/property/floatpropertywidget.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include <QMenu>
#include <QMouseEvent>

namespace voreen {

FloatPropertyWidget::FloatPropertyWidget(FloatProperty* prop, QWidget* parent, bool addVisibilityControl)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , widget_(new DoubleSliderSpinBoxWidget(this))
{
    tgtAssert(property_, "no property");
    widget_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    widget_->setSliderTracking(property_->hasTracking());
    QPropertyWidget::layout_->addWidget(widget_, Qt::AlignLeft);

    updateFromPropertySlot();

    connect(widget_, SIGNAL(valueChanged(double)), this, SLOT(setProperty(double)));
    connect(widget_, SIGNAL(sliderPressedChanged(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(widget_, SIGNAL(valueChanged(double)), this, SIGNAL(widgetChanged()));

    precisionMenu_ = new QMenu(this);
    highAction_ = precisionMenu_->addAction("High Precision");
    highAction_->setCheckable(true);
    precisionMenu_->addSeparator();
    instantValueChangeAction_ = precisionMenu_->addAction("Tracking Mode");
    instantValueChangeAction_->setCheckable(true);
    if(property_->hasTracking())
        instantValueChangeAction_->toggle();
}

FloatPropertyWidget::~FloatPropertyWidget() {
    widget_->disconnect();
    delete widget_;
}

void FloatPropertyWidget::updateFromPropertySlot() {
    if (property_ != 0) {
        widget_->blockSignals(true);
        widget_->setDecimals(property_->getNumDecimals());
        widget_->setMinValue(property_->getMinValue());
        widget_->setMaxValue(property_->getMaxValue());
        widget_->setSingleStep(property_->getStepping());
        widget_->setValue(property_->get());
        widget_->setSliderTracking(property_->hasTracking());
        widget_->blockSignals(false);
    }
}

void FloatPropertyWidget::setWidget(const float value, const float minValue,
                                  const float maxValue, const float stepping)
{
    widget_->blockSignals(true);
    widget_->setMinValue(static_cast<const double>(minValue));
    widget_->setMaxValue(static_cast<const double>(maxValue));
    widget_->setSingleStep(static_cast<const double>(stepping));
    widget_->setValue(static_cast<const double>(value));
    widget_->blockSignals(false);
}

void FloatPropertyWidget::setProperty(double value) {
    if ((property_ != 0) && (!disconnected_)) {
        property_->set(value);
        emit valueModifiedByUser();
    }

    emit valueChanged(value);
}

void FloatPropertyWidget::mousePressEvent(QMouseEvent* event) {
    if(event->button() == Qt::RightButton) {

        QAction* prec = precisionMenu_->exec(mapToGlobal(event->pos()));
        if (prec == highAction_){
            if(highAction_->isChecked()) {
                property_->setStepping(0.0001f);
                property_->setNumDecimals(4);
            }
            else {
                property_->setStepping(0.05f);
                property_->setNumDecimals(2);
            }
        }
        else if(prec == instantValueChangeAction_) {
            property_->setTracking(!property_->hasTracking());
        }
        updateFromProperty();
    }
    QWidget::mousePressEvent(event);
}

} // namespace
