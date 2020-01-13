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

#include "voreen/qt/widgets/property/intpropertywidget.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include <QHBoxLayout>
#include <QMouseEvent>
#include <QMenu>

namespace voreen {

IntPropertyWidget::IntPropertyWidget(IntProperty* prop, QWidget* parent, bool addVisibilityControl)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , widget_(new SliderSpinBoxWidget(this))
{
    tgtAssert(property_, "no property");
    widget_->setSliderTracking(property_->hasTracking());
    addWidget(widget_);

    updateFromPropertySlot();

    connect(widget_, SIGNAL(valueChanged(int)), this, SLOT(setProperty(int)));
    connect(widget_, SIGNAL(sliderPressedChanged(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(widget_, SIGNAL(valueChanged(int)), this, SIGNAL(widgetChanged()));

    instantValueChangeMenu_ = new QMenu(this);
    instantValueChangeAction_ = instantValueChangeMenu_->addAction("Tracking Mode");
    instantValueChangeAction_->setCheckable(true);
    if (property_->hasTracking())
        instantValueChangeAction_->toggle();
}

IntPropertyWidget::~IntPropertyWidget() {
    widget_->disconnect();
    delete widget_;
}

void IntPropertyWidget::updateFromPropertySlot() {
    if (property_ != 0) {
        widget_->blockSignals(true);
        widget_->setMinValue(property_->getMinValue());
        widget_->setMaxValue(property_->getMaxValue());
        widget_->setSingleStep(property_->getStepping());
        widget_->setValue(property_->get());
        widget_->setSliderTracking(property_->hasTracking());
        widget_->blockSignals(false);
    }
}

void IntPropertyWidget::setWidget(const int value, const int minValue,
                                  const int maxValue, const int stepping)
{
    widget_->blockSignals(true);
    widget_->setMinValue(minValue);
    widget_->setMaxValue(maxValue);
    widget_->setSingleStep(stepping);
    widget_->setValue(value);
    widget_->blockSignals(false);
}

void IntPropertyWidget::setProperty(int value) {
    if ((property_ != 0) && (!disconnected_)) {
        property_->set(value);
        emit valueModifiedByUser();
    }

    emit valueChanged(value);
}

void IntPropertyWidget::mousePressEvent(QMouseEvent* event) {
    if(event->button() == Qt::RightButton) {
        QAction* prec = instantValueChangeMenu_->exec(mapToGlobal(event->pos()));
        if(prec == instantValueChangeAction_) {
            property_->setTracking(!property_->hasTracking());
        }
        updateFromPropertySlot();
    }
    QWidget::mousePressEvent(event);
}

} // namespace
