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

#include "voreen/qt/widgets/property/numeric/intintervalpropertywidget.h"

#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include <QHBoxLayout>
#include <QMouseEvent>
#include <QMenu>

namespace voreen {

IntIntervalPropertyWidget::IntIntervalPropertyWidget(IntIntervalProperty* prop, QWidget* parent, bool addVisibilityControl)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , minWidget_(new SliderSpinBoxWidget(this))
    , maxWidget_(new SliderSpinBoxWidget(this))
{
    tgtAssert(property_, "no property");

    vboxLayout_ = new QVBoxLayout();
    vboxLayout_->addWidget(minWidget_);
    vboxLayout_->addWidget(maxWidget_);
    addLayout(vboxLayout_);

    updateFromPropertySlot();

    connect(minWidget_, SIGNAL(valueChanged(int)), this, SLOT(setPropertyMin(int)));
    connect(minWidget_, SIGNAL(sliderPressedChanged(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(minWidget_, SIGNAL(valueChanged(int)), this, SIGNAL(widgetChanged()));

    connect(maxWidget_, SIGNAL(valueChanged(int)), this, SLOT(setPropertyMax(int)));
    connect(maxWidget_, SIGNAL(sliderPressedChanged(bool)), this, SLOT(toggleInteractionMode(bool)));
    connect(maxWidget_, SIGNAL(valueChanged(int)), this, SIGNAL(widgetChanged()));

    instantValueChangeMenu_ = new QMenu(this);
    instantValueChangeAction_ = instantValueChangeMenu_->addAction("Tracking Mode");
    instantValueChangeAction_->setCheckable(true);
}

IntIntervalPropertyWidget::~IntIntervalPropertyWidget() {
    minWidget_->disconnect();
    maxWidget_->disconnect();
    delete minWidget_;
    delete maxWidget_;
}

void IntIntervalPropertyWidget::updateFromPropertySlot() {
    if (property_ != 0) {
        minWidget_->blockSignals(true);
        minWidget_->setMinValue(property_->getMinValue());
        minWidget_->setMaxValue(property_->getMaxValue());
        minWidget_->setValue(property_->get().x);
        minWidget_->setSingleStep(property_->getStepping());
        minWidget_->blockSignals(false);

        maxWidget_->blockSignals(true);
        maxWidget_->setMinValue(property_->getMinValue());
        maxWidget_->setMaxValue(property_->getMaxValue());
        maxWidget_->setValue(property_->get().y);
        maxWidget_->setSingleStep(property_->getStepping());
        maxWidget_->blockSignals(false);
    }
}

void IntIntervalPropertyWidget::setPropertyMin(int newValue) {
    if ((property_ == 0) && disconnected_){
        return;
    }

    tgt::Interval<int> interval = property_->getInterval();
    int min = newValue;
    int max = maxWidget_->getValue();

    tgt::ivec2 value = tgt::ivec2(min, max);
    // if possible only fixup max
    value = interval.fixupMax(value);
    // if not als min
    value = interval.fixupMin(value);

    property_->set(value);

    minWidget_->blockSignals(true);
    minWidget_->setValue(value.x);
    minWidget_->blockSignals(false);
    maxWidget_->blockSignals(true);
    maxWidget_->setValue(value.y);
    maxWidget_->blockSignals(false);

    emit valueModifiedByUser();
    emit valueChanged(value);
}

void IntIntervalPropertyWidget::setPropertyMax(int newValue) {
    if ((property_ == 0) && disconnected_){
        return;
    }

    tgt::Interval<int> interval = property_->getInterval();
    int min = minWidget_->getValue();
    int max = newValue;

    tgt::ivec2 value = tgt::ivec2(min, max);
    // if possible only fixup min
    value = interval.fixupMin(value);
    // if not als max
    value = interval.fixupMax(value);

    property_->set(value);

    minWidget_->blockSignals(true);
    minWidget_->setValue(value.x);
    minWidget_->blockSignals(false);
    maxWidget_->blockSignals(true);
    maxWidget_->setValue(value.y);
    maxWidget_->blockSignals(false);

    emit valueModifiedByUser();
    emit valueChanged(value);
}

} // namespace
