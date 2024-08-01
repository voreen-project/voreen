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

#include "voreen/qt/widgets/property/intboundingboxpropertywidget.h"

#include "voreen/core/properties/vectorproperty.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include <QMouseEvent>
#include <QMenu>

namespace voreen {

IntBoundingBoxPropertyWidget::IntBoundingBoxPropertyWidget(IntBoundingBoxProperty* const prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , prop_(prop)
    , xProperty_("", "", tgt::vec2(prop->get().getLLF().x, prop->get().getURB().x), prop->getMinValue().x,
                 prop->getMaxValue().x, prop->getMinRange().x, prop->getMaxRange().x)
    , yProperty_("", "", tgt::vec2(prop->get().getLLF().y, prop->get().getURB().y), prop->getMinValue().y,
                 prop->getMaxValue().y, prop->getMinRange().y, prop->getMaxRange().y)
    , zProperty_("", "", tgt::vec2(prop->get().getLLF().z, prop->get().getURB().z), prop->getMinValue().z,
                 prop->getMaxValue().z, prop->getMinRange().z, prop->getMaxRange().z)
{
    xWidget_ = new IntIntervalPropertyWidget(&xProperty_, this);
    yWidget_ = new IntIntervalPropertyWidget(&yProperty_, this);
    zWidget_ = new IntIntervalPropertyWidget(&zProperty_, this);

    xProperty_.addWidget(xWidget_);
    yProperty_.addWidget(yWidget_);
    zProperty_.addWidget(zWidget_);

    groupBoxX_ = new QGroupBox(this);
    groupBoxY_ = new QGroupBox(this);
    groupBoxZ_ = new QGroupBox(this);
    vboxX_ = new QVBoxLayout();
    vboxY_ = new QVBoxLayout();
    vboxZ_ = new QVBoxLayout();
    vbox_     = new QVBoxLayout();

    xLabel_ = new QLabel(tr("x component:"));
    yLabel_ = new QLabel(tr("y component:"));
    zLabel_ = new QLabel(tr("z component:"));

    connect(xWidget_, SIGNAL(valueChanged(tgt::ivec2)), this, SLOT(valueChanged(tgt::ivec2)));
    connect(yWidget_, SIGNAL(valueChanged(tgt::ivec2)), this, SLOT(valueChanged(tgt::ivec2)));
    connect(zWidget_, SIGNAL(valueChanged(tgt::ivec2)), this, SLOT(valueChanged(tgt::ivec2)));


    addLayout(vbox_);

    vbox_->addWidget(groupBoxX_);
    vbox_->addWidget(groupBoxY_);
    vbox_->addWidget(groupBoxZ_);

    groupBoxX_->setLayout(vboxX_);
    groupBoxY_->setLayout(vboxY_);
    groupBoxZ_->setLayout(vboxZ_);

    vboxX_->addWidget(xLabel_);
    vboxX_->addWidget(xWidget_);
    vboxY_->addWidget(yLabel_);
    vboxY_->addWidget(yWidget_);
    vboxZ_->addWidget(zLabel_);
    vboxZ_->addWidget(zWidget_);

    updateFromPropertySlot();

}

IntBoundingBoxPropertyWidget::~IntBoundingBoxPropertyWidget() {
    delete xWidget_;
    delete yWidget_;
    delete zWidget_;
    delete xLabel_;
    delete yLabel_;
    delete zLabel_;
    delete vbox_;
}

void IntBoundingBoxPropertyWidget::updateFromPropertySlot() {
    tgt::IntBounds bounds = prop_->get();
    tgt::ivec3 min = bounds.getLLF();
    tgt::ivec3 max = bounds.getURB();

    tgt::vec2 xminmax = tgt::vec2(min.x, max.x);
    tgt::vec2 yminmax = tgt::vec2(min.y, max.y);
    tgt::vec2 zminmax = tgt::vec2(min.z, max.z);

    tgt::ivec3 minValue = prop_->getMinValue();
    tgt::ivec3 maxValue = prop_->getMaxValue();
    tgt::ivec3 stepping = prop_->getStepping();

    bool tracking = prop_->hasTracking();
    int numDecimals = prop_->getNumDecimals();

    xProperty_.setMinValue(minValue.x);
    xProperty_.setMaxValue(maxValue.x);
    xProperty_.setStepping(stepping.x);
    xProperty_.setTracking(tracking);
    xProperty_.setNumDecimals(numDecimals);

    yProperty_.setMinValue(minValue.y);
    yProperty_.setMaxValue(maxValue.y);
    yProperty_.setStepping(stepping.y);
    yProperty_.setTracking(tracking);
    yProperty_.setNumDecimals(numDecimals);

    zProperty_.setMinValue(minValue.z);
    zProperty_.setMaxValue(maxValue.z);
    zProperty_.setStepping(stepping.z);
    zProperty_.setTracking(tracking);
    zProperty_.setNumDecimals(numDecimals);

    xProperty_.set(xminmax);
    yProperty_.set(yminmax);
    zProperty_.set(zminmax);
}

void IntBoundingBoxPropertyWidget::valueChanged(tgt::ivec2){
    tgt::vec2 xminmax = xProperty_.get();
    tgt::vec2 yminmax = yProperty_.get();
    tgt::vec2 zminmax = zProperty_.get();

    tgt::ivec3 min = tgt::ivec3(xminmax.x, yminmax.x, zminmax.x);
    tgt::ivec3 max = tgt::ivec3(xminmax.y, yminmax.y, zminmax.y);

    tgt::IntBounds bounds;
    bounds.addPoint(min);
    bounds.addPoint(max);

    prop_->set(bounds);
    emit valueModifiedByUser();
}
}   // namespace
