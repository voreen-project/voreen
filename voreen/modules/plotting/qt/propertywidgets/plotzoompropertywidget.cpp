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

#include "voreen/qt/widgets/property/plotzoompropertywidget.h"

#include "voreen/core/plotting/plotpredicate.h"
#include "voreen/core/properties/plotselectionproperty.h"

#include <QLabel>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QVBoxLayout>
#include <QGridLayout>
#include <sstream>
#include <limits>

namespace voreen {

PlotZoomPropertyWidget::PlotZoomPropertyWidget(PlotSelectionProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
    , lblX_(0)
    , lblY_(0)
    , lblZ_(0)
    , btnPush_(0)
    , btnPop_(0)
    , spXlow_(0)
    , spXhigh_(0)
    , spYlow_(0)
    , spYhigh_(0)
    , spZlow_(0)
    , spZhigh_(0)
    {
    QWidget* mainWidget = new QWidget(this);
    QGridLayout* gl = new QGridLayout();

    lblX_ = new QLabel(tr("x Axis:"));
    spXlow_ = new QDoubleSpinBox();
    spXlow_->setDecimals(4);
    spXhigh_ = new QDoubleSpinBox();
    spXhigh_->setDecimals(4);

    gl->addWidget(lblX_, 0, 0);
    gl->addWidget(spXlow_, 0, 1);
    gl->addWidget(spXhigh_, 0, 2);

    lblY_ = new QLabel(tr("y Axis:"));
    spYlow_ = new QDoubleSpinBox();
    spYlow_->setDecimals(4);
    spYhigh_ = new QDoubleSpinBox();
    spYhigh_->setDecimals(4);

    gl->addWidget(lblY_, 1, 0);
    gl->addWidget(spYlow_, 1, 1);
    gl->addWidget(spYhigh_, 1, 2);

    lblZ_ = new QLabel(tr("z Axis:"));
    spZlow_ = new QDoubleSpinBox();
    spZlow_->setDecimals(4);
    spZhigh_ = new QDoubleSpinBox();
    spZhigh_->setDecimals(4);

    gl->addWidget(lblZ_, 2, 0);
    gl->addWidget(spZlow_, 2, 1);
    gl->addWidget(spZhigh_, 2, 2);

    btnPush_ =  new QPushButton(tr("Apply"));
    btnPop_ =  new QPushButton(tr("Undo"));
    gl->addWidget(btnPush_, 3, 1);
    gl->addWidget(btnPop_, 3, 2);

    connect(btnPop_, SIGNAL(clicked()), this, SLOT(btnPopClicked()));
    connect(btnPush_, SIGNAL(clicked()), this, SLOT(btnPushClicked()));

    connect(spXlow_, SIGNAL(editingFinished()), this, SLOT(setSpinBoxRanges()));
    connect(spXhigh_, SIGNAL(editingFinished()), this, SLOT(setSpinBoxRanges()));
    connect(spYlow_, SIGNAL(editingFinished()), this, SLOT(setSpinBoxRanges()));
    connect(spYhigh_, SIGNAL(editingFinished()), this, SLOT(setSpinBoxRanges()));

    updateFromProperty();

    mainWidget->setLayout(gl);
    addWidget(mainWidget);

    addVisibilityControls();
}

void PlotZoomPropertyWidget::updateFromProperty() {
    if (! property_->get().empty()) {
        spXlow_->setValue(property_->top().xZoom_.getLeft());
        spXhigh_->setValue(property_->top().xZoom_.getRight());
        spYlow_->setValue(property_->top().yZoom_.getLeft());
        spYhigh_->setValue(property_->top().yZoom_.getRight());
        spZlow_->setValue(property_->top().zZoom_.getLeft());
        spZhigh_->setValue(property_->top().zZoom_.getRight());
        setSpinBoxRanges();
    }
    if (property_->getShowThirdDimension()) {
        lblZ_->setVisible(true);
        spZlow_->setVisible(true);
        spZhigh_->setVisible(true);
    }
    else {
        lblZ_->setVisible(false);
        spZlow_->setVisible(false);
        spZhigh_->setVisible(false);
    }
}

void PlotZoomPropertyWidget::btnPushClicked() {
    if (checkSanityOfInput()) {
        property_->push(PlotZoomState(Interval<plot_t>(spXlow_->value(), spXhigh_->value()),
                                      Interval<plot_t>(spYlow_->value(), spYhigh_->value()),
                                      Interval<plot_t>(spZlow_->value(), spZhigh_->value())));
    }
}

void PlotZoomPropertyWidget::btnPopClicked() {
    if (property_->size() > 1)
        property_->pop();
}

bool PlotZoomPropertyWidget::checkSanityOfInput() {
    if (spXlow_->value() > spXhigh_->value())
        return false;
    if (spYlow_->value() > spYhigh_->value())
        return false;
    if (spZlow_->value() > spZhigh_->value())
        return false;
    return true;
}

void PlotZoomPropertyWidget::setSpinBoxRanges() {
    // REMARK: std::numeric_limits<double>::min() does not work here :( hopefully -1000000 is low enough
    spXlow_->setMinimum(-1000000);
    spXlow_->setMaximum(spXhigh_->value());

    spXhigh_->setMinimum(spXlow_->value());
    spXhigh_->setMaximum(std::numeric_limits<double>::max());

    spYlow_->setMinimum(-1000000);
    spYlow_->setMaximum(spYhigh_->value());

    spYhigh_->setMinimum(spYlow_->value());
    spYhigh_->setMaximum(std::numeric_limits<double>::max());

    spZlow_->setMinimum(-1000000);
    spZlow_->setMaximum(spZhigh_->value());

    spZhigh_->setMinimum(spZlow_->value());
    spZhigh_->setMaximum(std::numeric_limits<double>::max());
}

} // namespace
