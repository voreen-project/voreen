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

#include "voreen/qt/widgets/property/floatmat2propertywidget.h"
#include "voreen/core/utils/serializationhelper.h"

#include <QDialog>
#include <QDoubleValidator>
#include <QGridLayout>
#include <QLineEdit>
#include <QMenu>

namespace voreen {

FloatMat2PropertyWidget::FloatMat2PropertyWidget(FloatMat2Property* prop, QWidget* parent)
    : MatrixPropertyWidget(prop, 2, parent)
    , prop_(prop)
{
    connect(matrixDialog_, SIGNAL(customContextMenuRequested(const QPoint &)) , this, SLOT(precisionMenu(const QPoint &)));
    x0_ = new QLineEdit(this);
    y0_ = new QLineEdit(this);
    x1_ = new QLineEdit(this);
    y1_ = new QLineEdit(this);

    x0_->setValidator(doubleValidator_);
    y0_->setValidator(doubleValidator_);
    x1_->setValidator(doubleValidator_);
    y1_->setValidator(doubleValidator_);

    connect(x0_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(y0_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(x1_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(y1_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));

    mainLayout_->addWidget(x0_, 1, 0);
    mainLayout_->addWidget(y0_, 1, 1);
    mainLayout_->addWidget(x1_, 2, 0);
    mainLayout_->addWidget(y1_, 2, 1);
    updateFromPropertySlot();
}

void FloatMat2PropertyWidget::updateFromPropertySlot() {
    doubleValidator_->setTop(prop_->getMaxValue().t00);
    doubleValidator_->setBottom(prop_->getMinValue().t00);

    // prevent property update during updating widgets
    x0_->blockSignals(true);
    y0_->blockSignals(true);
    x1_->blockSignals(true);
    y1_->blockSignals(true);

    x0_->setText(QString::number(prop_->get().t00));
    y0_->setText(QString::number(prop_->get().t01));
    x1_->setText(QString::number(prop_->get().t10));
    y1_->setText(QString::number(prop_->get().t11));

    // reactivate widgets
    x0_->blockSignals(false);
    y0_->blockSignals(false);
    x1_->blockSignals(false);
    y1_->blockSignals(false);
}

void FloatMat2PropertyWidget::updateValue() {
    tgt::mat2 matrix;
    matrix.t00 = x0_->text().toDouble();
    matrix.t01 = y0_->text().toDouble();
    matrix.t10 = x1_->text().toDouble();
    matrix.t11 = y1_->text().toDouble();
    prop_->set(matrix);
    emit valueModifiedByUser();
}

void FloatMat2PropertyWidget::identity() {
    x0_->setText("1");
    y0_->setText("0");
    x1_->setText("0");
    y1_->setText("1");
}

void FloatMat2PropertyWidget::saveMatrix(const std::string& filename) const {
    tgt::mat2 matrix = prop_->get();
    SerializationHelper::save<tgt::mat2>(filename, "Matrix", matrix);
}

void FloatMat2PropertyWidget::loadMatrix(const std::string& filename) {
    tgt::mat2 matrix;
    SerializationHelper::load<tgt::mat2>(filename, "Matrix", matrix);
    prop_->set(matrix);
    emit valueModifiedByUser();
}

} //namespace voreen
