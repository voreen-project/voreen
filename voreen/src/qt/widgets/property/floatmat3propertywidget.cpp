/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/qt/widgets/property/floatmat3propertywidget.h"
#include "voreen/core/utils/serializationhelper.h"

#include <QGridLayout>
#include <QLineEdit>
#include <QMenu>
#include <QDialog>

namespace voreen {

FloatMat3PropertyWidget::FloatMat3PropertyWidget(FloatMat3Property* prop, QWidget* parent)
    : MatrixPropertyWidget(prop, 3, parent)
    , prop_(prop)
{
    connect(matrixDialog_, SIGNAL(customContextMenuRequested(const QPoint &)) , this, SLOT(precisionMenu(const QPoint &)));
    x0_ = new QLineEdit(this);
    y0_ = new QLineEdit(this);
    z0_ = new QLineEdit(this);
    x1_ = new QLineEdit(this);
    y1_ = new QLineEdit(this);
    z1_ = new QLineEdit(this);
    x2_ = new QLineEdit(this);
    y2_ = new QLineEdit(this);
    z2_ = new QLineEdit(this);

    x0_->setValidator(doubleValidator_);
    y0_->setValidator(doubleValidator_);
    z0_->setValidator(doubleValidator_);
    x1_->setValidator(doubleValidator_);
    y1_->setValidator(doubleValidator_);
    z1_->setValidator(doubleValidator_);
    x2_->setValidator(doubleValidator_);
    y2_->setValidator(doubleValidator_);
    z2_->setValidator(doubleValidator_);

    connect(x0_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(y0_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(z0_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(x1_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(y1_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(z1_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(x2_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(y2_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));
    connect(z2_, SIGNAL(textEdited(const QString &)), this, SLOT(autoUpdate()));

    mainLayout_->addWidget(x0_, 1, 0);
    mainLayout_->addWidget(y0_, 1, 1);
    mainLayout_->addWidget(z0_, 1, 2);
    mainLayout_->addWidget(x1_, 2, 0);
    mainLayout_->addWidget(y1_, 2, 1);
    mainLayout_->addWidget(z1_, 2, 2);
    mainLayout_->addWidget(x2_, 3, 0);
    mainLayout_->addWidget(y2_, 3, 1);
    mainLayout_->addWidget(z2_, 3, 2);

    updateFromPropertySlot();
}

void FloatMat3PropertyWidget::updateFromPropertySlot() {
    doubleValidator_->setTop(prop_->getMaxValue().t00);
    doubleValidator_->setBottom(prop_->getMinValue().t00);

    // prevent property update during updating widgets
    x0_->blockSignals(true);
    y0_->blockSignals(true);
    z0_->blockSignals(true);
    x1_->blockSignals(true);
    y1_->blockSignals(true);
    z1_->blockSignals(true);
    x2_->blockSignals(true);
    y2_->blockSignals(true);
    z2_->blockSignals(true);

    x0_->setText(QString::number(prop_->get().t00));
    y0_->setText(QString::number(prop_->get().t01));
    z0_->setText(QString::number(prop_->get().t02));
    x1_->setText(QString::number(prop_->get().t10));
    y1_->setText(QString::number(prop_->get().t11));
    z1_->setText(QString::number(prop_->get().t12));
    x2_->setText(QString::number(prop_->get().t20));
    y2_->setText(QString::number(prop_->get().t21));
    z2_->setText(QString::number(prop_->get().t22));

    // reactivate widgets
    x0_->blockSignals(false);
    y0_->blockSignals(false);
    z0_->blockSignals(false);
    x1_->blockSignals(false);
    y1_->blockSignals(false);
    z1_->blockSignals(false);
    x2_->blockSignals(false);
    y2_->blockSignals(false);
    z2_->blockSignals(false);
}

void FloatMat3PropertyWidget::updateValue() {
    tgt::mat3 matrix;
    matrix.t00 = x0_->text().toDouble();
    matrix.t01 = y0_->text().toDouble();
    matrix.t02 = z0_->text().toDouble();
    matrix.t10 = x1_->text().toDouble();
    matrix.t11 = y1_->text().toDouble();
    matrix.t12 = z1_->text().toDouble();
    matrix.t20 = x2_->text().toDouble();
    matrix.t21 = y2_->text().toDouble();
    matrix.t22 = z2_->text().toDouble();
    prop_->set(matrix);
    emit valueModifiedByUser();
}

void FloatMat3PropertyWidget::identity() {
    x0_->setText("1");
    y0_->setText("0");
    z0_->setText("0");
    x1_->setText("0");
    y1_->setText("1");
    z1_->setText("0");
    x2_->setText("0");
    y2_->setText("0");
    z2_->setText("1");
}

void FloatMat3PropertyWidget::saveMatrix(const std::string& filename) const {
    tgt::mat3 matrix = prop_->get();
    SerializationHelper::save<tgt::mat3>(filename, "Matrix", matrix);
}

void FloatMat3PropertyWidget::loadMatrix(const std::string& filename) {
    tgt::mat3 matrix;
    SerializationHelper::load<tgt::mat3>(filename, "Matrix", matrix);
    prop_->set(matrix);
    emit valueModifiedByUser();
}

} //namespace voreen
