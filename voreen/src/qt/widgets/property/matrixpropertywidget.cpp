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

#include "voreen/core/properties/property.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"

#include "voreen/qt/widgets/property/matrixpropertywidget.h"

#include <QCheckBox>
#include <QDialog>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QFileDialog>
#include <QUrl>
#include <QDesktopServices>
#include <QErrorMessage>

namespace voreen {

MatrixPropertyWidget::MatrixPropertyWidget(Property* prop, int numColumns, QWidget* parent)
    : QPropertyWidget(prop, parent)
{
    editButton_ = new QPushButton(tr("Edit"), this);
    updateButton_ = new QPushButton(tr("Update"), this);
    identityButton_ = new QPushButton(tr("Identity"), this);
    loadButton_ = new QPushButton(tr("Load"), this);
    saveButton_ = new QPushButton(tr("Save"), this);
    autoUpdate_ = new QCheckBox("Auto Update", this);
    QFontInfo fontInfo(font());
    editButton_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
    layout_->addWidget(editButton_);
    matrixDialog_ = new QDialog(this);
    mainLayout_ = new QGridLayout(matrixDialog_);

    connect(editButton_, SIGNAL(clicked()), this, SLOT(executeDialog()));
    connect(updateButton_, SIGNAL(clicked()), this, SLOT(updateValue()));
    connect(identityButton_, SIGNAL(clicked()), this, SLOT(identity()));
    connect(saveButton_, SIGNAL(clicked()), this, SLOT(save()));
    connect(loadButton_, SIGNAL(clicked()), this, SLOT(load()));

    QHBoxLayout* buttonLayout = new QHBoxLayout;
    buttonLayout->addWidget(autoUpdate_);
    buttonLayout->addWidget(updateButton_);
    buttonLayout->addWidget(identityButton_);
    buttonLayout->addWidget(loadButton_);
    buttonLayout->addWidget(saveButton_);
    mainLayout_->addLayout(buttonLayout, 0, 0, 1, numColumns);
}

void MatrixPropertyWidget::executeDialog() {
    matrixDialog_->move(QPoint(QCursor::pos().x() - 500, QCursor::pos().y()));
    matrixDialog_->show();
}

void MatrixPropertyWidget::save() {
    QFileDialog fileDialog(this, tr("Save Matrix As..."), VoreenApplication::app()->getUserDataPath().c_str());
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setConfirmOverwrite(true);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QStringList filters;
    filters << (QString("Voreen Matrix (*.") + QString::fromStdString(getExtension()) + ")");
    fileDialog.setNameFilters(filters);
    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    if (fileDialog.exec()) {
        QString name = fileDialog.selectedFiles().at(0);
        std::string filename = name.toStdString();
        if (!name.endsWith(QString(".") + QString::fromStdString(getExtension())))
            filename += "." + getExtension();

        try {
            saveMatrix(filename);
        }
        catch (SerializationException& e) {
            QErrorMessage* errorMessageDialog = new QErrorMessage(this);
            errorMessageDialog->setWindowTitle(tr("Matrix serialization failed"));
            errorMessageDialog->showMessage(e.what());
        }
    }
}

void MatrixPropertyWidget::load() {
    QFileDialog fileDialog(this, tr("Load Matrix"), VoreenApplication::app()->getUserDataPath().c_str());
    fileDialog.setFileMode(QFileDialog::ExistingFile);
    fileDialog.setAcceptMode(QFileDialog::AcceptOpen);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QStringList filters;
    filters << QString("Voreen Matrix (*.") + QString::fromStdString(getExtension()) + ")";
    fileDialog.setNameFilters(filters);
    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    if (fileDialog.exec()) {
        std::string filename = fileDialog.selectedFiles().at(0).toStdString();

        try {
            loadMatrix(filename);
        }
        catch (SerializationException& e) {
            QErrorMessage* errorMessageDialog = new QErrorMessage(this);
            errorMessageDialog->setWindowTitle(tr("Matrix deserialization failed"));
            errorMessageDialog->showMessage(e.what());
        }
    }
}

void MatrixPropertyWidget::autoUpdate() {
    if(autoUpdate_->isChecked())
        updateValue();
}

} //namespace voreen
