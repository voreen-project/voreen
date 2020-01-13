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

#include "voreen/qt/widgets/networkmodeconfigdialog.h"

#include "voreen/qt/networkeditor/networkeditor.h"

#include "voreen/qt/widgets/property/propertyownerwidget.h"

#include "voreen/qt/voreenapplicationqt.h"

#include <QApplication>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QScrollArea>
#include <QGridLayout>
#include <QPushButton>
#include <QMessageBox>

namespace voreen {

NetworkModeConfigDialog::NetworkModeConfigDialog(NetworkEditor* nwe, QWidget* parent)
    : QDialog(parent, static_cast<Qt::WindowFlags>(Qt::Tool | Qt::Window))
    ,nwe_(nwe)
{
    setWindowTitle(tr("Network Mode Configuration"));

    QVBoxLayout* mainLayout = new QVBoxLayout();
    setLayout(mainLayout);

    // scroll area
    QScrollArea* scrollArea = new QScrollArea(this);
    scrollArea->setWidgetResizable(true);
    scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scrollArea->setContentsMargins(0, 0, 0, 0);
    mainLayout->addWidget(scrollArea);

    // container widget for PropertyOwnerWidgets
    QWidget* containerWidget = new QWidget();
    scrollArea->setWidget(containerWidget);
    QVBoxLayout* widgetLayout = new QVBoxLayout();
    containerWidget->setLayout(widgetLayout);

    // application widget
    PropertyOwnerWidget* appPow = new PropertyOwnerWidget(nwe_, this, "Network Editor", true, false);
    widgetLayout->addWidget(appPow);

    widgetLayout->addStretch();

    // button row
    mainLayout->addSpacerItem(new QSpacerItem(1, 4));
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    resetButton_ = new QPushButton("Reset");
    buttonLayout->addWidget(resetButton_);
    buttonLayout->addStretch();
    closeButton_ = new QPushButton("Close");
    buttonLayout->addWidget(closeButton_);
    mainLayout->addLayout(buttonLayout);

    connect(resetButton_, SIGNAL(clicked()), this, SLOT(resetSettings()));
    connect(closeButton_, SIGNAL(clicked()), this, SIGNAL(closeSettings()));
}

void NetworkModeConfigDialog::resetSettings() {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle(QApplication::tr("Reset Settings"));
    msgBox.setIcon(QMessageBox::Question);
    msgBox.setText(QApplication::tr("This will reset all network mode settings to their default values."));
    msgBox.setInformativeText(QApplication::tr("Do you want to proceed?"));
    msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
    msgBox.setDefaultButton(QMessageBox::Cancel);
    int ret = msgBox.exec();
    if (ret == QMessageBox::Ok) {
        nwe_->resetAllProperties();
        //HACK: forces items to repaint
        QApplication::processEvents();
        QApplication::setActiveWindow(nwe_);
        QApplication::setActiveWindow(this);
    }
}

} // namespace
