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

#include "modules/opencl/properties/openclproperty.h"
#include "modules/opencl/qt/openclpropertywidget.h"
#include "modules/opencl/qt/openclplugin.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/shaderproperty.h"

#include <QPushButton>

namespace voreen {

OpenCLPropertyWidget::OpenCLPropertyWidget(OpenCLProperty* prop, QWidget* parent)
    : QPropertyWidgetWithToolWindow(prop, parent)
    , plugin_(0)
    , property_(prop)
    , editBt_(new QPushButton(tr("edit")))
{

    if (isToolWindowVisibleOnStartup())
        createToolWindow(Qt::LeftDockWidgetArea, QString::fromStdString(" (original source: " + property_->get().programFilename_ + ")"), 700, 700);

    addWidget(editBt_);

    connect(editBt_, SIGNAL(clicked()), this, SLOT(setProperty()));
    connect(editBt_, SIGNAL(clicked()), this, SIGNAL(widgetChanged()));

    QFontInfo fontInfo(font());
    editBt_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
}

void OpenCLPropertyWidget::updateFromPropertySlot() {
    if (plugin_) {
        plugin_->updateFromProperty();
        plugin_->update();
    }
}

void OpenCLPropertyWidget::setProperty() {
    if (!disconnected_) {
        // lazy instantiation of shader editor window
        if (!toolWindow_) {
            createToolWindow(Qt::LeftDockWidgetArea, QString::fromStdString(" (original source: " + property_->get().programFilename_ + ")"), 700, 700);

            tgtAssert(toolWindow_, "OpenCL editor not instantiated");
        }

        if (toolWindow_->isVisible()) {
            //close widget
            toolWindow_->close();
        }
        else {
            //open Widget
            toolWindow_->showNormal();
        }
    }
}

void OpenCLPropertyWidget::disconnect() {
    disconnected_ = true;
    if (plugin_)
        plugin_->disconnect();
}

QWidget* OpenCLPropertyWidget::createToolWindowWidget() {
    plugin_ = new OpenCLPlugin(property_, parentWidget());
    plugin_->createWidgets();
    plugin_->createConnections();
    //connect(plugin_, SIGNAL(modified()), this, SIGNAL(modified()));

    return plugin_;
}

void OpenCLPropertyWidget::customizeToolWindow() {
}

} // namespace voreen
