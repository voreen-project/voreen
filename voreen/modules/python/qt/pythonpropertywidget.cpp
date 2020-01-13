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

#include "modules/python/properties/pythonproperty.h"
#include "modules/python/qt/pythonpropertywidget.h"
#include "modules/python/qt/pythonplugin.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/core/processors/processor.h"

#include <QPushButton>

namespace voreen {

PythonPropertyWidget::PythonPropertyWidget(PythonProperty* prop, QWidget* parent)
    : QPropertyWidgetWithToolWindow(prop, parent)
    , plugin_(nullptr)
    , property_(prop)
    , editBt_(new QPushButton(tr("edit")))
{

    if (isToolWindowVisibleOnStartup())
        createToolWindow(Qt::LeftDockWidgetArea, QString::fromStdString(" (original source: " + property_->get().getFilename() + ")"), 700, 700);

    addWidget(editBt_);

    connect(editBt_, SIGNAL(clicked()), this, SLOT(setProperty()));
    connect(editBt_, SIGNAL(clicked()), this, SIGNAL(widgetChanged()));

    QFontInfo fontInfo(font());
    editBt_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
}

void PythonPropertyWidget::updateFromPropertySlot() {
    if (plugin_) {
        plugin_->updateFromProperty();
        plugin_->update();
    }
}

void PythonPropertyWidget::setProperty() {
    if (!disconnected_) {
        // lazy instantiation of shader editor window
        if (!toolWindow_) {
            createToolWindow(Qt::LeftDockWidgetArea, QString::fromStdString(" (original source: " + property_->get().getFilename() + ")"), 700, 700);

            tgtAssert(toolWindow_, "Python editor not instantiated");
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

void PythonPropertyWidget::disconnect() {
    disconnected_ = true;
    if (plugin_)
        plugin_->disconnect();
}

QWidget* PythonPropertyWidget::createToolWindowWidget() {
    plugin_ = new PythonPlugin(property_, parentWidget());
    return plugin_;
}

void PythonPropertyWidget::customizeToolWindow() {
}

} // namespace voreen
