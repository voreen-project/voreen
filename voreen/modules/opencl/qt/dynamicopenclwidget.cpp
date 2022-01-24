/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "dynamicopenclwidget.h"
#include "openclpropertywidget.h"
#include "voreen/qt/voreenapplicationqt.h"

#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QCheckBox>


namespace voreen {

const std::string DynamicOpenCLWidget::loggerCat_("voreen.DynamicOpenCLWidget");

DynamicOpenCLWidget::DynamicOpenCLWidget(QWidget* parent, DynamicCLProcessor* clProcessor)
    : QProcessorWidget(clProcessor, parent)
    , plugin_(nullptr)
{
    tgtAssert(clProcessor, "No DynamicCLProcessor processor");

    setWindowTitle(QString::fromStdString(clProcessor->getID()));
    resize(800, 480);
}

DynamicOpenCLWidget::~DynamicOpenCLWidget() {
    delete plugin_;
}

void DynamicOpenCLWidget::initialize() {
    QProcessorWidget::initialize();

    DynamicCLProcessor* clProcessor = dynamic_cast<DynamicCLProcessor*>(processor_);

    plugin_ = new OpenCLPlugin(clProcessor->getOpenCLProperty());
    //OpenCLPropertyWidget* clWidget = (OpenCLPropertyWidget*)(*(clProcessor->getOpenCLProperty()->getPropertyWidgets().begin()));
    //plugin_ = clWidget->getPlugin();
    plugin_->createWidgets();
    plugin_->createConnections();

    QGridLayout* layout = new QGridLayout();
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(plugin_);
    setLayout(layout);

    if(isVisible())
        show();

    initialized_ = true;
}

void DynamicOpenCLWidget::updateFromProcessor() {
    plugin_->updateFromProperty();
}

} //namespace voreen

