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

#include "dynamicglslwidget.h"
#include "voreen/qt/voreenapplicationqt.h"

#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QCheckBox>


namespace voreen {

const std::string DynamicGLSLWidget::loggerCat_("voreen.DynamicGLSLWidget");

DynamicGLSLWidget::DynamicGLSLWidget(QWidget* parent, DynamicGLSLProcessor* glslProcessor)
    : QProcessorWidget(glslProcessor, parent)
{
    tgtAssert(glslProcessor, "No DynamicGLSLProcessor processor");

    setWindowTitle(QString::fromStdString(glslProcessor->getGuiName()));
    resize(800, 480);
}

DynamicGLSLWidget::~DynamicGLSLWidget() {
    delete plugin_;
}

void DynamicGLSLWidget::initialize() {
    QProcessorWidget::initialize();

    DynamicGLSLProcessor* glslProcessor = dynamic_cast<DynamicGLSLProcessor*>(processor_);

    plugin_ = new ShaderPlugin(glslProcessor->getShader(), parentWidget());
    plugin_->createWidgets();
    plugin_->createConnections();
    //connect(plugin_, SIGNAL(modified()), this, SIGNAL(modified()));

    QGridLayout* layout = new QGridLayout();
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(plugin_);
    setLayout(layout);

    //show();

    initialized_ = true;
}

void DynamicGLSLWidget::updateFromProcessor() {
    plugin_->updateFromProperty();
}

} //namespace voreen

