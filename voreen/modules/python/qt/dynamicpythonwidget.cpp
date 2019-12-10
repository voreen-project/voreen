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

#include "dynamicpythonwidget.h"
#include "voreen/qt/voreenapplicationqt.h"

#include "../pythonmodule.h"
#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QCheckBox>


namespace voreen {

const std::string DynamicPythonWidget::loggerCat_("voreen.DynamicPythonWidget");

DynamicPythonWidget::DynamicPythonWidget(QWidget* parent, DynamicPythonProcessor* pyProcessor)
    : QProcessorWidget(pyProcessor, parent)
    , plugin_(nullptr)
{
    tgtAssert(pyProcessor, "No DynamicPythonProcessor processor");

    setWindowTitle(QString::fromStdString(pyProcessor->getID()));
    resize(800, 480);
}

DynamicPythonWidget::~DynamicPythonWidget() {
    delete plugin_;
}

void DynamicPythonWidget::initialize() {
    QProcessorWidget::initialize();

    DynamicPythonProcessor* pyProcessor = dynamic_cast<DynamicPythonProcessor*>(processor_);

    plugin_ = new PythonPlugin(pyProcessor->getPythonProperty(), parentWidget());

    QGridLayout* layout = new QGridLayout();
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(plugin_);
    setLayout(layout);

    if(isVisible())
        show();

    initialized_ = true;
}

void DynamicPythonWidget::updateFromProcessor() {
    plugin_->updateFromProperty();
}

} //namespace voreen

