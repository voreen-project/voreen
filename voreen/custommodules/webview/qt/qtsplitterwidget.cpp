/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "qtsplitterwidget.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/processornetwork.h"

#include "voreen/qt/voreenapplicationqt.h"

#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QUrl>


namespace voreen {

const std::string QtSplitterWidget::loggerCat_("voreen.QtSplitterWidget");

QtSplitterWidget::QtSplitterWidget(QWidget* parent, QtSplitter* processor)
    : QProcessorWidget(processor, parent)
    , processor_(processor)
{
    tgtAssert(processor, "No QtSplitter processor");

    setWindowTitle(QString::fromStdString(processor->getID()));
    resize(800, 480);
}

void QtSplitterWidget::initialize() {
    QProcessorWidget::initialize();

    splitter_ = new QSplitter();
    QGridLayout* layout = new QGridLayout();
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(splitter_);
    setLayout(layout);

    connect(splitter_, SIGNAL(splitterMoved(int, int )), this, SLOT(splitterMoved()));

    if(isVisible())
        show();

    initialized_ = true;
}

void QtSplitterWidget::updateFromProcessor() {

    auto* app = VoreenApplication::app();
    if(!app) {
        LERROR("No GUI application");
        return;
    }

    // First update size. At this point, there still is a match between processor and widget.
    if(processor_->getSizes().size() == static_cast<size_t>(splitter_->count())) {
        splitter_->setSizes(QList<int>::fromStdList(processor_->getSizes()));
    }

    std::set<QWidget*> orphanedWidgets;
    for(int i=0; i<splitter_->count(); i++) {
        orphanedWidgets.insert(splitter_->widget(i));
    }

    auto* network = app->getNetworkEvaluator()->getProcessorNetwork();

    auto instances = processor_->getInstances();
    int numInstances = static_cast<int>(instances.size());
    for(int i=0; i<numInstances; i++) {

        auto* widget = network->getProcessor(instances[i])->getProcessorWidget();
        if(!widget) {
            LWARNING("No widget");
            continue;
        }

        auto* qwidget = dynamic_cast<QWidget*>(widget);
        tgtAssert(qwidget, "no QWidget");

        // We keep this widget.
        if(qwidget->parentWidget() == splitter_) {
            orphanedWidgets.erase(qwidget);
        }

        splitter_->insertWidget(i, qwidget);
    }

    // Reset parent of orphaned widgets.
    tgtAssert(VoreenApplicationQt::qtApp(), "No VoreenApplicationQt");
    for(auto* widget : orphanedWidgets) {
        // Settings hereby need to match the original parenting of the QProcessorWidget (see respective constructor).
        widget->setParent(VoreenApplicationQt::qtApp()->getMainWindow(), Qt::Tool);
        widget->show();
    }

    // Set / Update orientation.
    switch (processor_->getOrientation()) {
        case QtSplitter::Orientation::HORIZONTAL:
            splitter_->setOrientation(Qt::Horizontal);
            break;
        case QtSplitter::Orientation::VERTICAL:
            splitter_->setOrientation(Qt::Vertical);
            break;
        default:
            tgtAssert(false, "Unhandled orientation");
            break;
    }

    update();
}

void QtSplitterWidget::splitterMoved() {
    processor_->setSizes(splitter_->sizes().toStdList());
}

} //namespace voreen

