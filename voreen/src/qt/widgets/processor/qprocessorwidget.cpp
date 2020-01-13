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

#include "tgt/vector.h"
#include "voreen/qt/widgets/processor/qprocessorwidget.h"
#include "voreen/qt/voreenapplicationqt.h"
#include "voreen/core/processors/processor.h"

#include <QDialog>
#include <QApplication>
#include <QMainWindow>
#include <QDesktopWidget>
#include <QMoveEvent>

namespace voreen {

QProcessorWidget::QProcessorWidget(Processor* processor, QWidget* parent)
    : QWidget(parent, Qt::Tool)
    , ProcessorWidget(processor)
{
    tgtAssert(processor, "No processor");
    setWindowTitle(QString::fromStdString(processor->getID()));
    QWidget::setVisible(false);
}

void QProcessorWidget::setVisible(bool visible) {
    if (isInitialized()) {
        if (visible) {
            initialized_ = false;
            restoreGeometryFromMeta();
            initialized_ = true;
        }
        else {
            saveGeometryToMeta();
        }
    }

    QWidget::setVisible(visible);
}

bool QProcessorWidget::isVisible() const {
    return QWidget::isVisible();
}

void QProcessorWidget::setSize(int height, int width) {
    QWidget::resize(height, width);
}

tgt::ivec2 QProcessorWidget::getSize() const {
    return tgt::ivec2(QWidget::width(), QWidget::height());
}

void QProcessorWidget::setPosition(int x, int y) {

    // compute position relative to mainwindow
    int xrel = x;
    int yrel = y;
    if (VoreenApplicationQt::qtApp() && VoreenApplicationQt::qtApp()->getMainWindow()) {
        QPoint mainPindowPos = VoreenApplicationQt::qtApp()->getMainWindow()->pos();
        xrel += mainPindowPos.x();
        yrel += mainPindowPos.y();
    }

    // check whether top-left corner lies inside the available screen geometry
    QRect screenGeometry = QApplication::desktop()->availableGeometry(QPoint(xrel+25,yrel+25));
    if (screenGeometry.contains(QPoint(xrel+25,yrel+25))) {
        QWidget::move(xrel, yrel);
    }
    else {
        //to remove the warning in fullscreen mode
        if(!isFullScreen() && processor_->isInitialized())
        LWARNINGC("voreenqt.QProcessorWidget",
            "setPosition(" << x << ", " << y << "): '" << windowTitle().toStdString() <<
            "' would be placed outside visible desktop area. Ignoring.");
    }

}

tgt::ivec2 QProcessorWidget::getPosition() const {
    // calculate position relative to mainwindow
    if (VoreenApplicationQt::qtApp() && VoreenApplicationQt::qtApp()->getMainWindow()) {
       QPoint mainPindowPos = VoreenApplicationQt::qtApp()->getMainWindow()->pos();
       return tgt::ivec2(QWidget::x() - mainPindowPos.x(), QWidget::y() - mainPindowPos.y());
    }
    else
        return tgt::ivec2(QWidget::x(), QWidget::y());
}

void QProcessorWidget::processorNameChanged() {
    setWindowTitle(QString::fromStdString(processor_->getGuiName()));
}

void QProcessorWidget::moveEvent(QMoveEvent* eve) {
    QWidget::moveEvent(eve);
    ProcessorWidget::onMove();
}

void QProcessorWidget::resizeEvent(QResizeEvent* eve) {
    QWidget::resizeEvent(eve);
    ProcessorWidget::onResize();
}

void QProcessorWidget::showEvent(QShowEvent* eve) {
    QWidget::showEvent(eve);
    ProcessorWidget::onShow();
}

void QProcessorWidget::hideEvent(QHideEvent* eve) {
    QWidget::hideEvent(eve);
    ProcessorWidget::onHide();
}

} // namespace voreen
