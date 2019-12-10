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

#include "voreen/qt/widgets/property/qpropertywidgetwithtoolwindow.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/voreenapplicationqt.h"

#include "voreen/core/datastructures/meta/windowstatemetadata.h"

#include <QAction>
#include <QApplication>
#include <QMainWindow>
#include <QDesktopWidget>

namespace {
    const std::string META_DATA_NAME = "ToolWindow";
}

namespace voreen {

QPropertyWidgetWithToolWindow::QPropertyWidgetWithToolWindow(Property* prop, QWidget* parent, bool showNameLabel, bool isToolWindowResizable)
    : QPropertyWidget(prop, parent, showNameLabel)
    , isToolWindowResizable_(isToolWindowResizable)
    , toolWindow_(0)
{
}

QPropertyWidgetWithToolWindow::~QPropertyWidgetWithToolWindow() {
    delete toolWindow_;
}

MetaDataBase* QPropertyWidgetWithToolWindow::getWidgetMetaData() const {

    // calculate position relative to mainwindow
    QPoint mainPindowPos(0,0);
    if (VoreenApplicationQt::qtApp() && VoreenApplicationQt::qtApp()->getMainWindow()) {
        mainPindowPos = VoreenApplicationQt::qtApp()->getMainWindow()->pos();
    }

    WindowStateMetaData* meta;
    if (toolWindow_ && isToolWindowResizable_) {
        meta = new WindowStateMetaData(
            toolWindow_->isVisible(),
            toolWindow_->pos().x() - mainPindowPos.x(),
            toolWindow_->pos().y() - mainPindowPos.y(),
            toolWindow_->width(),
            toolWindow_->height());
    }
    else if (toolWindow_ && !isToolWindowResizable_) {
        meta = new WindowStateMetaData(
            toolWindow_->isVisible(),
            toolWindow_->pos().x() - mainPindowPos.x(),
            toolWindow_->pos().y() - mainPindowPos.y());
    }
    else
        meta = new WindowStateMetaData(false);

    const_cast<QPropertyWidgetWithToolWindow*>(this)->getProperty()->getMetaDataContainer().addMetaData(META_DATA_NAME, meta);

    return meta;
}

void QPropertyWidgetWithToolWindow::createToolWindow(Qt::DockWidgetArea area, const QString& titlePostfix, const int& initialWidth, const int& initialHeight) {
    tgtAssert(!toolWindow_, "Tool window already instantiated");

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

    QString title;
    if (getProperty()->getOwner()) {
        title.append(getProperty()->getOwner()->getID().c_str());
        title.append(" - ");
    }
    title.append(QString::fromStdString(getProperty()->getGuiName()));
    title.append(titlePostfix);
    // replace line breaks by spaces
    title.replace("\r\n", " ");
    title.replace("\r", " ");
    title.replace("\n", " ");

    QMainWindow* mainWindow = VoreenApplicationQt::qtApp()->getMainWindow();
    toolWindow_ = new VoreenToolWindow(new QAction(title, mainWindow), mainWindow, createToolWindowWidget(), title, false);
    if (mainWindow && area != Qt::NoDockWidgetArea)
        mainWindow->addDockWidget(area, toolWindow_);
    customizeToolWindow();

    // set default size, might be overwritten by meta data
    if (initialWidth >= 0 && initialHeight >= 0)
        toolWindow_->resize(initialWidth, initialHeight);

    WindowStateMetaData* meta = dynamic_cast<WindowStateMetaData*>(getProperty()->getMetaDataContainer().getMetaData(META_DATA_NAME));
    // restore serialized geometry
    if (meta) {
        // compute position relative to mainwindow
        int xrel = meta->getX();
        int yrel = meta->getY();
        if (VoreenApplicationQt::qtApp() && VoreenApplicationQt::qtApp()->getMainWindow()) {
            QPoint mainPindowPos = VoreenApplicationQt::qtApp()->getMainWindow()->pos();
            xrel += mainPindowPos.x();
            yrel += mainPindowPos.y();
        }

        // check whether serialized left-top corner lies inside the available screen geometry
        QRect screenGeometry = QApplication::desktop()->availableGeometry(QPoint(xrel+25, yrel+25));
        if (screenGeometry.contains(QPoint(xrel+25, yrel+25))) {
            toolWindow_->move(xrel, yrel);
        }
        else {
            LWARNINGC("voreenqt.QPropertyWidgetWithEditorWindow",
                getPropertyGuiName() << " editor: Serialized position (" << meta->getX() << ", " << meta->getY() << ") " <<
                "outside visible desktop area. Ignoring.");
        }

        // size
        if (isToolWindowResizable_ && meta->getWidth() > 0 && meta->getHeight() > 0)
            toolWindow_->resize(meta->getWidth(), meta->getHeight());

        // visibility
        toolWindow_->setVisible(meta->getVisible());
    }

    QApplication::restoreOverrideCursor();
}


void QPropertyWidgetWithToolWindow::toggleToolWindow() {
    if(toolWindow_)
        toolWindow_->setVisible(!toolWindow_->isVisible());
}

bool QPropertyWidgetWithToolWindow::isToolWindowVisibleOnStartup() const {
    WindowStateMetaData* meta = dynamic_cast<WindowStateMetaData*>(
        const_cast<QPropertyWidgetWithToolWindow*>(this)->getProperty()->getMetaDataContainer().getMetaData(META_DATA_NAME));

    if (!meta)
        return false;

    return meta->getVisible();
}

} // namespace
