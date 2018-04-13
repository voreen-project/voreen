/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/qt/mainwindow/menuentities/voreenqtmenuentity.h"

#include "voreen/core/network/networkevaluator.h"
#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

#include <QAction>
#include <QWidget>

namespace voreen {

VoreenQtMenuEntity::VoreenQtMenuEntity()
    : initialized_(false)
    , menuAction_(0)
    , toolbarAction_(0)
    , toolWindow_(0)
    , mainWindow_(0)
{}


VoreenQtMenuEntity::~VoreenQtMenuEntity() {
    delete menuAction_;
    menuAction_ = 0;
    delete toolbarAction_;
    toolbarAction_ = 0;
    delete toolWindow_;
    toolWindow_ = 0;
}

bool VoreenQtMenuEntity::isInitialized() const {
    return initialized_;
}

QAction* VoreenQtMenuEntity::getMenuAction() const {
    tgtAssert(initialized_,"getMenuBarAction called on not initialized menu entity!");
    return menuAction_;
}

QAction* VoreenQtMenuEntity::getToolBarAction() const {
    tgtAssert(initialized_,"getToolBarAction called on not initialized menu entity!");
    return toolbarAction_;
}

VoreenToolWindow* VoreenQtMenuEntity::getToolWindow() const {
    tgtAssert(initialized_,"getToolWidget called on not initialized menu entity!");
    return toolWindow_;
}

void VoreenQtMenuEntity::initialize() {
    tgtAssert(!isInitialized(),"initialize() is called on already initialized entity!");

    if(!mainWindow_)
        throw VoreenException("No main window assigned before initialization!");

    QWidget* widget = createWidget();
    //if no widget is needed, no tool window will be created
    if(widget) {
        toolWindow_ = createToolWindow(widget);
        toolWindow_->updateTitle(QString::fromStdString(getName()),getIcon());
    }

    menuAction_ = createMenuAction();
    toolbarAction_ = createToolbarAction();

    initialized_ = true;
}

void VoreenQtMenuEntity::deinitialize() {
    initialized_ = false;
}

QAction* VoreenQtMenuEntity::createMenuAction() {
    QAction* action = new QAction(getIcon(),QWidget::tr(getName().c_str()),0);
    action->setShortcut(QWidget::tr(getShortCut().c_str()));
    action->setCheckable(true);
    action->setChecked(false);
    if(toolWindow_) {
        QObject::connect(action, SIGNAL(triggered(bool)), toolWindow_, SLOT(setVisible(bool)));
        QObject::connect(toolWindow_->toggleViewAction(), SIGNAL(toggled(bool)), action, SLOT(setChecked(bool)));
    }
    return action;
}

QAction* VoreenQtMenuEntity::createToolbarAction() {
    QAction* act = createMenuAction();
    act->setShortcut(QWidget::tr("")); //remove shortcut to prevent ambiguous shortcuts
    return act;
}

VoreenToolWindow* VoreenQtMenuEntity::createToolWindow(QWidget* widget) {
    if(!widget)
        return 0;

    bool dockable = (getAllowedDockWidgetAreas() != Qt::NoDockWidgetArea);

    VoreenToolWindow* window = new VoreenToolWindow(0, mainWindow_, widget, QWidget::tr(getName().c_str()), dockable);
    window->setAllowedAreas(getAllowedDockWidgetAreas());

    if (dockable) {
        if (getInitialDockWidgetArea() == Qt::NoDockWidgetArea) {
            mainWindow_->addDockWidget(Qt::LeftDockWidgetArea, window);
            window->setFloating(true);
        }
        else {
            mainWindow_->addDockWidget(getInitialDockWidgetArea(), window);
        }
    }

    return window;
}

void VoreenQtMenuEntity::setMainWindow(VoreenQtMainWindow* mainWindow) {
    tgtAssert(mainWindow, "null pointer passed");
    mainWindow_ = mainWindow;
}

}
