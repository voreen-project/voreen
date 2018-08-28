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

#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/screenshotmenuentity.h"
#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/networkscreenshotplugin.h"
#include "voreen/qt/networkeditor/networkeditor.h"
#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/canvasrendererscreenshotplugin.h"
#include "modules/core/processors/output/canvasrenderer.h"

#include "modules/core/qt/processor/canvasrendererwidget.h"

#include <QAction>
#include <QMenu>

namespace voreen {

ScreenshotMenuEntity::ScreenshotMenuEntity()
    : VoreenQtMenuEntity()
    , networkScreenshotPlugin_(0)
    , currentWorkspace_(0)
{
    WsHndlr.registerWorkspaceUsingWidget(this);
}

ScreenshotMenuEntity::~ScreenshotMenuEntity() {
    WsHndlr.unregisterWorkspaceUsingWidget(this);
    delete networkScreenshotPlugin_;
    for(auto it = canvasRendererScreenshotPluginMap_.begin(); it != canvasRendererScreenshotPluginMap_.end(); it++) {
        delete it->second;
    }
}

void ScreenshotMenuEntity::setMainWindow(VoreenQtMainWindow* mainWindow) {
    tgtAssert(mainWindow, "null pointer passed");
    mainWindow_ = mainWindow;
    //remove network screenshot in application mode
    connect(mainWindow_,SIGNAL(guiModeUpdated(VoreenQtMainWindow::GuiMode)),this,SLOT(adjustScreenshotMenuSlot(VoreenQtMainWindow::GuiMode)));
}

    //------------------------------
    //  Menu and Main-Functions
    //------------------------------
QAction* ScreenshotMenuEntity::createMenuAction() {
    QAction* action = new QAction(getIcon(),QWidget::tr(getName().c_str()),0);
    action->setShortcut(QWidget::tr(getShortCut().c_str()));
    action->setData(0);
    QMenu* screenshotMenu = new QMenu();
    action->setMenu(screenshotMenu);
    action->setEnabled(false);
    connect(action, SIGNAL(triggered()), this, SLOT(screenshotActionTriggeredSlot()));
    return action;
}

void ScreenshotMenuEntity::adjustScreenshotMenuSlot(VoreenQtMainWindow::GuiMode mode) {
#ifdef VRN_MODULE_DEVIL
    if (!menuAction_ || !toolbarAction_)
        return;

    //clear old actions
    foreach(QAction* act, menuAction_->menu()->actions()) {
        delete act;
    }
    menuAction_->menu()->clear();
    toolbarAction_->menu()->clear();

    if (!mainWindow_->getNetworkEditor()->getProcessorNetwork()) {
        menuAction_->setEnabled(false);
        toolbarAction_->setEnabled(false);
        return;
    }

    menuAction_->setEnabled(true);
    toolbarAction_->setEnabled(true);
    std::vector<CanvasRenderer*> canvasRenderers = mainWindow_->getNetworkEditor()->getProcessorNetwork()->getProcessorsByType<CanvasRenderer>();

    for (size_t i=0; i < canvasRenderers.size(); ++i) {
        QAction* menuAction = new QAction(QString::fromStdString(canvasRenderers[i]->getID()), this);
        // store index of canvasrenderer in action for identification
        // TODO: replace by name, when processors' names are unique
        menuAction->setData((int)i);
        menuAction_->menu()->addAction(menuAction);
        toolbarAction_->menu()->addAction(menuAction);
        connect(menuAction, SIGNAL(triggered()), this, SLOT(screenshotActionTriggeredSlot()));
    }

    // add network screenshot functionality
    if(mode == VoreenQtMainWindow::MODE_NETWORK) {
        QAction* networkAction = new QAction(tr("Network Graph"), this);
        networkAction->setData(-1);
        menuAction_->menu()->addSeparator();
        toolbarAction_->menu()->addSeparator();
        menuAction_->menu()->addAction(networkAction);
        toolbarAction_->menu()->addAction(networkAction);
        connect(networkAction, SIGNAL(triggered()), this, SLOT(screenshotActionTriggeredSlot()));
    } else { //clean up
        delete networkScreenshotPlugin_;
        networkScreenshotPlugin_ = 0;
    }
#endif
}

void ScreenshotMenuEntity::screenshotActionTriggeredSlot() {
    if (dynamic_cast<QAction*>(QObject::sender())) {
        int rendererIndex = static_cast<QAction*>(QObject::sender())->data().toInt();
        std::vector<CanvasRenderer*> canvasRenderers
            = mainWindow_->getNetworkEditor()->getProcessorNetwork()->getProcessorsByType<CanvasRenderer>();

        if (rendererIndex >= 0 && (size_t)rendererIndex < canvasRenderers.size()) {
            auto it = canvasRendererScreenshotPluginMap_.find(canvasRenderers[rendererIndex]);
            if(it == canvasRendererScreenshotPluginMap_.end()) {
                it = canvasRendererScreenshotPluginMap_.insert(std::make_pair(canvasRenderers[rendererIndex],
                                                                              new CanvasRendererScreenshotPlugin(mainWindow_, currentWorkspace_, canvasRenderers[rendererIndex])
                                                                              )).first;
                it->second->resolutionComboBoxChangedSlot(0); //initialize
            }
            it->second->show();
        }
        else if (rendererIndex == -1) {
            // screenshot of network graph
            if(!networkScreenshotPlugin_) {
                networkScreenshotPlugin_ = new NetworkScreenshotPlugin(mainWindow_, currentWorkspace_, mainWindow_->getNetworkEditor());
                networkScreenshotPlugin_->resolutionComboBoxChangedSlot(0); //initialize
            }
            networkScreenshotPlugin_->show();
        }
    }
}

    //------------------------------
    //  UsesWorkspace
    //------------------------------
void ScreenshotMenuEntity::setWorkspace(Workspace* workspace) {
    //handle old workspace
    for(auto it = canvasRendererScreenshotPluginMap_.begin(); it != canvasRendererScreenshotPluginMap_.end(); it++) {
        delete it->second;
    }
    canvasRendererScreenshotPluginMap_.clear();

    //handle new workspace
    if(networkScreenshotPlugin_)
        networkScreenshotPlugin_->updateWorkspace(workspace);

    //adjust menu
    adjustScreenshotMenuSlot(mainWindow_->getCurrentGuiMode());

    //update current workspace and observe network
    currentWorkspace_ = workspace;
    if(workspace && workspace->getProcessorNetwork())
        workspace->getProcessorNetwork()->addObserver(this);
}

    //------------------------------
    //  ProcessorNetworkObserver
    //------------------------------
void ScreenshotMenuEntity::processorAdded(const Processor* processor) {
    if(dynamic_cast<const CanvasRenderer*>(processor)) {
        //update menu
        adjustScreenshotMenuSlot(mainWindow_->getCurrentGuiMode());
    }
}

void ScreenshotMenuEntity::processorRemoved(const Processor* processor) {
    if(const CanvasRenderer* cr = dynamic_cast<const CanvasRenderer*>(processor)) {
        //update menu
        adjustScreenshotMenuSlot(mainWindow_->getCurrentGuiMode());
        //remove meta and plugin
        std::map<const CanvasRenderer*,CanvasRendererScreenshotPlugin*>::iterator it = canvasRendererScreenshotPluginMap_.find(cr);
        if(it != canvasRendererScreenshotPluginMap_.end()) { //plugin exists
            it->second->removeMetaData();
            delete it->second;
            canvasRendererScreenshotPluginMap_.erase(it);
        }
    }
}
    /** Update plugin name and meta data, if a CanvasRenderer has been renamed. */
void ScreenshotMenuEntity::processorRenamed(const Processor* processor, const std::string& prevName) {
    if(const CanvasRenderer* cr = dynamic_cast<const CanvasRenderer*>(processor)) {
        //update menu
        adjustScreenshotMenuSlot(mainWindow_->getCurrentGuiMode());
        //rename meta and plugin
        std::map<const CanvasRenderer*,CanvasRendererScreenshotPlugin*>::iterator it = canvasRendererScreenshotPluginMap_.find(cr);
        if(it != canvasRendererScreenshotPluginMap_.end()) { //plugin exists
            it->second->canvasRendererHasBeenRenamed();
        }
    }
}
} //namespace
