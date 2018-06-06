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

#include "voreenbiologymainwindow.h"
#include "voreenbiologysplashscreen.h"
#include "voreenbiologystartupwizard.h"
#include "voreenbiologyapplication.h"
#include "voreenbiologyaboutbox.h"

#include "voreen/qt/mainwindow/menuentities/voreenqtmenuentity.h"

#include "voreen/qt/networkeditor/networkeditor.h"
#include "voreen/qt/widgets/processorlistwidget.h"
#include "voreen/qt/widgets/propertylistwidget.h"
#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/widgets/volumeviewer.h"
#include "voreen/qt/widgets/consoleplugin.h"

#include "voreen/qt/helpbrowser.h"

namespace voreen {

namespace {

const int MAX_RECENT_FILES = 8;

// Version number of restoring state of the main window.
// Increase when incompatible changes happen.
const int WINDOW_STATE_VERSION = 15;  // V4.0

} // namespace

////////// VoreenVEMainWindow ////////////////////////////////////////////////////////////

namespace {

// Add glass highlight effect to standard menubar
class FancyMenuBar : public QMenuBar {
protected:
    void paintEvent(QPaintEvent* event) {
        QMenuBar::paintEvent(event);

        // draw semi-transparent glass highlight over upper half of menubar
        QPainter painter(this);
        painter.setBrush(QColor(255, 255, 255, 76));
        painter.setPen(Qt::NoPen);
        painter.drawRect(0, 0, rect().width(), rect().height() / 2);
    }
};

} // namespace

const std::string VoreenBiologyMainWindow::loggerCat_("voreenbiology.VoreenBiologyMainWindow");

VoreenBiologyMainWindow::VoreenBiologyMainWindow(const std::string& workspace, bool noInitialWorkspace, bool resetSettings)
    : VoreenQtMainWindow("Voreen Biology", workspace, noInitialWorkspace, resetSettings)
{
    //general settings
    setWindowIcon(QIcon(":/voreenbiology/icons/logo.png"));
}

VoreenBiologyMainWindow::~VoreenBiologyMainWindow() {
}

//**********************************************************************************************
//      switch bewteen network and application mode
//**********************************************************************************************
void VoreenBiologyMainWindow::setGuiMode(GuiMode guiMode) {
    switch(guiMode) {
    case MODE_APPLICATION:
        //if we are in network mode, save all window settings
        if (currentGuiMode_ == MODE_NETWORK) {
            networkModeState_ = saveState(WINDOW_STATE_VERSION);
            networkEditorWindowState_ = networkEditorWindow_->saveGeometry();
        }

        //restore all windows
        setUpdatesEnabled(false); //why?
        // hide all first to prevent some flicker
        networkEditorWindow_->hide();
        networkEditorWidget_->setVisible(false);
        if (getToolWindow(processorListWidget_))
            getToolWindow(processorListWidget_)->hide();
        setUpdatesEnabled(true);

        if (!restoreState(applicationModeState_, WINDOW_STATE_VERSION)) {
            if (getToolWindow(processorListWidget_) && getToolWindow(processorListWidget_)->isEnabled()) {
                if (getToolWindow(propertyListWidget_))
                    getToolWindow(propertyListWidget_)->show();
            }
            if (getToolWindow(volumeViewer_) && getToolWindow(volumeViewer_)->isEnabled())
                getToolWindow(volumeViewer_)->hide();
            if (getToolWindow(consolePlugin_))
                getToolWindow(consolePlugin_)->hide();
            if (getToolWindow(wsdWidget_))
                getToolWindow(wsdWidget_)->show();
        }

        if (getToolWindow(propertyListWidget_))
            getToolWindow(propertyListWidget_)->show();

       
        // resize canvas after gui has been settled
        adjustCanvasWidgets(MODE_APPLICATION);
        //update property list widget
        setUpdatesEnabled(false);
        if(propertyListWidget_) propertyListWidget_->setWidgetMode(PropertyListWidget::APPLICATION);
        setUpdatesEnabled(true);
        break;
    case MODE_NETWORK:
        tgtAssert(false, "MODE_NETWORK is not supported by Voreen Biology");
        break;
    default:
        tgtAssert(false, "unknown gui mode state!");
        break;
    }
    //update menu entities
    adjustMenuEntitiesAndToolBarsVisibility(guiMode);
    //set current gui mode
    currentGuiMode_ = guiMode;
}

//**********************************************************************************************
//      create GUI elements (Menu and Toolrbar)
//**********************************************************************************************
void VoreenBiologyMainWindow::createQMenus() {
    VoreenQtMainWindow::createQMenus();
    //remove view tab in main window
    mainMenu_->removeAction(viewMenu_->menuAction());

    //
    // Help menu
    //
    QAction* helpFirstStepsAct = new QAction(QIcon(":/qt/icons/help.png"), tr("&Getting Started..."), this);
    helpFirstStepsAct->setShortcut(tr("F1"));
    connect(helpFirstStepsAct, SIGNAL(triggered()), this, SLOT(helpFirstSteps()));
    helpMenu_->addAction(helpFirstStepsAct);
    helpMenu_->addSeparator();
    // Add some web links to the menu. Use the redirects (in the "go" directory) to be
    // independent of website reorganization.
    QAction* websiteAct = new QAction(tr("Voreen Website..."), this);
    websiteAct->setData(tr("http://voreen.uni-muenster.de"));
    connect(websiteAct, SIGNAL(triggered()), this, SLOT(helpWebsite()));
    helpMenu_->addAction(websiteAct);
    helpMenu_->addSeparator();
    QAction* aboutAction = new QAction(QIcon(":/qt/icons/about.png"), tr("&About VoreenBiology..."), this);
    connect(aboutAction, SIGNAL(triggered()), this, SLOT(helpAbout()));
    helpMenu_->addAction(aboutAction);
}

//**********************************************************************************************
//      help menu functions
//**********************************************************************************************
void VoreenBiologyMainWindow::helpFirstSteps() {
    //QDesktopServices::openUrl(QUrl("http://voreen.uni-muenster.de/?q=user-interface"));
    QDesktopServices::openUrl(QUrl(("file:///" + tgt::FileSystem::cleanupPath(VoreenApplication::app()->getBasePath() + "/doc/voreenbiology/main.pdf")).c_str()));
}

void VoreenBiologyMainWindow::helpWebsite() {
    // URL is stored in data portion of the sender
    if (QAction* act = dynamic_cast<QAction*>(sender()))
        QDesktopServices::openUrl(QUrl(act->data().toString()));
}

void VoreenBiologyMainWindow::helpAbout() {
    VoreenBiologyAboutBox about(this);
    about.exec();
}


void VoreenBiologyMainWindow::selectInitialWorkspace() {
    if (!VoreenApplication::app()->getShowStartupWizard()){
        workspaceHandler_->newWorkspace();
        return;
    }

    QStringList recentFileList = settings_.value("recentFileList").toStringList();
    QStringList standardWorkspaceList = getTemplateWorkspaces();
    VoreenBiologyStartupWizard w(recentFileList, standardWorkspaceList, this);
    
    QApplication::setOverrideCursor(Qt::ArrowCursor);
    bool shouldOpenWorkspace = w.exec();
    QApplication::restoreOverrideCursor();


    if (shouldOpenWorkspace){
        // load selected workspace
        QString selectedWorkspace = w.getSelectedWorkspace();
        workspaceHandler_->openWorkspace(selectedWorkspace);
    }else{
        // new workspace
        workspaceHandler_->newWorkspace();
    }

    workspaceLoadedSuccessfully();
}


} // namespace
