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

#include "voreenvemainwindow.h"
#include "voreenvesplashscreen.h"
#include "voreenveapplication.h"
#include "voreenvestartupwizard.h"

#include "voreen/qt/mainwindow/menuentities/voreenqtmenuentity.h"

#include "voreen/qt/networkeditor/networkeditor.h"
#include "voreen/qt/widgets/processorlistwidget.h"
#include "voreen/qt/widgets/propertylistwidget.h"
#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/widgets/volumeviewer.h"
#include "voreen/qt/widgets/consoleplugin.h"

#include "voreenveaboutbox.h"
#include "voreen/qt/helpbrowser.h"
#include "voreen/core/utils/commandlineparser.h"

#include <QDesktopServices>
#include <QMenuBar>
#include <QToolBar>
#include <QMdiSubWindow>

namespace voreen {

namespace {

const int MAX_RECENT_FILES = 8;

// Version number of restoring state of the main window.
// Increase when incompatible changes happen.
const int WINDOW_STATE_VERSION = 15;  // V4.0

} // namespace

////////// VoreenVEMainWindow ////////////////////////////////////////////////////////////

const std::string VoreenVEMainWindow::loggerCat_("voreenve.VoreenVEMainWindow");

VoreenVEMainWindow::VoreenVEMainWindow(const std::string& workspace, bool noInitialWorkspace, bool resetSettings)
    : VoreenQtMainWindow("VoreenVE", workspace, noInitialWorkspace, resetSettings)
    , viewToolBar_(0)
{
    //general settings
    setWindowIcon(QIcon(":/qt/icons/voreen-logo_64x64.png"));
}

VoreenVEMainWindow::~VoreenVEMainWindow() {
}

//**********************************************************************************************
//      switch bewteen network and application mode
//**********************************************************************************************
void VoreenVEMainWindow::guiModeChanged() {
    if (modeApplicationAction_->isChecked() && currentGuiMode_ != MODE_APPLICATION)
        setGuiMode(MODE_APPLICATION);
    else if (modeNetworkAction_->isChecked() && currentGuiMode_ != MODE_NETWORK)
        setGuiMode(MODE_NETWORK);
}

void VoreenVEMainWindow::setGuiMode(GuiMode guiMode) {
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

        modeApplicationAction_->setChecked(true);

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
                getToolWindow(volumeViewer_)->show();
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
        //if we are in application mode, store all window settings
        if (currentGuiMode_ == MODE_APPLICATION) {
            applicationModeState_ = saveState(WINDOW_STATE_VERSION);
        }
        // first update canvas widget
        adjustCanvasWidgets(MODE_NETWORK);

        if (!restoreState(networkModeState_, WINDOW_STATE_VERSION)) {
            if (getToolWindow(processorListWidget_))
                getToolWindow(processorListWidget_)->show();
            if (getToolWindow(propertyListWidget_))
                getToolWindow(propertyListWidget_)->show();
            if (getToolWindow(volumeViewer_))
                getToolWindow(volumeViewer_)->show();
            if (getToolWindow(consolePlugin_))
                getToolWindow(consolePlugin_)->show();
            if (getToolWindow(wsdWidget_))
                getToolWindow(wsdWidget_)->show();
        }
        //restore editor
        setUpdatesEnabled(false);
        if (networkEditorWindow_->restoreGeometry(networkEditorWindowState_))
            networkEditorWindow_->show();
        else
             networkEditorWindow_->showMaximized();
        networkEditorWidget_->setVisible(true); // only show now, so it immediately gets the correct size
        modeNetworkAction_->setChecked(true);

        //update property list widget
        if(propertyListWidget_) propertyListWidget_->setWidgetMode(PropertyListWidget::NETWORK);
        setUpdatesEnabled(true);
        break;
    default:
        tgtAssert(false, "unknown gui mode state!");
    }
    //update menu entities
    adjustMenuEntitiesAndToolBarsVisibility(guiMode);
    //set current gui mode
    currentGuiMode_ = guiMode;
}

//**********************************************************************************************
//      create GUI elements (Menu and Toolrbar)
//**********************************************************************************************
void VoreenVEMainWindow::createQMenus() {
    VoreenQtMainWindow::createQMenus();
    //
    // View menu
    //
    viewMenu_->addSeparator();
    modeNetworkAction_ = new QAction(QIcon(":/voreenve/icons/network-mode.png"),
                                     tr("&Network Mode"), this);
    modeNetworkAction_->setCheckable(true);
    modeNetworkAction_->setShortcut(tr("F3"));
    modeApplicationAction_ = new QAction(QIcon(":/voreenve/icons/application-mode.png"),
                                           tr("&Application Mode"), this);
    modeApplicationAction_->setCheckable(true);
    modeApplicationAction_->setShortcut(tr("F4"));
    QActionGroup* guiModeGroup = new QActionGroup(this);
    guiModeGroup->addAction(modeApplicationAction_);
    guiModeGroup->addAction(modeNetworkAction_);
    modeApplicationAction_->setChecked(true);
    connect(guiModeGroup, SIGNAL(triggered(QAction*)), this, SLOT(guiModeChanged()));
    viewMenu_->addAction(modeNetworkAction_);
    viewMenu_->addAction(modeApplicationAction_);

    //
    // Help menu
    //
    QAction* helpFirstStepsAct = new QAction(QIcon(":/qt/icons/help.png"), tr("&Getting Started..."), this);
    helpFirstStepsAct->setShortcut(tr("F1"));
    connect(helpFirstStepsAct, SIGNAL(triggered()), this, SLOT(helpFirstSteps()));
    helpMenu_->addAction(helpFirstStepsAct);
    QAction* helpNetworkEditorAct = new QAction(QIcon(":/qt/icons/help.png"), tr("&Network Editor..."), this);
    connect(helpNetworkEditorAct, SIGNAL(triggered()), this, SLOT(helpNetworkEditor()));
    helpMenu_->addAction(helpNetworkEditorAct);
    QAction* helpAnimationAct = new QAction(QIcon(":/qt/icons/video_export.png"), tr("&Animation Manual..."), this);
    connect(helpAnimationAct, SIGNAL(triggered()), this, SLOT(helpAnimation()));
    helpMenu_->addAction(helpAnimationAct);
    QAction* helpTutorialSlidesAct = new QAction(QIcon(":/qt/icons/pdf.png"), tr("&Tutorial Slides..."), this);
    connect(helpTutorialSlidesAct, SIGNAL(triggered()), this, SLOT(helpTutorialSlides()));
    helpMenu_->addAction(helpTutorialSlidesAct);
    helpMenu_->addSeparator();
    // Add some web links to the menu. Use the redirects (in the "go" directory) to be
    // independent of website reorganization.
    QAction* websiteAct = new QAction(tr("Voreen Website..."), this);
    websiteAct->setData(tr("http://voreen.uni-muenster.de"));
    connect(websiteAct, SIGNAL(triggered()), this, SLOT(helpWebsite()));
    helpMenu_->addAction(websiteAct);
    helpMenu_->addSeparator();
    QAction* aboutAction = new QAction(QIcon(":/qt/icons/about.png"), tr("&About VoreenVE..."), this);
    connect(aboutAction, SIGNAL(triggered()), this, SLOT(helpAbout()));
    helpMenu_->addAction(aboutAction);
}

void VoreenVEMainWindow::createQToolbars() {
    VoreenQtMainWindow::createQToolbars();

    // view toolbar
    viewToolBar_ = addToolBar(tr("View"));
#ifdef __APPLE__
    viewToolBar_->setIconSize(iconSize());
#endif
    viewToolBar_->setObjectName("view-toolbar");
    QLabel* viewLabel = new QLabel(tr("   View  "));
    viewLabel->setObjectName("toolBarLabel"); //if you change the name, the color changes to black. Y???????
    viewToolBar_->addWidget(viewLabel);
    viewToolBar_->addAction(modeNetworkAction_);
    viewToolBar_->addAction(modeApplicationAction_);
}

//**********************************************************************************************
//      help menu functions
//**********************************************************************************************
void VoreenVEMainWindow::helpFirstSteps() {
    QDesktopServices::openUrl(QUrl("http://www.uni-muenster.de/Voreen/documentation/user_interface.html"));
}

void VoreenVEMainWindow::helpNetworkEditor() {
    QDesktopServices::openUrl(QUrl("http://www.uni-muenster.de/Voreen/documentation/network_editor.html"));
}

void VoreenVEMainWindow::helpTutorialSlides() {
    QString path(VoreenApplication::app()->getBasePath("doc/voreen_tutorial_slides.pdf").c_str());
    QDesktopServices::openUrl(QUrl(QString::fromStdString("file:///") + path, QUrl::TolerantMode));
}

void VoreenVEMainWindow::helpAnimation() {
    QString path(VoreenApplication::app()->getBasePath("doc/animation/animation_manual.pdf").c_str());
    QDesktopServices::openUrl(QUrl(QString::fromStdString("file:///") + path, QUrl::TolerantMode));
    // old animation had a html version. Keep the code if we create a html version of the new documentation in the future
    /*QString path(VoreenApplication::app()->getBasePath("doc/animation/animation.html").c_str());
    HelpBrowser* help = new HelpBrowser(QUrl::fromLocalFile(path), tr("VoreenVE Animation Manual"));
    help->resize(1050, 700);
    help->show();
    connect(this, SIGNAL(closeMainWindow()), help, SLOT(close()));*/
}

void VoreenVEMainWindow::helpWebsite() {
    // URL is stored in data portion of the sender
    if (QAction* act = dynamic_cast<QAction*>(sender()))
        QDesktopServices::openUrl(QUrl(act->data().toString()));
}

void VoreenVEMainWindow::helpAbout() {
    VoreenVEAboutBox about(this);
    about.exec();
}

void VoreenVEMainWindow::selectInitialWorkspace(){
    CommandLineParser* cmdParser = VoreenApplication::app()->getCommandLineParser();


    if (!VoreenApplication::app()->getShowStartupWizard() || cmdParser->isOptionSet("no-workspace")){
        workspaceHandler_->newWorkspace();
        return;
    }

    if (cmdParser->isOptionSet("workspace")){
        std::string initialWorkspace;
        VoreenApplication::app()->getCommandLineParser()->getOptionValue<std::string>("workspace", initialWorkspace);
        workspaceHandler_->openWorkspace(QString::fromStdString(initialWorkspace));
        return;
    }


    QStringList recentFileList = settings_.value("recentFileList").toStringList();
    QStringList standardWorkspaceList = getTemplateWorkspaces();
    VoreenVeStartupWizard w(recentFileList, standardWorkspaceList, this);

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
}

} // namespace
