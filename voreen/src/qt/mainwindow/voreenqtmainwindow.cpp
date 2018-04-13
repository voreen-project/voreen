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

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"
#include "voreen/qt/mainwindow/voreenqtsplashscreen.h"
#include "voreen/qt/mainwindow/voreenqtstartupwizard.h"

#include "voreen/qt/mainwindow/menuentities/voreenqtmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/animationmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/applicationmodeconfigmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/consolemenuentity.h"
#include "voreen/qt/mainwindow/menuentities/generalsettingsmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/inputmappingmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/openworkspacemenuentity.h"
#include "voreen/qt/mainwindow/menuentities/networkmodesettingsmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/newworkspacemenuentity.h"
#include "voreen/qt/mainwindow/menuentities/processorlistmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/propertylistmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/rebuildshadermenuentity.h"
#include "voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewermenuentity.h"
#include "voreen/qt/mainwindow/menuentities/saveworkspaceasmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/saveworkspacemenuentity.h"
#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/screenshotmenuentity.h"
#include "voreen/qt/mainwindow/menuentities/volumeviewermenuentity.h"
#include "voreen/qt/mainwindow/menuentities/workspacedescriptionmenuentity.h"

#include "tgt/gpucapabilities.h"
#include "tgt/filesystem.h"
#include "tgt/logmanager.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/network/processornetwork.h"

#include "voreen/core/network/workspace.h"
#include "voreen/core/network/networkevaluator.h"

// core module is always available
#include "modules/core/processors/output/canvasrenderer.h"
#include "modules/core/qt/processor/canvasrendererwidget.h"

#include "voreen/qt/widgets/consoleplugin.h"
#include "voreen/qt/widgets/inputmappingdialog.h"
#include "voreen/qt/widgets/volumeviewer.h"
#include "voreen/qt/widgets/voreensettingsdialog.h"
#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/animation/animationeditor.h"
#include "voreen/qt/widgets/processorlistwidget.h"
#include "voreen/qt/widgets/propertylistwidget.h"
#include "voreen/qt/networkeditor/networkeditor.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/qt/voreenapplicationqt.h"
#include "voreen/qt/voreenmoduleqt.h"

#include "voreen/core/version.h"
#include "voreen/core/utils/stringutils.h"

#include <QApplication>
#include <QDesktopServices>
#include <QVariant>
#include <QMenuBar>
#include <QMessageBox>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QFileInfo>
#include <QToolBar>
#include <QDesktopWidget>
#include <QScreen>

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

const std::string VoreenQtMainWindow::loggerCat_("voreenqt.VoreenQtMainWindow");

VoreenQtMainWindow::VoreenQtMainWindow(const std::string& title, const std::string& workspace, bool noInitialWorkspace, bool resetSettings)
    : QMainWindow()
    , currentGuiMode_(MODE_NONE)

    , applicationTitle_(title)

    , networkEvaluator_(0)
    , networkEditorWindow_(0)
    , networkEditorWidget_(0)
    , workspaceHandler_(0)

    , mainMenu_(0)
    , toolsToolBar_(0)

    , processorListWidget_(0)
    , propertyListWidget_(0)
    , volumeViewer_(0)
    , consolePlugin_(0)
    , animationEditor_(0)
    , noInitialWorkspace_(noInitialWorkspace)
    , initialized_(false)
{
    //general settings
    setMinimumSize(600, 400);
    setAcceptDrops(true);
    // show tooltips (for icons in toolbar) even if the window is inactive, as often the case
    // when a canvas window is active
    setAttribute(Qt::WA_AlwaysShowToolTips);
    setDockOptions(QMainWindow::AnimatedDocks); // disallow tabbed docks

    // initialialize the console early so it gets all the interesting messages
    consolePlugin_ = new ConsolePlugin(this, VoreenApplication::app()->getLogLevel());

    // if we have a stylesheet we want the fancy menu bar, please
    if (!qApp->styleSheet().isEmpty())
        setMenuBar(new FancyMenuBar());

    // clear session settings (window states, paths, ...), if specified by cmd line parameter
    if (resetSettings) {
        settings_.clear();
        LWARNING("Restored session settings");
    }
    else {
        loadSettings();
    }


    updateWindowTitle();
}

VoreenQtMainWindow::~VoreenQtMainWindow() {
    //delete visible actions
    foreach(QAction* act, toolBarVisibleMenu_->actions())
        delete act;
    //delete menu entities
    /*for(int i = 0; i < menuEntities_.size(); i++)
        delete menuEntities_[i];*/

    try {
        VoreenApplication::app()->deinitialize();
    }
    catch (VoreenException& e) {
        if (tgt::LogManager::isInited())
            LERROR("Failed to deinitialize VoreenApplication: " << e.what());
        std::cerr << "Failed to deinitialize VoreenApplication: " << e.what() << std::endl;
    }
}

void VoreenQtMainWindow::initialize(VoreenQtMainWindow::GuiMode initMode, VoreenQtSplashScreen* splash) {
    if (splash)
        splash->updateProgressMessage("Initializing OpenGL...", 0.50f);

    if (tgt::Singleton<tgt::LogManager>::isInited())
        LINFO("Log file: " << tgt::FileSystem::cleanupPath(LogMgr.getLogDir() + "/" + VoreenApplication::app()->getLogFile()));

    // initialize OpenGL
    try {
        VoreenApplication::app()->initializeGL();
    }
    catch(VoreenException& e) {
        if (tgt::LogManager::isInited())
            LFATALC("voreenve.MainWindow", "OpenGL initialization failed: " << e.what());
        else
            std::cerr << "OpenGL initialization failed: " << e.what();

        if (splash)
            splash->close();
        qApp->processEvents();
        QMessageBox::critical(this, tr("Initialization Error"), tr("OpenGL initialization failed. Quit."));

        exit(EXIT_FAILURE);
    }

    // some hardware/driver checks
    if (!GpuCaps.isOpenGlVersionSupported(tgt::GpuCapabilities::GlVersion::TGT_GL_VERSION_2_0)) {
        if (splash)
            splash->close();
        qApp->processEvents();
        std::ostringstream glVersion;
        glVersion << GpuCaps.getGlVersion();
        QMessageBox::critical(this, tr("Incompatible OpenGL Version"),
                              tr("Voreen requires OpenGL version 2.0 or higher, which does not seem to be "
                                 "supported on this system (reported version: %1). Therefore, the application "
                                 "will most likely not work properly.").arg(glVersion.str().c_str()));
        qApp->processEvents();
    }
    // Deactivated on intel until we have a reliable detection. (stefan)
    else if ( !GpuCaps.isShaderModelSupported(tgt::GpuCapabilities::SHADER_MODEL_3) && (GpuCaps.getVendor() != tgt::GpuCapabilities::GPU_VENDOR_INTEL) ) {
        if (splash)
            splash->close();
        qApp->processEvents();
        QMessageBox::critical(this, tr("Incompatible Shader Model"),
                              tr("Voreen requires Shader Model 3 or higher, which does not seem to be "
                                 "supported on this system. Therefore, the application will most likely not "
                                 "work properly."));
        qApp->processEvents();
    }
    else if (!GpuCaps.areFramebufferObjectsSupported()) {
        if (splash)
            splash->close();
        qApp->processEvents();
        QMessageBox::critical(this, tr("Framebuffer Objects Missing"),
                              tr("Voreen uses OpenGL framebuffer objects, which do not seem to be supported "
                                 "on this system. Therefore, the application will most likely not work properly."));
        qApp->processEvents();
    }

    if (splash)
        splash->updateProgressMessage("Creating visualization...", 0.60f);

    // create visualization object
    VoreenQtWorkspaceHandler::init();
    workspaceHandler_ = &WsHndlr;
    tgtAssert(workspaceHandler_, "no workspace handler");
    workspaceHandler_->registerMainWindow(this);
    // create network evaluator
    networkEvaluator_ = new NetworkEvaluator(true);
    VoreenApplication::app()->registerNetworkEvaluator(networkEvaluator_);
    workspaceHandler_->registerNetworkEvaluator(networkEvaluator_);

    if (splash)
        splash->updateProgressMessage("Creating user interface...", 0.80f);

    // mdi area
    mdiArea_ = new QMdiArea(this);
    mdiArea_->setOption(QMdiArea::DontMaximizeSubWindowOnActivation, true);
    setCentralWidget(mdiArea_);
    // put network editor in mdi area
    networkEditorWidget_ = new NetworkEditor(this, networkEvaluator_);
    networkEditorWidget_->setWindowTitle(tr("Processor Network"));
    networkEditorWindow_ = new QMdiSubWindow(this, Qt::FramelessWindowHint); // Deprecated? new VoreenMdiSubWindow(networkEditorWidget_, this, Qt::FramelessWindowHint);
    networkEditorWindow_->setWidget(networkEditorWidget_);
    networkEditorWindow_->setAttribute(Qt::WA_DeleteOnClose, false);
    networkEditorWindow_->setWindowState(networkEditorWindow_->windowState() | Qt::WindowFullScreen);
    networkEditorWindow_->hide(); // hide initially to prevent flicker
    mdiArea_->addSubWindow(networkEditorWindow_);
    workspaceHandler_->registerNetworkEditorWidget(networkEditorWidget_);
    qApp->processEvents();

    // create menu and toolbar
    createMainWindowMenusAndToolbars();
    qApp->processEvents();

    // signals indicating a change in network
    connect(workspaceHandler_, SIGNAL(workspaceHasBeenModifiedSignal()), this, SLOT(updateWindowTitle()));

    GuiMode lastMode;
    loadWindowSettings(&lastMode);
    loadToolbarSettings();

    if(initMode == MODE_NONE)
        setGuiMode(lastMode);
    else
        setGuiMode(initMode);

    // hide splash
    if (splash) {
        splash->updateProgressMessage("Initialization complete.", 1.f);
        qApp->processEvents();
        splash->close(); // hide to enable workspace interaction, e.g., regressiontest checkbox
    }

    selectInitialWorkspace();

    //
    // now the GUI is complete => load initial workspace
    //
    adjustMenuEntitiesAndToolBarsVisibility(currentGuiMode_, true);
    initialized_ = true;
}

void VoreenQtMainWindow::selectInitialWorkspace(){
    workspaceHandler_->newWorkspace();
}

void VoreenQtMainWindow::deinitialize() {
    if(!initialized_)
        return;

    // save widget settings first
    saveSettings();

    // deinitialize plugins (deleting not necessary, since done by the Qt parent)
    for (size_t i=0; i<customMenuEntities_.size(); i++) {
        VoreenQtMenuEntity* plugin = customMenuEntities_.at(i);
        if (plugin->isInitialized()) {
            try {
                LINFO("Deinitializing Voreen plugin '" << plugin->getName() << "'");
                plugin->deinitialize();
                plugin->initialized_ = false;
            }
            catch (tgt::Exception& e) {
                LERROR("Failed to deinitialize VoreenVE plugin '" << plugin->getName() << "': " << e.what());
            }
        }
    }

    // deinitialize plugins (deleting not necessary, since done by the Qt parent)
    for (size_t i=0; i<menuEntities_.size(); i++) {
        VoreenQtMenuEntity* entity = menuEntities_.at(i);
        if (entity->isInitialized()) {
            try {
                entity->deinitialize();
                delete entity;
            }
            catch (tgt::Exception& e) {
                LERROR("Failed to deinitialize VoreenVEMenuEntity '" << entity->getName() << "': " << e.what());
            }
        }
    }
    propertyListWidget_ = 0;
    processorListWidget_ = 0;
    volumeViewer_ = 0;
    animationEditor_ = 0;

    // free workspace, unregister network from widgets
    VoreenQtWorkspaceHandler::deinit();

    // delete network evaluator
    networkEvaluator_->setProcessorNetwork(0);
    delete networkEvaluator_;
    networkEvaluator_ = 0;

    // finalize OpenGL
    try {
        VoreenApplication::app()->deinitializeGL();
    }
    catch (VoreenException& e) {
        if (tgt::LogManager::isInited())
            LERROR("VoreenApplication::deinitializeGL failed: " << e.what());
        std::cerr << "VoreenApplication::deinitializeGL failed: " << e.what() << std::endl;
    }
}

NetworkEditor* VoreenQtMainWindow::getNetworkEditor() const {
    return networkEditorWidget_;
}

NetworkEvaluator* VoreenQtMainWindow::getNetworkEvaluator() const {
    return networkEvaluator_;
}

VoreenQtMainWindow::GuiMode VoreenQtMainWindow::getCurrentGuiMode() const {
    return currentGuiMode_;
}

//**********************************************************************************************
//      window title
//**********************************************************************************************
void VoreenQtMainWindow::updateWindowTitle() {
    QString title = tr(applicationTitle_.c_str());

    if (workspaceHandler_ &&workspaceHandler_->getWorkspace() && workspaceHandler_->getWorkspace()->isModified())
        title += " *";

    if (workspaceHandler_ && workspaceHandler_->getWorkspace() && !workspaceHandler_->getWorkspace()->getFilename().empty()) {
        QFileInfo f(QString::fromStdString(workspaceHandler_->getWorkspace()->getFilename())); // get filename without path
        title += " - " + f.fileName();
    }

    setWindowTitle(title);
}

std::string VoreenQtMainWindow::getApplicationTitle() const {
    return applicationTitle_;
}

void VoreenQtMainWindow::changeEvent(QEvent* event) {
    // Filter out window title changes which were done outside setWindowTitle (non-virtual) of
    // this class. This is used to prevent MDI windows from adding their title to the main
    // window title when maximized.
    if (event->type() == QEvent::WindowTitleChange) {
        if (windowTitle() != originalWindowTitle_)
            setWindowTitle(originalWindowTitle_);
    }
}

void VoreenQtMainWindow::setWindowTitle(const QString& title) {
    originalWindowTitle_ = title;
    QMainWindow::setWindowTitle(title);
}

//**********************************************************************************************
//      application settings
//**********************************************************************************************
void VoreenQtMainWindow::saveSettings() {
    // write version number of the config file format (might be useful someday)
    settings_.setValue("ConfigVersion", 1);

    //get the newest changes of window configurations
    switch(currentGuiMode_) {
    case MODE_APPLICATION:
        applicationModeState_ = saveState(WINDOW_STATE_VERSION);
        break;
    case MODE_NETWORK:
        networkModeState_ = saveState(WINDOW_STATE_VERSION);
        networkEditorWindowState_ = networkEditorWindow_->saveGeometry();
        break;
    default:
        tgtAssert(false,"save settings has been called without proper gui mode!");
        break;
    }

    //save main window settings
    settings_.beginGroup("MainWindow");
        settings_.setValue("size", size());
        settings_.setValue("pos", pos());
        settings_.setValue("maximized", (windowState() & Qt::WindowMaximized) != 0);
        settings_.setValue("lastUsedWorkspace", QString::fromStdString(workspaceHandler_->getWorkspace() ? workspaceHandler_->getWorkspace()->getFilename() : ""));
        settings_.setValue("applicationModeState", applicationModeState_);
        settings_.setValue("networkModeState", networkModeState_);
        settings_.setValue("networkEditorWindowState", networkEditorWindowState_);
        settings_.setValue("lastUsedGuiMode", currentGuiMode_);
    settings_.endGroup();
    //save pathes
    //save all tool windows
    settings_.beginGroup("Windows");
        for (int i=0; i < toolWindows_.size(); ++i) {
            if (!toolWindows_[i]->objectName().isEmpty()) {
                settings_.beginGroup(toolWindows_[i]->objectName());
                settings_.setValue("visible", toolWindows_[i]->isVisible());
                settings_.setValue("pos", toolWindows_[i]->pos());
                settings_.setValue("size", toolWindows_[i]->size());
                //special handling of property lod in network mode
                if(PropertyListWidget* widget = dynamic_cast<PropertyListWidget*>(toolWindows_[i]->child()))
                    settings_.setValue("propertyLevelOfDetail", widget->getLevelOfDetail());
                settings_.endGroup();
            }
        }
    settings_.endGroup();
    //store toolbar settings
    saveToolbarSettings();
}

void VoreenQtMainWindow::loadSettings() {
    // set defaults
    QSize windowSize = QSize(0, 0);
    QPoint windowPosition = QPoint(0, 0);
    bool windowMaximized = true;

    // restore settings
    settings_.beginGroup("MainWindow");
        windowSize = settings_.value("size", windowSize).toSize();
        windowPosition = settings_.value("pos", windowPosition).toPoint();
        windowMaximized = settings_.value("maximized", windowMaximized).toBool();
        applicationModeState_ = settings_.value("applicationModeState").toByteArray();
        networkModeState_ = settings_.value("networkModeState").toByteArray();
        networkEditorWindowState_ = settings_.value("networkEditorWindowState").toByteArray();
    settings_.endGroup();

    settings_.beginGroup("Startup");
    // load last startup values
    startupWorkspace_ = settings_.value("workspace", true).toBool();
    // set default values for the current startup
    settings_.setValue("workspace", false);
    settings_.endGroup();

    if (windowSize.isNull()) {
        windowMaximized = true;
    } else {
        resize(windowSize);
    }

    // ensure that the main window is restored on a visible screen
    // particular, when switching between different multi desktops modes
    QRect screenGeometry = QApplication::desktop()->screen()->geometry();

    // modify screen geometry to account maximized windows having negative position
    screenGeometry.setRect(screenGeometry.x() - 10, screenGeometry.y() - 10,
                           screenGeometry.width() + 20, screenGeometry.height() + 20);

    if (screenGeometry.contains(windowPosition) &&
        screenGeometry.contains(QPoint(windowPosition.x()+windowSize.width(),
                                       windowPosition.y()+windowSize.height())))
    {
        move(windowPosition);
    }

    if (windowMaximized)
        setWindowState(windowState() | Qt::WindowMaximized);
}

void VoreenQtMainWindow::loadWindowSettings(GuiMode* lastUsedMode) {
    // Restore visibility, position and size of tool windows from settings
    for (int i=0; i < toolWindows_.size(); ++i) {
        if (!toolWindows_[i]->objectName().isEmpty()) {
            toolWindows_[i]->setVisible(false);
        }
    }

    settings_.beginGroup("Windows");
    for (int i=0; i < toolWindows_.size(); ++i) {
        if (!toolWindows_[i]->objectName().isEmpty()) {
            settings_.beginGroup(toolWindows_[i]->objectName());
            if (settings_.contains("size"))
                toolWindows_[i]->resize(settings_.value("size").toSize());

            // Ignore position (0, 0) for invisible windows as otherwise all previously
            // invisible windows would be placed at (0, 0) after restarting the application.
            if (settings_.contains("pos") &&
                (settings_.value("pos").toPoint() != QPoint(0, 0) || settings_.value("visible").toBool()))
            {
                toolWindows_[i]->move(settings_.value("pos").toPoint());
            }

            if (settings_.contains("visible"))
                toolWindows_[i]->setVisible(settings_.value("visible").toBool());

            //special handling of property lod in network mode
            if(PropertyListWidget* widget = dynamic_cast<PropertyListWidget*>(toolWindows_[i]->child())) {
                widget->setLevelOfDetail(static_cast<Property::LevelOfDetail>(settings_.value("propertyLevelOfDetail", Property::LOD_DEFAULT).toInt()));
            }

            settings_.endGroup();
        }
    }
    settings_.endGroup();

    settings_.beginGroup("MainWindow");
    int guiMode = settings_.value("lastUsedGuiMode",MODE_NETWORK).toInt();
    settings_.endGroup();

    *lastUsedMode = static_cast<GuiMode>(guiMode);
    //setGuiMode  is called in initialize afterwards
}

void VoreenQtMainWindow::workspaceLoadedSuccessfully() {
    settings_.beginGroup("Startup");
        settings_.setValue("workspace", true);
    settings_.endGroup();
}

void VoreenQtMainWindow::saveToolbarSettings() {
    settings_.beginGroup("Toolbar");
        for(size_t i = 0; i < menuEntities_.size(); i++) {
            if(QAction* action = menuEntities_[i]->getToolBarAction()) {
                settings_.setValue(action->text(),action->isVisible());
            }
        }
        for(size_t i = 0; i < customMenuEntities_.size(); i++) {
            if(QAction* action = customMenuEntities_[i]->getToolBarAction()) {
                settings_.setValue(action->text(),action->isVisible());
            }
        }
    settings_.endGroup();
}

void VoreenQtMainWindow::loadToolbarSettings() {
    settings_.beginGroup("Toolbar");
        for(size_t i = 0; i < menuEntities_.size(); i++) {
            if(QAction* action = menuEntities_[i]->getToolBarAction()) {
                action->setVisible(settings_.value(action->text(), true).toBool());
            }
        }
        for(size_t i = 0; i < customMenuEntities_.size(); i++) {
            if(QAction* action = customMenuEntities_[i]->getToolBarAction()) {
                action->setVisible(settings_.value(action->text(), true).toBool());
            }
        }
    settings_.endGroup();
}

void VoreenQtMainWindow::adjustCanvasWidgets(GuiMode guiMode) {
    // adjust canvas widgets
    if (workspaceHandler_->getWorkspace() && workspaceHandler_->getWorkspace()->getProcessorNetwork()) {
        const ProcessorNetwork* network = workspaceHandler_->getWorkspace()->getProcessorNetwork();
        const ApplicationModeConfiguration& config = workspaceHandler_->getWorkspace()->getApplicationModeConfig();
        const std::vector<CanvasRenderer*>& canvasProcessors = network->getProcessorsByType<CanvasRenderer>();
        //store old meta data
        writeCanvasMetaData();

        //Step1: restore all data
        for (size_t i=0; i < canvasProcessors.size(); ++i) {
            QProcessorWidget* qpw = dynamic_cast<QProcessorWidget*>(canvasProcessors[i]->getProcessorWidget());
            if (qpw) {
                switch(guiMode) {
                case MODE_APPLICATION: {
                    //restore canvas size using the meta data
                    if(currentGuiMode_ == MODE_NETWORK){ //otherwise nothing to do
                        IVec2MetaData* positionMeta = dynamic_cast<IVec2MetaData*>(canvasProcessors[i]->getMetaDataContainer().getMetaData("preferedApplicationModeCanvasPosition"));
                        if (!positionMeta)
                            LDEBUG("adjustCanvasWidgets(): No meta data object returned");
                        else
                            qpw->setPosition(positionMeta->getValue().x, positionMeta->getValue().y);
                        IVec2MetaData* sizeMeta = dynamic_cast<IVec2MetaData*>(canvasProcessors[i]->getMetaDataContainer().getMetaData("preferedApplicationModeCanvasSize"));
                        if (!sizeMeta)
                            LDEBUG("adjustCanvasWidgets(): No meta data object returned");
                        else
                            qpw->setSize(sizeMeta->getValue().x, sizeMeta->getValue().y);
                        BoolMetaData* fullscreenMeta = dynamic_cast<BoolMetaData*>(canvasProcessors[i]->getMetaDataContainer().getMetaData("preferedApplicationModeCanvasFS"));
                        if (!sizeMeta)
                            LDEBUG("adjustCanvasWidgets(): No meta data object returned");
                        else
                            if(fullscreenMeta->getValue())
                                qpw->setWindowState(windowState() | Qt::WindowFullScreen);
                    }

                                       } break;
                case MODE_NETWORK: {
                    //restore canvas. if started in app mode, old values are loaded
                    if(currentGuiMode_ == MODE_APPLICATION){
                        IVec2MetaData* positionMeta = dynamic_cast<IVec2MetaData*>(canvasProcessors[i]->getMetaDataContainer().getMetaData("preferedNetworkModeCanvasPosition"));
                        if (!positionMeta)
                            LDEBUG("adjustCanvasWidgets(): No meta data object returned");
                        else
                            qpw->setPosition(positionMeta->getValue().x, positionMeta->getValue().y);
                        IVec2MetaData* sizeMeta = dynamic_cast<IVec2MetaData*>(canvasProcessors[i]->getMetaDataContainer().getMetaData("preferedNetworkModeCanvasSize"));
                        if (!sizeMeta)
                            LDEBUG("adjustCanvasWidgets(): No meta data object returned");
                        else
                            qpw->setSize(sizeMeta->getValue().x, sizeMeta->getValue().y);
                        BoolMetaData* fullscreenMeta = dynamic_cast<BoolMetaData*>(canvasProcessors[i]->getMetaDataContainer().getMetaData("preferedNetworkModeCanvasFS"));
                        if (!sizeMeta)
                            LDEBUG("adjustCanvasWidgets(): No meta data object returned");
                        else
                            if(fullscreenMeta->getValue())
                                qpw->setWindowState(windowState() | Qt::WindowFullScreen);
                    }
                    qpw->show();
                                   } break;
                default:
                    tgtAssert(false,"Unknown gui mode");
                    break;
                }
            } //if qpw
        } //for

        //Step2: handle main canvas (after restore to handle canvas-size.links correctly)
        for (size_t i=0; i < canvasProcessors.size(); ++i) {
            QProcessorWidget* qpw = dynamic_cast<QProcessorWidget*>(canvasProcessors[i]->getProcessorWidget());
            if (qpw) {
                switch(guiMode) {
                case MODE_APPLICATION:
                    //make selected or first canvas to main canvas (maximize)
                    if((config.getMainCanvas().empty() && i == 0) || !config.getMainCanvas().compare(canvasProcessors[i]->getClassName() + "::" + canvasProcessors[i]->getGuiName())) {
                        qpw->showMaximized();
                        QMdiSubWindow* subw = mdiArea_->addSubWindow(qpw, Qt::WindowTitleHint); // Qt::WindowTitleHint removes close button and menu icon
                        subw->showMaximized();
                        //disable bool property
                        Property* prop = canvasProcessors[i]->getProperty("showFullScreen");
                        if(prop) prop->setReadOnlyFlag(true);
                        return; //we are done here
                    }
                    break;
                case MODE_NETWORK:
                    //retore selected or first "main" canvas
                    if((config.getMainCanvas().empty() && i == 0) || !config.getMainCanvas().compare(canvasProcessors[i]->getClassName() + "::" + canvasProcessors[i]->getGuiName())) {
                        QObject* par = qpw->parent();
                        qpw->setParent(this);
                        if (dynamic_cast<QMdiSubWindow*>(par)) {
                            mdiArea_->removeSubWindow(qpw);
                            delete par;
                        }
                        qpw->setParent(this);
                        qpw->setWindowFlags(Qt::Tool);
                        qpw->showNormal();
                        static_cast<QWidget*>(qpw)->setVisible(true); //otherwise, the resizing does not work
                        //enable bool property
                        Property* prop = canvasProcessors[i]->getProperty("showFullScreen");
                        if(prop) prop->setReadOnlyFlag(false);
                        return; //we are done here
                    }
                    break;
                default:
                    tgtAssert(false,"Unknown gui mode");
                    break;
                }
            } //if qpw
        } //for

    } //if workspace
}

void VoreenQtMainWindow::writeCanvasMetaData() {
    if (workspaceHandler_->getWorkspace() && workspaceHandler_->getWorkspace()->getProcessorNetwork()) {
        const ProcessorNetwork* network = workspaceHandler_->getWorkspace()->getProcessorNetwork();
        const std::vector<CanvasRenderer*>& canvasProcessors = network->getProcessorsByType<CanvasRenderer>();

        for (size_t i=0; i < canvasProcessors.size(); ++i) {
            ProcessorWidget* pw = canvasProcessors[i]->getProcessorWidget();
            if (pw) {
                QProcessorWidget* qpw = dynamic_cast<QProcessorWidget*>(pw);
                if (qpw) {
                    switch(currentGuiMode_) {
                    case MODE_NETWORK: {
                        IVec2MetaData* sizeMeta = new IVec2MetaData(qpw->getSize());
                        canvasProcessors[i]->getMetaDataContainer().addMetaData("preferedNetworkModeCanvasSize",sizeMeta);
                        IVec2MetaData* positionMeta = new IVec2MetaData(qpw->getPosition());
                        canvasProcessors[i]->getMetaDataContainer().addMetaData("preferedNetworkModeCanvasPosition",positionMeta);
                        BoolMetaData* fullscreenMeta = new BoolMetaData(qpw->isFullScreen());
                        canvasProcessors[i]->getMetaDataContainer().addMetaData("preferedNetworkModeCanvasFS",fullscreenMeta);
                        } break;
                    case MODE_APPLICATION: {
                        IVec2MetaData* sizeMeta = new IVec2MetaData(qpw->getSize());
                        canvasProcessors[i]->getMetaDataContainer().addMetaData("preferedApplicationModeCanvasSize",sizeMeta);
                        IVec2MetaData* positionMeta = new IVec2MetaData(qpw->getPosition());
                        canvasProcessors[i]->getMetaDataContainer().addMetaData("preferedApplicationModeCanvasPosition",positionMeta);
                        BoolMetaData* fullscreenMeta = new BoolMetaData(qpw->isFullScreen());
                        canvasProcessors[i]->getMetaDataContainer().addMetaData("preferedApplicationModeCanvasFS",fullscreenMeta);
                        } break;
                    default:
                        tgtAssert(false, "Unknown gui mode");
                    }
                }
            }
        }

    }
}

//**********************************************************************************************
//      create GUI elements (Menu and Toolrbar)
//**********************************************************************************************
std::vector<VoreenQtMenuEntity*> VoreenQtMainWindow::getAllMenuEntities() {
    std::vector<VoreenQtMenuEntity*> tmp;
    tmp.insert(tmp.end(),menuEntities_.begin(),menuEntities_.end());
    tmp.insert(tmp.end(),customMenuEntities_.begin(),customMenuEntities_.end());
    return tmp;
}

void VoreenQtMainWindow::createMainWindowMenusAndToolbars() {
    //inizializes the qt container for the entities
    createQMenus();
    createQToolbars();
    //initializes the basic entities and user defined plug-ins
    createBasicMenuEntities();
    createCustomMenuEntities();
    initializeMenuEntities();
    LGL_ERROR;
}

void VoreenQtMainWindow::createQMenus() {
    //create menu
    mainMenu_ = menuBar();
    mainMenu_->setVisible(true);
    fileMenu_ = mainMenu_->addMenu(tr("&File"));
    viewMenu_ = mainMenu_->addMenu(tr("&View"));
    toolsMenu_ = mainMenu_->addMenu(tr("&Tools"));
    settingsMenu_ = mainMenu_->addMenu(tr("&Settings"));
    helpMenu_ = mainMenu_->addMenu(tr("&Help"));

    //
    // File menu
    //
    //action used as reference to insert file actions before it
    quitAction_ = new QAction(QIcon(":/qt/icons/exit.png"), tr("&Quit"), this);
    quitAction_->setShortcut(tr("Ctrl+Q"));
    quitAction_->setStatusTip(tr("Exit the application"));
    quitAction_->setToolTip(tr("Exit the application"));
    connect(quitAction_, SIGNAL(triggered()), this, SLOT(close()));
    fileMenu_->addAction(quitAction_);

    fileMenu_->addSeparator();

    // Recent files
    for (int i = 0; i < MAX_RECENT_FILES; i++) {
        recentFileActs_.append(new QAction(this));
        connect(recentFileActs_[i], SIGNAL(triggered()), this, SLOT(openRecentFile()));
        fileMenu_->addAction(recentFileActs_[i]);
    }
    updateRecentFiles();

    QAction* toolBarConfigAction = new QAction(tr("Configure Toolbar"), this);
    toolBarVisibleMenu_ = new QMenu();
    toolBarConfigAction->setMenu(toolBarVisibleMenu_);
    viewMenu_->addAction(toolBarConfigAction);
}

void VoreenQtMainWindow::createQToolbars() {
#ifdef __APPLE__
    const QSize iconSize = QSize(23,23);
#endif
    // file toolbar
    fileToolBar_ = addToolBar(tr("File"));
#ifdef __APPLE__
    fileToolBar_->setIconSize(iconSize);
#endif
    fileToolBar_->setObjectName("file-toolbar");

    // tools toolbar
    toolsToolBar_ = addToolBar(tr("Tools"));
#ifdef __APPLE__
    toolsToolBar_->setIconSize(iconSize);
#endif
    toolsToolBar_->setObjectName("tools-toolbar");
    QLabel* toolLabel = new QLabel(tr("   Tools "));
    toolLabel->setObjectName("toolBarLabel");
    toolsToolBar_->addWidget(toolLabel);
#ifdef __APPLE__ // we are on a mac system
    // HACK (Workaround) for Qt Mac Bug, makes MainWindow reappear
    // see for details:
    // http://bugreports.qt.nokia.com/browse/QTBUG-5069?page=com.atlassian.jira.plugin.system.issuetabpanels%3Aall-tabpanel
    show();
#endif
}

void VoreenQtMainWindow::createBasicMenuEntities() {
    //file
    menuEntities_.push_back(new NewWorkspaceMenuEntity());
    menuEntities_.push_back(new OpenWorkspaceMenuEntity());
    menuEntities_.push_back(new SaveWorkspaceMenuEntity());
    menuEntities_.push_back(new SaveWorkspaceAsMenuEntity());
    //basic sidebar elements
    menuEntities_.push_back(new ProcessorListMenuEntity());
    menuEntities_.push_back(new PropertyListMenuEntity());
    menuEntities_.push_back(new ConsoleMenuEntity(consolePlugin_)); //has been created in initialize
    menuEntities_.push_back(new VolumeViewerMenuEntity());
    menuEntities_.push_back(new WorkspaceDescriptionMenuEntity());
    //basic tools
    menuEntities_.push_back(new ScreenshotMenuEntity());
    menuEntities_.push_back(new AnimationMenuEntity());
    menuEntities_.push_back(new InputMappingMenuEntity());
    menuEntities_.push_back(new RenderTargetViewerMenuEntity());
    menuEntities_.push_back(new RebuildShaderMenuEntity());
    menuEntities_.push_back(new ApplicationModeConfigMenuEntity());
    //settings
    menuEntities_.push_back(new GeneralSettingsMenuEntity());
    menuEntities_.push_back(new NetworkModeSettingsMenuEntity());
}

void VoreenQtMainWindow::createCustomMenuEntities() {
    // retrieve VoreenVE plugins from application/modules
    VoreenApplicationQt* appQt = VoreenApplicationQt::qtApp();
    if (!appQt) {
        LERROR("VoreenApplicationQt not instantiated");
    }
    else {
        std::vector<std::string> cmeNames;
        for (size_t i=0; i<appQt->getQtModules().size(); i++) {
            VoreenModuleQt* qtModule = appQt->getQtModules().at(i);
            const std::vector<VoreenQtMenuEntity*>& moduleCustomME = qtModule->getCustomMenuEntities();
            for (size_t j=0; j<moduleCustomME.size(); j++) {
                customMenuEntities_.push_back(moduleCustomME.at(j));
                cmeNames.push_back(moduleCustomME.at(j)->getName());
            }
        }
        if(!customMenuEntities_.empty()) {
            LINFO("VoreenQt custom menu entities: " << strJoin(cmeNames, ", "));
        }
    }
}


void VoreenQtMainWindow::initializeMenuEntities() {
    //get all entities
    std::vector<VoreenQtMenuEntity*> entitiesToAdd;
    entitiesToAdd.insert(entitiesToAdd.end(),menuEntities_.begin(),menuEntities_.end());
    entitiesToAdd.insert(entitiesToAdd.end(),customMenuEntities_.begin(),customMenuEntities_.end());

    //initialize and set members
    for (size_t i=0; i<entitiesToAdd.size(); i++) {
        VoreenQtMenuEntity* entity = entitiesToAdd.at(i);

        entity->setMainWindow(this);
        try {
            if(!entity->isInitialized()) {
                entity->initialize();
                LGL_ERROR;
            }
        }
        catch (tgt::Exception& e) {
            LERROR(e.what());
            continue;
        }
        //add toolwindow
        if(entity->getToolWindow())
            toolWindows_ << entity->getToolWindow();
        //TODO: remove
        if(ProcessorListMenuEntity* proc = dynamic_cast<ProcessorListMenuEntity*>(entity))
            processorListWidget_ = static_cast<ProcessorListWidget*>(proc->getToolWindow()->child());
        if(PropertyListMenuEntity* prop = dynamic_cast<PropertyListMenuEntity*>(entity))
            propertyListWidget_ = static_cast<PropertyListWidget*>(prop->getToolWindow()->child());
        if(VolumeViewerMenuEntity* volview = dynamic_cast<VolumeViewerMenuEntity*>(entity))
            volumeViewer_ = static_cast<VolumeViewer*>(volview->getToolWindow()->child());
        if(WorkspaceDescriptionMenuEntity* wsd = dynamic_cast<WorkspaceDescriptionMenuEntity*>(entity))
            wsdWidget_ = static_cast<QWidget*>(wsd->getToolWindow()->child());
        if(AnimationMenuEntity* anim = dynamic_cast<AnimationMenuEntity*>(entity))
            animationEditor_ = static_cast<AnimationEditor*>(anim->getToolWindow()->child());
    }
    //handle file entities
    for (std::vector<VoreenQtMenuEntity*>::iterator it = entitiesToAdd.begin(); it != entitiesToAdd.end(); /* empty */) {
        if((*it)->getMenuCategory() == VoreenQtMenuEntity::MC_FILE) {
            //add menu action
            if((*it)->getMenuAction())
                fileMenu_->insertAction(quitAction_,(*it)->getMenuAction()); //HACK: insert to front
            //add tool bar action
            if((*it)->getToolBarAction()) {
                //handle action
                QAction* action = new QAction(tr((*it)->getName().c_str()),viewMenu_);
                action->setCheckable(true);
                action->setChecked(true);
                connect(action,SIGNAL(toggled(bool)),(*it)->getToolBarAction(),SLOT(setVisible(bool)));
                (*it)->getToolBarAction()->setData(qVariantFromValue((void*)action));
                connect((*it)->getToolBarAction(),SIGNAL(changed()),this,SLOT(changeToolbarVisibilityItem()));
                //add to toolbar
                toolBarVisibleMenu_->addAction(action);
                fileToolBar_->addAction((*it)->getToolBarAction());
            }
            //remove entity
            it = entitiesToAdd.erase(it);
        } else
            it++;
    }
    fileMenu_->insertSeparator(quitAction_);
    toolBarVisibleMenu_->addSeparator();

    //handle settings entities
    for (std::vector<VoreenQtMenuEntity*>::iterator it = entitiesToAdd.begin(); it != entitiesToAdd.end(); /* empty */) {
        if((*it)->getMenuCategory() == VoreenQtMenuEntity::MC_SETTINGS) {
            //add menu action
            if((*it)->getMenuAction())
                settingsMenu_->addAction((*it)->getMenuAction());
            //add toolbar entry
             if((*it)->getToolBarAction()) {
                //handle action
                QAction* action = new QAction(tr((*it)->getName().c_str()),viewMenu_);
                action->setCheckable(true);
                action->setChecked(true);
                connect(action,SIGNAL(toggled(bool)),(*it)->getToolBarAction(),SLOT(setVisible(bool)));
                (*it)->getToolBarAction()->setData(qVariantFromValue((void*)action));
                connect((*it)->getToolBarAction(),SIGNAL(changed()),this,SLOT(changeToolbarVisibilityItem()));
                //add to toolbar
                toolBarVisibleMenu_->addAction(action);
                fileToolBar_->addAction((*it)->getToolBarAction());
            }
            //remove entity
            it = entitiesToAdd.erase(it);
        } else
            it++;
    }

    //handle basic sidebar entities
    for (std::vector<VoreenQtMenuEntity*>::iterator it = entitiesToAdd.begin(); it != entitiesToAdd.end(); /* empty */) {
        //TODO: add to front
        if((*it)->getMenuCategory() == VoreenQtMenuEntity::MC_BASIC_SIDEBAR) {
            //add menu action
            if((*it)->getMenuAction())
                toolsMenu_->addAction((*it)->getMenuAction());
            //add tool bar action
            if((*it)->getToolBarAction()) {
                //handle action
                QAction* action = new QAction(tr((*it)->getName().c_str()),viewMenu_);
                action->setCheckable(true);
                action->setChecked(true);
                connect(action,SIGNAL(toggled(bool)),(*it)->getToolBarAction(),SLOT(setVisible(bool)));
                (*it)->getToolBarAction()->setData(qVariantFromValue((void*)action));
                connect((*it)->getToolBarAction(),SIGNAL(changed()),this,SLOT(changeToolbarVisibilityItem()));
                //add to toolbar
                toolBarVisibleMenu_->addAction(action);
                toolsToolBar_->addAction((*it)->getToolBarAction());
            }
            //remove entity
            it = entitiesToAdd.erase(it);
        } else
            it++;
    }
    toolsMenu_->addSeparator();
    toolBarVisibleMenu_->addSeparator();

    //handle basic tools entities
    for (std::vector<VoreenQtMenuEntity*>::iterator it = entitiesToAdd.begin(); it != entitiesToAdd.end(); /* empty */) {
        //TODO: add to front
        if((*it)->getMenuCategory() == VoreenQtMenuEntity::MC_BASIC_TOOL) {
            //add menu action
            if((*it)->getMenuAction())
                toolsMenu_->addAction((*it)->getMenuAction());
            //add tool bar action
            if((*it)->getToolBarAction()) {
                //handle action
                QAction* action = new QAction(tr((*it)->getName().c_str()),viewMenu_);
                action->setCheckable(true);
                action->setChecked(true);
                connect(action,SIGNAL(toggled(bool)),(*it)->getToolBarAction(),SLOT(setVisible(bool)));
                (*it)->getToolBarAction()->setData(qVariantFromValue((void*)action));
                connect((*it)->getToolBarAction(),SIGNAL(changed()),this,SLOT(changeToolbarVisibilityItem()));
                //add to toolbar
                toolBarVisibleMenu_->addAction(action);
                toolsToolBar_->addAction((*it)->getToolBarAction());
            }
            //remove entity
            it = entitiesToAdd.erase(it);
        } else
            it++;
    }
    toolsMenu_->addSeparator();
    toolBarVisibleMenu_->addSeparator();

    //handle VE plugins entities
    for (std::vector<VoreenQtMenuEntity*>::iterator it = entitiesToAdd.begin(); it != entitiesToAdd.end(); /* empty */) {
        //TODO: add to front
        if((*it)->getMenuCategory() == VoreenQtMenuEntity::MC_CUSTOM_PLUGIN) {
            //add menu action
            if((*it)->getMenuAction())
                toolsMenu_->addAction((*it)->getMenuAction());
            //add tool bar action
            if((*it)->getToolBarAction()) {
                //handle action
                QAction* action = new QAction(tr((*it)->getName().c_str()),viewMenu_);
                action->setCheckable(true);
                action->setChecked(true);
                connect(action,SIGNAL(toggled(bool)),(*it)->getToolBarAction(),SLOT(setVisible(bool)));
                (*it)->getToolBarAction()->setData(qVariantFromValue((void*)action));
                connect((*it)->getToolBarAction(),SIGNAL(changed()),this,SLOT(changeToolbarVisibilityItem()));
                //add to toolbar
                toolBarVisibleMenu_->addAction(action);
                toolsToolBar_->addAction((*it)->getToolBarAction());
            }
            //remove entity
            it = entitiesToAdd.erase(it);
        } else
            it++;
    }
    tgtAssert(entitiesToAdd.empty(),"Not all entities have been added!");
}

void VoreenQtMainWindow::adjustMenuEntitiesAndToolBarsVisibility(GuiMode modeToBeSet, bool forced) {
    //return, if we are in the current mode
    if(currentGuiMode_ == modeToBeSet && !forced)
        return;

    //get all entities
    std::vector<VoreenQtMenuEntity*> entitiesToAdd;
    entitiesToAdd.insert(entitiesToAdd.end(),menuEntities_.begin(),menuEntities_.end());
    entitiesToAdd.insert(entitiesToAdd.end(),customMenuEntities_.begin(),customMenuEntities_.end());

    switch(modeToBeSet) {
    case MODE_NETWORK: {
        // configure toolbar and menus
        for (size_t i=0; i<entitiesToAdd.size(); i++) {
            VoreenQtMenuEntity* entity = entitiesToAdd.at(i);
            if(entity->getMenuAction())
                entity->getMenuAction()->setVisible(true);
            if(entity->getToolBarAction()) {
                QAction* visAct = reinterpret_cast<QAction*>(entity->getToolBarAction()->data().value<void*>());
                entity->getToolBarAction()->setVisible(visAct->isChecked());
            }
        }
        //restore visibility first, then activate menu
        settingsMenu_->menuAction()->setVisible(true);
        toolsMenu_->menuAction()->setVisible(true);
        toolBarVisibleMenu_->menuAction()->setVisible(true);

        bool filesVisible = false;
        foreach(QAction* act, fileToolBar_->actions()) {
            if(act->isVisible() && !act->isSeparator()) {
                filesVisible = true;
                break;
            }
        }
        fileToolBar_->setVisible(filesVisible);

        bool toolsVisible = false;
        foreach(QAction* act, toolsToolBar_->actions()) {
            if(act->isVisible() && !act->isSeparator()) {
                toolsVisible = true;
                break;
            }
        }
        toolsToolBar_->setVisible(toolsVisible);

                       } break;
    case MODE_APPLICATION: {
        const ApplicationModeConfiguration* config = 0;
        config = (workspaceHandler_->getWorkspace() ? &(workspaceHandler_->getWorkspace()->getApplicationModeConfig()) : 0);
        //general updates
        toolBarVisibleMenu_->menuAction()->setVisible(false); //deactivate before changing the entity visibility

        // configure toolbar and menus
        for (size_t i=0; i<entitiesToAdd.size(); i++) {
            VoreenQtMenuEntity* entity = entitiesToAdd.at(i);
            bool visible = entity->getDefaultVisibilityInApplicationMode();
            if(config && config->isMenuEntityVisible(entity->getName()).first) {
                visible = config->isMenuEntityVisible(entity->getName()).second;
            }
            if(entity->getMenuAction())
                entity->getMenuAction()->setVisible(visible);
            if(entity->getToolBarAction()) {
                entity->getToolBarAction()->blockSignals(true);
                entity->getToolBarAction()->setVisible(visible);
                entity->getToolBarAction()->blockSignals(false);
            }
        }

        //check for visibility
        bool filesVisible = false;
        foreach(QAction* act, fileToolBar_->actions()) {
            if(act->isVisible() && !act->isSeparator()) {
                filesVisible = true;
                break;
            }
        }
        fileToolBar_->setVisible(filesVisible);

        bool settingsVisible = false;
        foreach(QAction* act, settingsMenu_->actions()) {
            if(act->isVisible() && !act->isSeparator()) {
                settingsVisible = true;
                break;
            }
        }
        settingsMenu_->menuAction()->setVisible(settingsVisible);

        bool toolsVisible = false;
        foreach(QAction* act, toolsMenu_->actions()) {
            if(act->isVisible() && !act->isSeparator()) {
                toolsVisible = true;
                break;
            }
        }
        toolsMenu_->menuAction()->setVisible(toolsVisible);
        toolsToolBar_->setVisible(toolsVisible);

                           } break;
    default:
        tgtAssert(false, "Unknown gui mode!");
    }

    emit guiModeUpdated(modeToBeSet);
}

VoreenToolWindow* VoreenQtMainWindow::getToolWindow(QWidget* childWidget) const {
    if (!childWidget)
        return 0;
    foreach(VoreenToolWindow* toolWindow, toolWindows_) {
        if (toolWindow->child() == childWidget)
            return toolWindow;
    }
    return 0;
}

void VoreenQtMainWindow::changeToolbarVisibilityItem() {
    //do nothing, if menu is invisible. see adjustMenuEntitiesAnd...
    if(toolBarVisibleMenu_->menuAction()->isVisible()) {
        QObject* obj = sender();
        QAction* act = dynamic_cast<QAction*>(obj);
        QAction* visAct = reinterpret_cast<QAction*>(act->data().value<void*>());
        visAct->setChecked(act->isVisible());
    }
}

//**********************************************************************************************
//      close event
//**********************************************************************************************
void VoreenQtMainWindow::closeEvent(QCloseEvent *event) {
    // allow for graceful exit when mainwindow is closed before initialization is complete
    if(!initialized_) {
        if (tgt::LogManager::isInited())
            LFATAL("Received close event before or during initialization! Shutting down.");
        else
            std::cerr << "MainWindow received close event before or during initialization! Shutting down." << std::endl;

        qApp->processEvents();
        exit(EXIT_FAILURE);
    }

    tgtAssert(workspaceHandler_, "no workspace handler");
    if (workspaceHandler_->getWorkspace() && workspaceHandler_->getWorkspace()->isModified()) {
        if (workspaceHandler_->askSave()) {
            event->accept();
            emit closeMainWindow();
        } else {
            event->ignore();
            return;
        }
    }

    deinitialize();

#ifdef __APPLE__
    // HACK to prevent exit time crash on OS X
    exit(0);
#endif

}

//**********************************************************************************************
//  workspace/path handling
//**********************************************************************************************
void VoreenQtMainWindow::setWorkspace(Workspace* workspace) {
    updateWindowTitle();

    // adjust canvas widgets (created during workspace load) to application mode
    if (currentGuiMode_ == MODE_APPLICATION) {
        adjustCanvasWidgets(MODE_APPLICATION);
        setGuiMode(MODE_APPLICATION);
    }

    if (workspace && !workspace->getFilename().empty())
        addToRecentFiles(QString::fromStdString(workspace->getFilename()));
    updateWindowTitle();

    networkEditorWidget_->selectPreviouslySelectedProcessors();
}

void VoreenQtMainWindow::openRecentFile() {
    QAction* action = qobject_cast<QAction*>(sender());
    if (action) {
        QString file(action->data().toString());
        workspaceHandler_->openWorkspace(file);
    }
}

void VoreenQtMainWindow::addToRecentFiles(const QString& filename) {
    QStringList files = settings_.value("recentFileList").toStringList();
    files.removeAll("");        // delete empty entries
    files.removeAll(filename);
    files.prepend(filename);
    while (files.size() > MAX_RECENT_FILES)
        files.removeLast();

    settings_.setValue("recentFileList", files);
    updateRecentFiles();
}

void VoreenQtMainWindow::updateRecentFiles() {
    QStringList files = settings_.value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), MAX_RECENT_FILES);
    for (int i = 0; i < numRecentFiles; ++i) {
        QString text = QString("&%1 %2").arg(i + 1).arg(QFileInfo(files[i]).fileName());
        recentFileActs_[i]->setText(text);
        recentFileActs_[i]->setData(files[i]);
        recentFileActs_[i]->setVisible(true);
    }
    for (int j = numRecentFiles; j < MAX_RECENT_FILES; ++j)
        recentFileActs_[j]->setVisible(false);
}


QStringList VoreenQtMainWindow::getTemplateWorkspaces(){
    std::string path = VoreenApplication::app()->getApplicationResourcePath()+"/workspaces/";
    std::vector<std::string> files = tgt::FileSystem::getRef().readDirectory(path, true);
    QStringList filesqt;
    for(size_t i = 0; i != files.size(); i++){
        std::string file = path+files[i];
        filesqt << QString(file.c_str());
    }
    return filesqt;
}

} // namespace
