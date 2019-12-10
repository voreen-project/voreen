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

#ifndef VRN_VOREENQTMAINWINDOW_H
#define VRN_VOREENQTMAINWINDOW_H

#include "voreen/qt/utils/voreenqtworkspacehandler.h"
#include "voreen/qt/voreenqtapi.h"

#include <QMainWindow>
#include <QSettings>
#include <vector>

class QMdiArea;
class QMdiSubWindow;

namespace voreen {

class VoreenQtSplashScreen;
class VoreenToolWindow;
class VoreenQtMenuEntity;

class ProcessorListWidget;
class VolumeViewer;
class PropertyListWidget;
class ConsolePlugin;
class InputMappingDialog;
class AnimationEditor;

//dialogs for the different settings
class VoreenSettingsDialog;
class ApplicationModeConfigDialog;
class NetworkModeConfigDialog;

/**
 * The main window that contains all the other widgets, tool bar, menu bar, etc.
 */
class VRN_QT_API VoreenQtMainWindow : public QMainWindow {
    Q_OBJECT
public:
    /**
     * The main window can be in three different states
     */
    enum GuiMode {
        MODE_NONE,          ///< initial state. Must be set to network or application
        MODE_NETWORK,       ///< network mode shows the processor items and is used by developer
        MODE_APPLICATION    ///< application mode is designed for end-user
    };
protected:
    GuiMode currentGuiMode_; ///< current mode of the application
public:
    /**
     * Constructor
     *
     * @param workspace if passed by the command line parser (main.cpp) this workspace will be loaded
     * @param noInitialWorkspace
     * @param resetSettings if true, all voreen settings are resetted
     */
    VoreenQtMainWindow(const std::string& title = "Voreen - Volume Rendering Engine", const std::string& workspace = "", bool noInitialWorkspace = false, bool resetSettings = false);

    /**
     * Destructor. Calls deinitialize.
     */
    ~VoreenQtMainWindow();

    /**
     * Initialize the entire main window
     *
     * - creates shared context
     * - checks openGL Version
     * - creates VoreenVisualization
     * - creates networkeditor and gui elements
     * - loads plugins
     *
     * @param initMode defines the startup mode. MODE_NONE meens last used mode.
     * @param splash if a splash screen is passed, the progress bar is beeing updated
     */
    virtual void initialize(GuiMode initMode, VoreenQtSplashScreen* splash = 0);

    /**
     * Deinitialize all plugins and openGL context
     */
    virtual void deinitialize();

    /** Returns the networkeditor */
    NetworkEditor* getNetworkEditor() const;
    /** Returns the NetworkEvaluator */
    NetworkEvaluator* getNetworkEvaluator() const;

    GuiMode getCurrentGuiMode() const;

    //-----------------------
    //  window title
    //-----------------------
public slots:
    /** sets the window title */
    void setWindowTitle(const QString& title);
    std::string getApplicationTitle() const;
protected slots:
    /** Updates the window after network modifications */
    void updateWindowTitle();
protected:
    /** Overwrites the window title. Needed? */
    void changeEvent(QEvent* event);
private:
    QString originalWindowTitle_;           ///< title used to restore title in changeEvent
    std::string applicationTitle_;          ///< base title of the application

    //-----------------------
    //  application settings
    //-----------------------
protected:
    /**
    * Saves all settings related to the main window and tool windows.
    * The function is called by deinitialize.
    */
    void saveSettings();
    /**
     * Stores the toolbar configuration.
     * Is called within saveSettings().
     */
    void saveToolbarSettings();
    /**
     * Loads all settings related to the main window.
     * The function is called by the constructor.
     *
     * @note currentGuiMode will be restored in loadWindowSettings()
     */
    void loadSettings();
    /**
     * Called in initialize.
     * Restores all windows and resets the currentGuiMode
     * @param lastUsedMode pointer to the last used mode, which has been loaded
     */
    void loadWindowSettings(GuiMode* lastUsedMode);
    /**
     * Called in initialize.
     * Restores the toolbar configuration.
     */
    void loadToolbarSettings();
    /**
     * Stores, if the workspace has been laoded successfully during the last run. Otherwise the standard.vws is loaded
     */
    void workspaceLoadedSuccessfully();
    /**
     * Let the user select a initial workspace.
     */
    virtual void selectInitialWorkspace();

public:
    /**
     * Function is used to write canvas meta data needed to restore canvas size.
     * @note will be called by save menu entities before the workspace is been saved.
     */
    void writeCanvasMetaData();

protected slots:
    /// Adjust the canvas widgets to the currently active gui mode.
    void adjustCanvasWidgets(GuiMode guiMode);

protected:
    QSettings settings_;                    ///< settings storing all at application related settings, like window position etc.
    QByteArray applicationModeState_;       ///< used to store and restore all windows in application mode
    QByteArray networkModeState_;           ///< used to store and restore all windows in network mode
    QByteArray networkEditorWindowState_;   ///< used to store and restore all "windows" of the network editor

    //-----------------------
    //  GUI mode changes
    //-----------------------
private:
    /**
     * Updates all widgets and menu entities / toolbars.
     * @see VoreenVEMainWindow
     */
     virtual void setGuiMode(GuiMode guiMode) = 0;

    //-----------------------
    //  GUI Menu/Toolbar
    //-----------------------
public:
    /** returns all menu entities and plugins */
    std::vector<VoreenQtMenuEntity*> getAllMenuEntities();
protected:
    /**
     * Creates and initializes all menus and toolbars.
     * Is called in VoreenVEMainWindow::initialie().
     * UpdateMenuEntities has to be called afterwards.
     */
    void createMainWindowMenusAndToolbars();

    /** Creates all basic menus. */
    virtual void createQMenus();

    /** Creates all basic toolbars. */
    virtual void createQToolbars();

    /**
     * Creates all basic entities.
     * @note have to be deleted in mainwindow destructor.
     */
    void createBasicMenuEntities();

    /**
     * Gets all plug-ins from all modules.
     * @note are deleted by each module
     */
    void createCustomMenuEntities();

    /** Inizializes the entities and assigns the main window. */
    void initializeMenuEntities();

    /**
     * Function is called in setGuiMode. It updates the menus and toolbars.
     * @param modeToBeSet updates the ME to this mode
     * @param forced always updates the visibility. Otherwise, currentMode_ equals modeToBeSet returns
     */
    void adjustMenuEntitiesAndToolBarsVisibility(GuiMode modeToBeSet, bool forced = false);

    /**
     * Returns the tool window that encloses the passed widget.
     * TODO: remove.
     */
    VoreenToolWindow* getToolWindow(QWidget* childWidget) const;

protected slots:
    /** Used to hide the toolbar item, if it was stored hidden in the settings. */
    void changeToolbarVisibilityItem();

protected:
    QList<VoreenToolWindow*> toolWindows_;                  ///< each tool is wrapped by a tool window
    std::vector<VoreenQtMenuEntity*> menuEntities_;         ///< all basic menu entities
    std::vector<VoreenQtMenuEntity*> customMenuEntities_;   ///< all custom menu entities (no delete is called)

    // menus
    QMenuBar* mainMenu_;        ///< main MenuBar, containing all menu points
    QMenu* fileMenu_;           ///< sub menu for file menu entities
    QMenu* viewMenu_;           ///< sub menu for view menu entities
    QMenu* toolsMenu_;          ///< sub menu for tool menu entities
    QMenu* settingsMenu_;       ///< sub menu for setting menu entities
    QMenu* toolBarVisibleMenu_; ///< menu used to configure toolbar visibility
    QMenu* helpMenu_;           ///< sub menu for help menu entities
    // tool bars
    QToolBar* fileToolBar_;     ///< toolbar containing the file enitities
    QToolBar* toolsToolBar_;    ///< customizable toolbar

    //-----------------------
    //  close event
    //-----------------------
signals:
    void closeMainWindow();
    void guiModeUpdated(VoreenQtMainWindow::GuiMode mode);

protected:
    void closeEvent(QCloseEvent* event);

    QAction* quitAction_;

    //-----------------------
    //  workspace handling
    //-----------------------
public slots:
    void setWorkspace(Workspace* workspace);
    void openRecentFile();
protected:
    void addToRecentFiles(const QString& filename);
    void updateRecentFiles();
    QStringList getTemplateWorkspaces();

    QList<QAction*> recentFileActs_;
    VoreenQtWorkspaceHandler* workspaceHandler_;

    //-----------------------
    //  members
    //-----------------------
protected:
    NetworkEvaluator* networkEvaluator_;
    // network editor
    QMdiArea* mdiArea_;
    QMdiSubWindow* networkEditorWindow_;
    NetworkEditor* networkEditorWidget_; ///< to allow access from rtv

    ProcessorListWidget* processorListWidget_;
    PropertyListWidget* propertyListWidget_;
    VolumeViewer* volumeViewer_;
    QWidget* wsdWidget_;
    ConsolePlugin* consolePlugin_;
    AnimationEditor* animationEditor_;

    bool noInitialWorkspace_;
    bool startupWorkspace_;
    bool initialized_; ///< to allow for graceful exit if program was closed before initialization was complete

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_VOREENQTMAINWINDOW_H

