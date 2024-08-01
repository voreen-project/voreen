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

#ifndef VRN_VOREENQTMENUENTITY_H
#define VRN_VOREENQTMENUENTITY_H

#include "voreen/qt/voreenqtapi.h"

#include "tgt/exception.h"

#include <QIcon>

class QAction;

namespace voreen {

class VoreenQtMainWindow;
class VoreenToolWindow;
class NetworkEvaluator;

/**
 * This is the base class for all menu entities and plugins used in the VoreenQtMainWindow.
 * Each entity can be associated with a (floating) tool window.
 */
class VRN_QT_API VoreenQtMenuEntity {

    friend class VoreenQtMainWindow;

public:
    /**
     * The category determines, to which submenu the entity is added
     */
    enum MenuCategory {
        MC_FILE,        ///< file menu, containing load/save etc.
        MC_SETTINGS,    ///< settings menu, containing options etc.
        MC_BASIC_SIDEBAR, ///< tool menu, processorlistwidget etc.
        MC_BASIC_TOOL,  ///< tool menu, render target vieweer etc. processorwidget etc.
        MC_CUSTOM_PLUGIN       ///< tool menu, plugins (define in modules)
    };

    VoreenQtMenuEntity();
    virtual ~VoreenQtMenuEntity();

    //------------------------------------------------------------------------------------------------
    // These functions should be overwritten by all subclasses
    //------------------------------------------------------------------------------------------------
    /** Icon used in menu and toolbar. Has to be ovrewritten. */
    virtual QIcon getIcon() const = 0;
    /** Name used in menu and toolbar. Has to be ovrewritten. */
    virtual std::string getName() const = 0;
    /** Shortcut. Default implementation is "no short cut". */
    virtual std::string getShortCut() const {return "";}
    /** Category. Default implementation is VEPlugin. */
    virtual MenuCategory getMenuCategory() const {return MC_CUSTOM_PLUGIN;}
    /** Defines where the window can be docked. Default is left and right. */
    virtual Qt::DockWidgetAreas getAllowedDockWidgetAreas() const {return static_cast<Qt::DockWidgetAreas>(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);}
    /** The initial window position. Default is no position. */
    virtual Qt::DockWidgetArea getInitialDockWidgetArea() const {return Qt::NoDockWidgetArea;}
    /** This flag determines, if this tool should be visible in application mode by default.*/
    virtual bool getDefaultVisibilityInApplicationMode() const {return false;}

protected:

    /**
     * OpenGL-dependent or time-consuming initializations should be performed here.
     * However, it is usually not necessary to override this function. It is called
     * by the VoreenMainWindow during application initialization.
     *
     * @throw tgt::Exception to indicate an initialization failure
     */
    virtual void initialize();

    /**
     * OpenGL-dependent or time-consuming deinitializations should be performed here.
     * However, it is usually not necessary to override this function. It is called
     * by the VoreenMainWindow during application initialization.
     *
     * @throw tgt::Exception to indicate an initialization failure
     */
    virtual void deinitialize();

    /**
     * Creates the action associated with this entity in the menu.
     * The default impementation opens/hides the toolwindow.
     * This function is called in initialize.
     */
    virtual QAction* createMenuAction();

    /**
     * Creates the action associated with this entity in the toolbar.
     * The default impementation copies the menu action. The toolbar action can be null.
     * This function is called in initialize.
     */
    virtual QAction* createToolbarAction();

    /**
     * Create the widget contained in the tool window.
     * The default implementation returns null.
     */
    virtual QWidget* createWidget() const {return 0;}

    /**
     * Creates the toolwindow and adds the widget(see createWidget).
     * If no widget has been created, no toolwindow is been created too.
     * This function is called in initialize. (should not be overwritten)
     */
    virtual VoreenToolWindow* createToolWindow(QWidget* widget);

    //----------------------------------
    //      getter and setter
    //----------------------------------
    /** Returns the menu action. */
    QAction* getMenuAction() const;
    /** Returns the toolbar action. (can be null) */
    QAction* getToolBarAction() const;
    /** Returns the tool window. (can be null) */
    VoreenToolWindow* getToolWindow() const;
    /** Returns whether the menu entity has already been initialized by the VoreenQtMainWindow. */
    bool isInitialized() const;
    /** Called by the VoreenQtMainWindow to assign itself. Should be called before initialized. */
    virtual void setMainWindow(VoreenQtMainWindow* mainWindow);

    //----------------------------------
    //      member
    //----------------------------------
    bool initialized_;                  ///< true, if initialize has been called successfully

    QAction*            menuAction_;    ///< action stored in the menu
    QAction*            toolbarAction_; ///< action stored in the toolbar (can be null)
    VoreenToolWindow*   toolWindow_;    ///< toolWindow, if createWidget has been implemented

    VoreenQtMainWindow* mainWindow_;    ///< reference to the main window
};

} // namespace

#endif // VRN_VOREENVEMENUENTITY_H
