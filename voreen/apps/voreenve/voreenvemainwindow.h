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

#ifndef VRN_VOREENVEMAINWINDOW_H
#define VRN_VOREENVEMAINWINDOW_H

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

namespace voreen {

/**
 * The main window that contains all the other widgets, tool bar, menu bar, etc.
 */
class VoreenVEMainWindow : public VoreenQtMainWindow {
    Q_OBJECT
public:
    /**
     * Constructor
     *
     * @param workspace if passed by the command line parser (main.cpp) this workspace will be loaded
     * @param noInitialWorkspace
     * @param resetSettings if true, all voreen settings are resetted
     */
    VoreenVEMainWindow(const std::string& workspace = "", bool noInitialWorkspace = false, bool resetSettings = false);

    /**
     * Destructor. Calls deinitialize.
     */
    ~VoreenVEMainWindow();

    //-----------------------
    //  GUI mode changes
    //-----------------------
public slots:
    void guiModeChanged();
private:
    /**
     * Updates all widgets and menu entities / toolbars.
     */
    virtual void setGuiMode(GuiMode guiMode);

    //-----------------------
    //  GUI Menu/Toolbar
    //-----------------------
protected:
    /** Creates all basic menus. */
    virtual void createQMenus();
    /**
     * Let the user select a initial workspace.
     */
    virtual void selectInitialWorkspace();

    /** Creates all basic toolbars. */
    virtual void createQToolbars();

    // tool bars
    QToolBar* viewToolBar_;     ///< if needed allows switch between network and application
    QAction* modeApplicationAction_;
    QAction* modeNetworkAction_;

    //-----------------------
    //  menu help
    //-----------------------
protected slots:
    void helpFirstSteps();
    void helpNetworkEditor();
    void helpAnimation();
    void helpTutorialSlides();
    void helpWebsite();
    void helpAbout();

protected:
    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_VOREENVEMAINWINDOW_H

