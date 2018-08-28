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

#ifndef VRN_VOREENBIOLOGYMAINWINDOW_H
#define VRN_VOREENBIOLOGYMAINWINDOW_H

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

namespace voreen {

/**
 * The main window that contains all the other widgets, tool bar, menu bar, etc.
 */
class VoreenBiologyMainWindow : public VoreenQtMainWindow {
    Q_OBJECT
public:
    /**
     * Constructor
     * 
     * @param workspace if passed by the command line parser (main.cpp) this workspace will be loaded
     * @param noInitialWorkspace 
     * @param resetSettings if true, all voreen settings are resetted
     */
    VoreenBiologyMainWindow(const std::string& workspace = "", bool noInitialWorkspace = false, bool resetSettings = false);

    /**
     * Destructor. Calls deinitialize.
     */
    ~VoreenBiologyMainWindow();


    //-----------------------
    //  GUI mode changes     
    //-----------------------
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

    //-----------------------
    //  menu help
    //-----------------------
protected slots:
    void helpFirstSteps();
    void helpWebsite();
    void helpAbout();

    //-----------------------
    //  workspace handling
    //-----------------------
    virtual void selectInitialWorkspace();

    //-----------------------
    //  members
    //-----------------------
protected:
    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_VOREENBIOLOGYMAINWINDOW_H

