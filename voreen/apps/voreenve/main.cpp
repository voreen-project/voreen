/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

//VOREEN_VE

#ifdef VRN_MODULE_PYTHON
// Must come first!
#include "modules/python/pythonmodule.h"
#endif

#include "voreenveapplication.h"
#include "voreenvemainwindow.h"
#include "voreenvesplashscreen.h"

#include "tgt/filesystem.h"

#include "tgt/logmanager.h"

#include "voreen/core/version.h"
#include "voreen/core/utils/exception.h"
#include "voreen/core/utils/commandlineparser.h"

#include <QFile>

#include <string>

#ifdef UNIX
#include <clocale>
#endif

using namespace voreen;

const std::string loggerCat_("voreenve.main");

int main(int argc, char** argv) {
    //disable argb visuals (Qt bug) fixes/works around invisible TF (etc) windows
#ifdef __unix__
    setenv ("XLIB_SKIP_ARGB_VISUALS", "1", 1);
#endif

    // create application
    VoreenVEApplication::setupPrerequisites();
    VoreenVEApplication vapp(argc, argv);
    vapp.setOverrideCursor(Qt::WaitCursor);

    // add command line options
    tgtAssert(vapp.getCommandLineParser(), "no command line parser")
    std::string workspaceFilename;
    vapp.getCommandLineParser()->addOption("workspace,w", workspaceFilename, CommandLineParser::MainOption,
        "Loads the specified workspace");

    bool noInitialWorkspace = false;
    vapp.getCommandLineParser()->addFlagOption("no-workspace", noInitialWorkspace, CommandLineParser::MainOption,
        "Disables loading of an initial workspace on startup");

#ifdef VRN_MODULE_PYTHON
    std::string scriptFilename;
    vapp.getCommandLineParser()->addOption("script", scriptFilename, CommandLineParser::MainOption,
        "Runs a Python script right after the initial workspace has been loaded");
#endif

    bool resetSettings = false;
    vapp.getCommandLineParser()->addFlagOption("resetSettings", resetSettings, CommandLineParser::MainOption,
        "Restores window settings and default paths");

    bool useStylesheet;
    vapp.getCommandLineParser()->addOption("useStylesheet", useStylesheet, CommandLineParser::AdditionalOption,
        "Use VoreenVE style sheet", true, "true");

    // initialize application (also loads modules and initializes them)
    try {
        vapp.initialize();
    }
    catch (VoreenException& e) {
        if (tgt::LogManager::isInited())
            LFATALC("voreenve.main", "Failed to initialize VoreenApplication: " << e.what());
        std::cerr << "Failed to initialize VoreenApplication: " << e.what();
        exit(EXIT_FAILURE);
    }

    // splash screen
    VoreenVESplashScreen* splash = 0;
    bool showSplash = vapp.getShowSplashScreen();
    if (showSplash) {
        splash = new VoreenVESplashScreen();
        splash->updateProgressMessage("Creating application...",0.15f);
        splash->show();
        qApp->processEvents();
    }

// fixes problems with locale settings due to Qt (see http://doc.qt.io/qt-5/qcoreapplication.html#locale-settings)
#ifdef UNIX
    std::setlocale(LC_NUMERIC, "C");
#endif

    // load and set style sheet
#if !defined(__APPLE__)
    if (useStylesheet) {
        QFile file(":/voreenve/widgetstyle/voreen.qss");
        file.open(QFile::ReadOnly);
        QString styleSheet = QLatin1String(file.readAll());
        vapp.setStyleSheet(styleSheet);
    }
#endif

    // create and show mainwindow
    if (showSplash)
        splash->updateProgressMessage("Creating main window...",0.30f);
    VoreenVEMainWindow mainWindow(workspaceFilename, noInitialWorkspace, resetSettings);
    vapp.setMainWindow(&mainWindow);
    mainWindow.show();

    // initialize mainwindow (also calls VoreenApplication::initializeGL())
    mainWindow.initialize(VoreenQtMainWindow::MODE_NONE, splash);

    vapp.restoreOverrideCursor();

    // hide splash
    if (showSplash){
        splash->updateProgressMessage("Initialization complete.",1.f);
        splash->close();
        delete splash;
    }
    // run Python script
#ifdef VRN_MODULE_PYTHON
    if (!scriptFilename.empty()) {
        // first make sure that all Qt events have been processed
        qApp->processEvents();
        if (PythonModule::getInstance()) {
            try {
                LINFO("Running Python script '" << scriptFilename << "' ...");
                PythonModule::getInstance()->runScript(scriptFilename, false);
                LINFO("Python script finished.");
            }
            catch (VoreenException& e) {
                LERROR(e.what());
            }
        }
        else
            LERROR("Failed to run Python script: PythonModule not instantiated");
    }
#endif

    return vapp.exec();
}
