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

//VOREEN BIOLOGY
#include "voreenbiologyapplication.h"
#include "voreenbiologymainwindow.h"
#include "voreenbiologysplashscreen.h"

#include "tgt/filesystem.h"
#include "tgt/logmanager.h"

#include "voreen/core/version.h"
#include "voreen/core/utils/exception.h"
#include "voreen/core/utils/commandlineparser.h"

#include "tgt/timer.h"
#include "tgt/event/eventhandler.h"
#include "tgt/event/eventlistener.h"

#ifdef WIN32
    #include "DongleInterface.h"
#endif

#include <string>
#include <QMessageBox>

#ifdef UNIX
#include <clocale>
#endif

using namespace voreen;

const std::string loggerCat_("voreenbiology.main");

#ifdef WIN32
/** 
 * Class handling the dongle checks
 */
class DongleChecker : public tgt::EventListener {
    //------------
    //  Members
    //------------
    tgt::EventHandler dongleHandler_;           ///< handler for the time events
    tgt::Timer* dongleTimer_;                   ///< timer for the checks
    static const int FREQUENCY_MSEC = 120000;    ///< time between two minutes

public:
    /**
     * Constructor
     * Makes the initial checks and starts the timer
     */
    DongleChecker() : dongleHandler_(), dongleTimer_(0)
    {
        //do the first check
        doDongleCheck();

        //start timer if everthing is fine
        dongleHandler_.addListenerToBack(this);
        dongleTimer_ = VoreenApplication::app()->createTimer(&dongleHandler_);
        dongleTimer_->start(FREQUENCY_MSEC,1);
    }

    /** Destructor */
    virtual ~DongleChecker() {
        delete dongleTimer_;
    }

private:
    /** Time event for the regular checks */
    virtual void timerEvent(tgt::TimeEvent* te)  {
        doDongleCheck();
        dongleTimer_->start(FREQUENCY_MSEC,1);
    }

    /** Main function doing the final check. */
    void doDongleCheck() {
        CDongleInterface dongle;
        bool dongleIsFine = false;
        if (dongle.ConnectDongle()) {
            if (dongle.CheckFeature(7 /* Voreen */)) {
                dongleIsFine = true;
            } else {
                QMessageBox::critical(0, "Requested feature not found", "Voreen Biology requires the connection of a dongle from LaVision BioTec.\n"\
                                 "The connected dongle does not support Voreen Biology. Please make sure that the correct dongle is connected and "\
                                 "the license is not expired. Otherwise VOreen Biology will exit and unsaved changes will get lost.");
            }
        } else {
            QMessageBox::critical(0, "Dongle not found", "Voreen Biology requires the connection of a dongle from LaVision BioTec.\n"\
                                 "Please make sure that the required dongle is connected to a USB port of this PC. Otherwise Voreen Biology will exit "\
                                 "and unsaved changes will get lost.");
        }

        //second check
        if (dongle.ConnectDongle()) {
            if (dongle.CheckFeature(7 /* Voreen */)) {
                dongleIsFine = true;
            }
        }

        //exit on problem
        if(!dongleIsFine) {
            exit(EXIT_FAILURE);
        }
    }
};

#endif

int main(int argc, char** argv) {
    //disable argb visuals (Qt bug) fixes/works around invisible TF (etc) windows
#ifdef __unix__
    setenv ("XLIB_SKIP_ARGB_VISUALS", "1", 1);
#endif

    // create application
    VoreenBiologyApplication vapp(argc, argv);

#ifdef WIN32
    //start timer for dongle. All done in this class
    DongleChecker dongleChecker;
#endif

    vapp.setOverrideCursor(Qt::WaitCursor);

    // initialize application (also loads modules and initializes them)
    try {
        vapp.initialize();
    }
    catch (VoreenException& e) {
        if (tgt::LogManager::isInited())
            LFATALC(loggerCat_, "Failed to initialize VoreenApplication: " << e.what());
        std::cerr << "Failed to initialize VoreenApplication: " << e.what();
        exit(EXIT_FAILURE);
    }

    // splash screen
    VoreenBiologySplashScreen* splash = 0;
    bool showSplash = vapp.getShowSplashScreen();
    if (showSplash) {
        splash = new VoreenBiologySplashScreen();
        splash->updateProgressMessage("Creating application...",0.15);
        splash->show();
        qApp->processEvents();
    }

// fixes problems with locale settings due to Qt (see http://doc.qt.io/qt-5/qcoreapplication.html#locale-settings)
#ifdef UNIX
    std::setlocale(LC_NUMERIC, "C");
#endif

    // load and set style sheet
#if !defined(__APPLE__)
    QFile file(":/voreenbiology/widgetstyle/voreenbiology.qss");
    file.open(QFile::ReadOnly);
    QString styleSheet = QLatin1String(file.readAll());
    vapp.setStyleSheet(styleSheet);
#endif

#ifndef VRN_SHARED_LIBS
    // init Qt resources, if voreen_qt has been built as static lib
    Q_INIT_RESOURCE(vrn_qt);
    Q_INIT_RESOURCE(voreenbiology);
#endif

    // create and show mainwindow
    if (showSplash)
        splash->updateProgressMessage("Creating main window...",0.30);
    VoreenBiologyMainWindow mainWindow("", false, false);
    vapp.setMainWindow(&mainWindow);
    mainWindow.show();

    // initialize mainwindow (also calls VoreenApplication::initializeGL())
    mainWindow.initialize(VoreenQtMainWindow::MODE_APPLICATION, splash);

    vapp.restoreOverrideCursor();

    // hide splash
    if (showSplash){
        splash->showMessage("Initialization complete.",1);
        delete splash;
    }

    return vapp.exec();
}
