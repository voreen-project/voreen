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

#ifdef WIN32
#include <windows.h>
#endif

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/networkevaluator.h"

#include "tgt/tgt_gl.h"
#include <GL/freeglut.h>

#include "tgt/init.h"
#include "tgt/shadermanager.h"
#include "tgt/logmanager.h"
#include "tgt/glcontextmanager.h"
#include "tgt/glut/glutcanvas.h"
#include "tgt/glcontextmanager.h"

#include <iostream>
#include <fstream>
#include <ctime>

using namespace voreen;

class StringLog : public tgt::Log {
protected:
    std::stringstream str_;
public:
    StringLog(bool dateStamping = true, bool timeStamping= true, bool showCat = true, bool showLevel = true) {
        timeStamping_ = timeStamping;
        dateStamping_ = dateStamping;
        showCat_ = showCat;
        showLevel_ = showLevel;
    }
    virtual ~StringLog() {}
    bool isOpen() { return true; }

    std::string getContent() { return str_.str(); }
    void reset() { str_.clear(); str_.str(""); }
protected:
    virtual void logFiltered(const std::string &cat, tgt::LogLevel level, const std::string &msg, const std::string& /*extendedInfo*/) {
        if (dateStamping_)
            str_ << "[" << getDateString() << "] ";
        if (timeStamping_)
            str_ << "[" << getTimeString() << "] ";
        if (showCat_)
            str_ << cat << " ";
        if (showLevel_)
            str_ << "(" << getLevelString(level) << ") ";
        str_ << msg << '\n';
    }
};

//---------------------------------------------------------------------------

tgt::GLUTCanvas* canvas_ = 0;

//---------------------------------------------------------------------------

int main(int argc, char** argv) {
    VoreenApplication app("processorinittest", "ProcessorInitTest",
        "Creates, initializes, deinitializes and deletes all registered processors (with OpenGL).", argc, argv, VoreenApplication::APP_ALL);
    app.initialize();
    LogMgr.getConsoleLog()->enableColors(false);

    glutInit(&argc, argv);

#ifndef VRN_OPENGL_COMPATIBILITY_PROFILE
    glutInitContextVersion(VRN_OPENGL_CONTEXT_MAJOR_VERSION, VRN_OPENGL_CONTEXT_MINOR_VERSION);
    glutInitContextProfile(GLUT_CORE_PROFILE);
#endif

    // initialize canvas
    canvas_ = new tgt::GLUTCanvas("Voreen - The Volume Rendering Engine", tgt::ivec2(512, 512), tgt::GLCanvas::RGBADD);
    canvas_->init();

    // Initialize GL and set canvas as main context.
    app.initializeGL();

    // initialize the network evaluator
    NetworkEvaluator* networkEvaluator = new NetworkEvaluator();
    ProcessorNetwork* network = new ProcessorNetwork();
    networkEvaluator->setProcessorNetwork(network);

    std::stringstream tmpOut;

    clock_t start,finish;
    double time;
    double totalTime = 0.0;

    size_t successful = 0;
    size_t skipped = 0;
    size_t numFailed = 0;
    std::vector<std::string> failed;
    const std::vector<const Processor*>& knownClasses = VoreenApplication::app()->getSerializableTypes<Processor>();

    StringLog* log = new StringLog(false, false, true, true);
    log->addCat("", true);
    LogMgr.addLog(log);
    for(size_t i=0; i<knownClasses.size(); ++i) {
        const Processor* processor = knownClasses[i];
        std::string classname = processor->getClassName();
        std::string module = processor->getModuleName();

        if (processor->getCodeState() == Processor::CODE_STATE_BROKEN) {
            tmpOut << "\t<testcase g=\"Init." << module << "\" time=\"" << 0.0f
                   << "\" group=\"" << module << "\" name=\"" << classname << "\">\n";
            tmpOut << "\t\t<system-out>\n";
            tmpOut << "skipped (Codestate: Broken)";
            tmpOut << "\t\t</system-out>\n";
            LWARNINGC("processorinittest", "Skipping class " << classname);
            tmpOut << "\t</testcase>\n";
            ++skipped;
            continue;
        }

        if (processor->getCodeState() == Processor::CODE_STATE_EXPERIMENTAL) {
            tmpOut << "\t<testcase g=\"Init." << module << "\" time=\"" << 0.0f
                   << "\" group=\"" << module << "\" name=\"" << classname << "\">\n";
            tmpOut << "\t\t<system-out>\n";
            tmpOut << "skipped (Codestate: Experimental)";
            tmpOut << "\t\t</system-out>\n";
            LWARNINGC("processorinittest", "Skipping class " << classname);
            tmpOut << "\t</testcase>\n";
            ++skipped;
            continue;
        }

        LINFOC("processorinittest", "Creating class " << classname << " [" << (i+1) << "/" << knownClasses.size() << "]");
        Processor* p = dynamic_cast<Processor*>(VoreenApplication::app()->createSerializableType(classname));
        if (!p) {
            std::cout << "Failed to create processor: " << classname << std::endl;
            return EXIT_FAILURE;
        }
        network->addProcessor(p);
        LINFOC("processorinittest", "Initializing class " << classname);

        log->reset();
        start = clock();
        bool result = networkEvaluator->initializeNetwork();
        finish = clock();
        time = (double(finish)-double(start))/CLOCKS_PER_SEC;
        totalTime += time;
        tmpOut << "\t<testcase g=\"Init." << module << "\" time=\"" << time << "\" group=\""
               << module << "\" name=\"" << classname << "\">\n";
        tmpOut << "\t\t<system-out>\n";
        tmpOut << log->getContent();
        tmpOut << "\t\t</system-out>\n";
        if (result) {
            //tmpOut << "/>\n";
            ++successful;
        }
        else {
            ++numFailed;
            //tmpOut << ">\n";
            tmpOut << "\t\t<failure type=\"xyz\" message=\"Failed to initialize!\"/>\n";
            //tmpOut << "\t\t\ttest123\n";
            //tmpOut << "\t\t</failure>\n";
            failed.push_back(classname);
            LERRORC("processorinittest", "Failed to initialize class " << classname);
        }
        tmpOut << "\t</testcase>\n";

        if (p->isInitialized()) {
            LINFOC("processorinittest", "Deinitializing class " << classname);
            networkEvaluator->deinitializeNetwork();
        }

        LINFOC("processorinittest", "Deleting class " << classname);
        network->removeProcessor(p);
    }
    if (numFailed == 0)
        LINFOC("processorinittest", "All class initialized successfully. (" << skipped << " skipped)");
    else {
        LERRORC("processorinittest", successful << " of " << knownClasses.size() << " classes initialized successfully. ("
                << skipped << " skipped)");
        LERRORC("processorinittest", "Failed init:");
        for(size_t i=0; i<failed.size(); ++i)
            LERRORC("processorinittest", failed[i]);
    }

    std::ofstream xmlLog;
    xmlLog.open ("memcheck.xml");
    xmlLog << "<testsuite name=\"valgrind-memcheckinit\" failures=\"" << numFailed << "\" time=\"" << totalTime << "\" errors=\""
           << numFailed << "\" tests=\"" << knownClasses.size() << "\" skipped=\"" << skipped << "\">\n";
    xmlLog << tmpOut.str();
    xmlLog << "</testsuite>\n";
    xmlLog.close();
    networkEvaluator->setProcessorNetwork(0, true);
    delete networkEvaluator;
    delete network;

    app.deinitializeGL();
    delete canvas_;

    app.deinitialize();

    if (numFailed == 0)
        return 0;
    else
        return EXIT_FAILURE;
}
