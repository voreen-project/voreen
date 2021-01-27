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

#ifdef WIN32
#include <windows.h>
#endif

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/voreenmodule.h"

#include "tgt/init.h"
#include "tgt/filesystem.h"
#include "tgt/shadermanager.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <set>

using namespace voreen;
using std::string;

int main(int argc, char** argv) {
    VoreenApplication app("descriptiontest", "descriptiontest", "", argc, argv, VoreenApplication::APP_ALL);
    app.initialize();
    LogMgr.getConsoleLog()->enableColors(false);

    std::stringstream tmpOut;
    std::string loggerCat_ = "descriptiontest";

    // specify modules whose processors should be included in the check
    std::set<std::string> MODULES_TO_CHECK;
    MODULES_TO_CHECK.insert("core");
    MODULES_TO_CHECK.insert("base");
    MODULES_TO_CHECK.insert("advancedlighting");
    MODULES_TO_CHECK.insert("connexe");
    //MODULES_TO_CHECK.insert("dcmtk");
    //MODULES_TO_CHECK.insert("deprecated");
    MODULES_TO_CHECK.insert("devil");
    MODULES_TO_CHECK.insert("dynamicglsl");
    //MODULES_TO_CHECK.insert("experimental");
    MODULES_TO_CHECK.insert("ffmpeg");
    //MODULES_TO_CHECK.insert("fontrendering");
    MODULES_TO_CHECK.insert("gdcm");
    MODULES_TO_CHECK.insert("itk");
    MODULES_TO_CHECK.insert("itk_generated");
    MODULES_TO_CHECK.insert("opencl");
    MODULES_TO_CHECK.insert("openmp");
    MODULES_TO_CHECK.insert("plotting");
    MODULES_TO_CHECK.insert("pvm");
    MODULES_TO_CHECK.insert("python");
    MODULES_TO_CHECK.insert("segy");
    MODULES_TO_CHECK.insert("staging");
    MODULES_TO_CHECK.insert("stereoscopy");
    MODULES_TO_CHECK.insert("tiff");
    MODULES_TO_CHECK.insert("volumelabeling");
    MODULES_TO_CHECK.insert("zip");

    size_t successful = 0;
    size_t skipped = 0;
    size_t numFailed = 0;
    //const ProcessorFactory::KnownClassesVector& knownClasses =  pf->getKnownClasses();
    const std::vector<const Processor*>& knownClasses = VoreenApplication::app()->getSerializableTypes<Processor>();

    for(size_t i=0; i<knownClasses.size(); ++i) {
        string classname = knownClasses[i]->getClassName();
        string processorModule = knownClasses[i]->getModuleName();
        std::string description = knownClasses[i]->getDescription();

        const VoreenModule* module = VoreenApplication::app()->getModule(processorModule);
        if (!module) {
            LERROR("Found no module " << processorModule << " for processor " << classname);
            return EXIT_FAILURE;
        }

        if (!MODULES_TO_CHECK.count(module->getDirName())) {
            tmpOut << "\t<testcase classname=\"description." << processorModule << "\" group=\"" << processorModule << "\" name=\"" << classname << "\">\n";
            tmpOut << "\t\t<system-out>\n";
            tmpOut << "skipped (module not tested: " << module->getDirName() << ")\n";
            tmpOut << "\t\t</system-out>\n";
            tmpOut << "\t</testcase>\n";
            ++skipped;
            continue;
        }

        if (knownClasses[i]->getCodeState() == Processor::CODE_STATE_BROKEN) {
            tmpOut << "\t<testcase classname=\"description." << processorModule << "\" group=\"" << processorModule << "\" name=\"" << classname << "\">\n";
            tmpOut << "\t\t<system-out>\n";
            tmpOut << "skipped (Codestate: Broken)\n";
            tmpOut << "\t\t</system-out>\n";
            tmpOut << "\t</testcase>\n";
            ++skipped;
            continue;
        }

        if (knownClasses[i]->getCodeState() == Processor::CODE_STATE_EXPERIMENTAL) {
            tmpOut << "\t<testcase classname=\"description." << processorModule << "\" time=\"" << 0.0f
                << "\" group=\"" << processorModule << "\" name=\"" << classname << "\">\n";
            tmpOut << "\t\t<system-out>\n";
            tmpOut << "skipped (Codestate: Experimental\n";
            tmpOut << "\t\t</system-out>\n";
            tmpOut << "\t</testcase>\n";
            ++skipped;
            continue;
        }

        tmpOut << "\t<testcase classname=\"description." << processorModule << "\" group=\""
               << processorModule << "\" name=\"" << classname << "\">\n";

        if(description != "") {
            ++successful;
        }
        else {
            ++numFailed;
            tmpOut << "\t\t<failure type=\"xyz\" message=\"No description!\"/>\n";
            LERROR("Class is testing/stable but has no description! " << classname);
        }

        tmpOut << "\t</testcase>\n";
    }

    std::ofstream xmlLog;
    xmlLog.open ("description.xml");
    xmlLog << "<testsuite name=\"description\" failures=\"" << numFailed << "\"  errors=\""
           << numFailed << "\" tests=\"" << knownClasses.size() << "\" skipped=\"" << skipped << "\">\n";
    xmlLog << tmpOut.str();
    xmlLog << "</testsuite>\n";
    xmlLog.close();

    app.deinitialize();

    if (numFailed == 0)
        return 0;
    else
        return EXIT_FAILURE;
}

