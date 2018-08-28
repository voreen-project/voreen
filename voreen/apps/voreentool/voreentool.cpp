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

#ifdef VRN_MODULE_PYTHON
// Must come first!
#include "modules/python/pythonmodule.h"
#endif

#ifdef VRN_MODULE_FLOWREEN
#include "modules/flowreen/processors/streamline/streamlinesave.h"
#endif

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/networkconfigurator.h"
#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/utils/commandlineparser.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/utils/stringutils.h"

#include "voreen/qt/voreenapplicationqt.h"

// core module is always available
#include "modules/core/processors/output/canvasrenderer.h"
#include "modules/core/processors/output/imagesequencesave.h"
#include "modules/core/processors/output/geometrysave.h"
#include "modules/core/processors/output/textsave.h"
#include "modules/core/processors/output/volumesave.h"
#include "modules/core/processors/output/volumelistsave.h"
#include "modules/core/processors/output/octreesave.h"
//#include "modules/core/processors/output/octreetovvdstorage.h"

#include "tgt/init.h"
#include "tgt/logmanager.h"
#include "tgt/filesystem.h"
#include "tgt/glcontextmanager.h"

#include <QApplication>

#ifdef UNIX
#include <clocale>
#endif

#include <locale>
#include <codecvt>
#include <string>

using namespace voreen;

const std::string APP_BINARY  = "voreentool";
const std::string APP_NAME    = "VoreenTool";
const std::string APP_DESC    = "Command-line interface for the manipulation and execution of Voreen networks.";

const std::string loggerCat_  = "voreentool.main";

/**
 * Evaluates the passed network using the passed evaluator.
 *
 * @param network The network to execute.
 * @param networkEvaluator The evaluator to use for executing the network.
 * @param configurationOptions Property configuration to apply *before* executing the network (see command line option).
 * @param actionOptions ButtonProperties that are to be triggered *after* the network has been evaluated (see command line option).
 * @param triggerAllAction Triggers all ButtonProperties after the network has been evaluated.
 * @param triggerImageSaves Triggers an 'image save' event on all Canvases and ImageSequenceSaves after network evaluation.
 * @param triggerVolumeSaves Triggers an 'volume save' event on all VolumeSave and VolumeListSaves after network evaluation.
 * @param triggerGeometrySaves Triggers an 'geometry save' event on all GeometrySave processors after network evaluation.
 * @param pythonScriptFilename Runs the passed Python script after the network configuration has been applied.
 *
 * @throw VoreenException If network evaluation failed.
 */
void executeNetwork(ProcessorNetwork* network, NetworkEvaluator* networkEvaluator,
    const std::vector<std::string>& configurationOptions, const std::vector<std::string>& actionOptions,
    bool triggerAllActions, bool triggerImageSaves, bool triggerVolumeSaves, bool triggerGeometrySaves,
    const std::string& pythonScriptFilename, const std::vector<std::string>& pythonScriptArgs);

/**
 * Exits the application with state EXIT_FAILURE.
 * Before application termination the passed error message is logged
 * and the global resources are deleted.
 */
void exitFailure(const std::string& errorMsg);

// global resources that are freed by exitFailure()
NetworkEvaluator* networkEvaluator_ = 0;
QApplication* qtApp_ = 0;

//-------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

    // create Voreen application
    VoreenApplicationQt::setupPrerequisites();

    qtApp_ = new QApplication(argc, argv);
    VoreenApplicationQt vrnApp(APP_BINARY, APP_NAME, APP_DESC, argc, argv,
        VoreenApplication::ApplicationFeatures(VoreenApplication::APP_ALL /*&~ VoreenApplication::APP_WIDGETS*/));

    // prepare command line parser
    CommandLineParser* cmdParser = vrnApp.getCommandLineParser();
    tgtAssert(cmdParser, "no CommandLineParser");
    int cmdStyle = po::command_line_style::default_style ^ po::command_line_style::allow_guessing;
    cmdParser->setCommandLineStyle(static_cast<po::command_line_style::style_t>(cmdStyle));

    std::string workspacePath;
    cmdParser->addOption<std::string>("workspace,w", workspacePath, CommandLineParser::MainOption,
        "The workspace to evaluate.");

    bool runEventLoop = false;
    cmdParser->addFlagOption("run-event-loop", runEventLoop, CommandLineParser::MainOption,
        "Do not terminate application after workspace execution. May be used for interacting with the workspace, "
        "or for remote controlling.");

    std::vector<std::string> configurationOptions;
    cmdParser->addMultiOption<std::string>("configuration,c", configurationOptions, CommandLineParser::MainOption,
        "Network configuration as a list of property value assignments in the format: \n"
        "ProcessorName.PropertyName=value");

    std::vector<std::string> actionOptions;
    cmdParser->addMultiOption<std::string>("actions,a", actionOptions, CommandLineParser::MainOption,
        "List of ButtonProperties to be triggered after the network has been evaluated, in the format: \n"
        "ProcessorName.ButtonPropertyName");

    bool glMode = false;
    cmdParser->addFlagOption("opengl", glMode, CommandLineParser::MainOption,
        "Run the workspace in OpenGL mode: initialize OpenGL context and create canvases.");

    bool triggerAllActions = false;
    cmdParser->addFlagOption("trigger-all-actions", triggerAllActions, CommandLineParser::MainOption,
        "Trigger all ButtonProperties after the network has been evaluated.");

    bool triggerImageSaves = false;
    cmdParser->addFlagOption("trigger-imagesaves", triggerImageSaves, CommandLineParser::MainOption,
        "Trigger a \"save file\" event on all Canvases and ImageSequenceSave processors "
        "after the network has been evaluated (requires OpenGL mode).");

    bool triggerVolumeSaves = false;
    cmdParser->addFlagOption("trigger-volumesaves", triggerVolumeSaves, CommandLineParser::MainOption,
        "Trigger a \"save file\" event on all VolumeSave and VolumeListSave processors after the network has been evaluated.");

    bool triggerGeometrySaves = false;
    cmdParser->addFlagOption("trigger-geometrysaves", triggerGeometrySaves, CommandLineParser::MainOption,
        "Trigger a \"save file\" event on all GeometrySave and TextSave processors after the network has been evaluated.");

    std::string scriptFilename;
    std::vector<std::string> scriptArgs;
#ifdef VRN_MODULE_PYTHON
    cmdParser->addOption("script", scriptFilename, CommandLineParser::MainOption,
        "Run a Python script after the network configuration has been applied.");
    cmdParser->addMultiOption("scriptArgs", scriptArgs, CommandLineParser::MainOption,
        "Arguments passed to the python script (possibly multiple).");
#endif

    std::string workingDirectory;
    cmdParser->addOption("workdir",  workingDirectory, CommandLineParser::MainOption,
        "Absolute work directory of the workspace. If unspecified, the workspace file location is used.");

    std::string tempDataPath;
    cmdParser->addOption("tempdir", tempDataPath, CommandLineParser::MainOption,
        "Path to directory in which temporary data will be placed.");

    // init application
    try {
        vrnApp.initialize();
    }
    catch (VoreenException& e) {
        exitFailure("Failed to initialize application: " + std::string(e.what()));
    }

    // Change temp data path before doing anything else after initialization
    if(!tempDataPath.empty()) {
        vrnApp.setTempDataPath(tempDataPath);
    }

    // check parameters
    if (workspacePath.empty() && !runEventLoop) {
        vrnApp.deinitialize();
        std::cout << "\nUsage: either specify workspace path (-w) or activate event loop (--run-event-loop)\n\n";
        return 0;
    }

    // initialize OpenGL context and initializeGL VoreenApplication
    if (glMode) {
        try {
            vrnApp.initializeGL();
        }
        catch (VoreenException& e) {
            exitFailure("OpenGL initialization failed: " + std::string(e.what()));
        }
    }

    // create network evaluator
    networkEvaluator_ = new NetworkEvaluator(glMode);
    vrnApp.registerNetworkEvaluator(networkEvaluator_);

    // load and execute workspace, if specified
    Workspace* workspace = 0;
    if (!workspacePath.empty()) {
        workspace = new Workspace();
        workspacePath = tgt::FileSystem::absolutePath(workspacePath);
        try {
            LINFO("Loading workspace " << workspacePath);
            if (!workingDirectory.empty()) {
                workingDirectory = tgt::FileSystem::absolutePath(workingDirectory);
                LINFO("Workspace working path: " << workingDirectory);
            }
            workspace->load(workspacePath, workingDirectory);
        }
        catch (SerializationException& e) {
            delete workspace;
            exitFailure(e.what());
        }
        // log errors occurred during workspace deserialization
        for (size_t i=0; i<workspace->getErrors().size(); i++)
            LERROR(workspace->getErrors().at(i));

        // execute network
        try {
            executeNetwork(workspace->getProcessorNetwork(), networkEvaluator_, configurationOptions, actionOptions,
                triggerAllActions, triggerImageSaves && glMode, triggerVolumeSaves, triggerGeometrySaves,
                scriptFilename, scriptArgs);
        }
        catch (VoreenException& e) {
            delete workspace;
            exitFailure(e.what());
        }
    }

    // start Qt main loop, if in interactive mode
    if (runEventLoop) {
        LINFO("Running event loop...");
        qtApp_->exec();
    }

    // clean up
    LINFO("Deinitializing network ...");
    networkEvaluator_->deinitializeNetwork();
    delete networkEvaluator_;
    delete workspace;

    if (glMode) {
        LDEBUG("Deinitializing OpenGL ...");
        vrnApp.deinitializeGL();
    }

    vrnApp.deinitialize();
    delete qtApp_;
    qtApp_ = 0;
    return 0;
}

//-------------------------------------------------------------------------------------------------

static std::wstring utf8_to_utf16(const std::string& str) {
    // Note: This of course assumes that the input is utf8 encoded.
    std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
    return converter.from_bytes(str);
}

void executeNetwork(ProcessorNetwork* network, NetworkEvaluator* networkEvaluator,
        const std::vector<std::string>& configurationOptions, const std::vector<std::string>& actionOptions,
        bool triggerAllActions, bool triggerImageSaves, bool triggerVolumeSaves, bool triggerGeometrySaves,
        const std::string& pythonScriptFilename, const std::vector<std::string>& pythonScriptArgs) {

    tgtAssert(network, "no network passed (null pointer)");
    tgtAssert(networkEvaluator, "no network evaluator passed (null pointer)") ;

    // initialize network evaluator and assign network to it, which also initializes the processors
    LINFO("Initializing network ...");
    networkEvaluator->setProcessorNetwork(network);

// fixes problems with locale settings due to Qt (see http://doc.qt.io/qt-5/qcoreapplication.html#locale-settings)
#ifdef UNIX
    std::setlocale(LC_NUMERIC, "C");
#endif

    // apply network configuration
    NetworkConfigurator configurator(network);
    LINFO("Applying network configuration ...");

    for (size_t i=0; i<configurationOptions.size(); i++) {
        std::string optionStr = configurationOptions.at(i);
        LDEBUG("applying configuraton option: " << optionStr);
        try {
            configurator.setPropertyValue(optionStr);
        }
        catch (VoreenException& e) {
            throw VoreenException("in config option '" + optionStr + "': " + e.what());
        }
    }

    // run Python script
#ifdef VRN_MODULE_PYTHON
    wchar_t* empty_argv = nullptr;
    if (!pythonScriptFilename.empty()) {
        if (!PythonModule::getInstance())
            throw VoreenException("Failed to run Python script: PythonModule not instantiated");
        LINFO("Running Python script '" << pythonScriptFilename << "' ...");

        std::vector<std::wstring> wstringargs;
        wstringargs.push_back(utf8_to_utf16(pythonScriptFilename));
        for(const std::string& s : pythonScriptArgs) {
            wstringargs.push_back(utf8_to_utf16(s));
        }

        std::vector<wchar_t*> wargs;
        for(std::wstring& s : wstringargs) {
            wargs.push_back(&s[0]); //s.data() only in c++17
        }
        PythonModule::getInstance()->setArgv(wargs.size(), wargs.data());
        PythonModule::getInstance()->runScript(pythonScriptFilename, false);  //< throws VoreenException on failure
        LINFO("Python script finished.");

        // Be sure to clear future dangling references
        PythonModule::getInstance()->setArgv(0, &empty_argv);
    }
#endif

    // evaluate network
    LINFO("Evaluating network ...");
    try {
        networkEvaluator->process();
        //TODO: hack to test stereoscopy module. Has to be a loop until network is valid.
        networkEvaluator->process();
        networkEvaluator->process();
    }
    catch (std::exception& e) {
        throw VoreenException("exception during network evaluation: " + std::string(e.what()));
    }

    // trigger actions
    LINFO("Triggering post-evaluation actions ...");
    if (triggerAllActions) {
        std::vector<ButtonProperty*> buttonProps = network->getPropertiesByType<ButtonProperty>();
        for (size_t i=0; i<buttonProps.size(); i++) {
            LDEBUG("triggering ButtonProperty: " << buttonProps.at(i)->getFullyQualifiedID());
            buttonProps.at(i)->clicked();
        }
    }
    else {
        for (size_t i=0; i<actionOptions.size(); i++) {
            std::string optionStr = actionOptions.at(i);
            LDEBUG("triggering action: " << optionStr);
            try {
                configurator.triggerButtonProperty(optionStr);
            }
            catch (VoreenException& e) {
                throw VoreenException("in action '" + optionStr + "': " + e.what());
            }
        }
    }

    // take Canvas screenshots and write out imagesequences
    if (triggerImageSaves) {
        std::vector<CanvasRenderer*> canvasRenderers = network->getProcessorsByType<CanvasRenderer>();
        for (size_t i=0; i<canvasRenderers.size(); i++) {
            ButtonProperty* screenshotProp = dynamic_cast<ButtonProperty*>(canvasRenderers.at(i)->getProperty("saveScreenshot"));
            if (!screenshotProp) {
                LERROR("'saveScreenshot' property of CanvasRenderer '" << canvasRenderers.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering CanvasRenderer screenshot: " << canvasRenderers.at(i)->getID());
            screenshotProp->clicked();
        }

        std::vector<ImageSequenceSave*> imageSequenceSaves = network->getProcessorsByType<ImageSequenceSave>();
        for (size_t i=0; i<imageSequenceSaves.size(); i++) {
            ButtonProperty* screenshotProp = dynamic_cast<ButtonProperty*>(imageSequenceSaves.at(i)->getProperty("save"));
            if (!screenshotProp) {
                LERROR("'save' property of ImageSequenceSave '" << imageSequenceSaves.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of ImageSequenceSave: " << imageSequenceSaves.at(i)->getID());
            screenshotProp->clicked();
        }
    }

    // trigger VolumeSave and VolumeListSave processors (also OctreeSave and OctreeToVVDStorage)
    if (triggerVolumeSaves) {
        std::vector<VolumeSave*> volumeSaves = network->getProcessorsByType<VolumeSave>();
        for (size_t i=0; i<volumeSaves.size(); i++) {
            ButtonProperty* saveProp = dynamic_cast<ButtonProperty*>(volumeSaves.at(i)->getProperty("save"));
            if (!saveProp) {
                LERROR("'save' property of VolumeSave '" << volumeSaves.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of VolumeSave: " << volumeSaves.at(i)->getID());
            saveProp->clicked();
        }

        std::vector<VolumeListSave*> volumeListSaves = network->getProcessorsByType<VolumeListSave>();
        for (size_t i=0; i<volumeListSaves.size(); i++) {
            ButtonProperty* saveProp = dynamic_cast<ButtonProperty*>(volumeListSaves.at(i)->getProperty("save"));
            if (!saveProp) {
                LERROR("'save' property of VolumeListSave '" << volumeListSaves.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of VolumeListSave: " << volumeListSaves.at(i)->getID());
            saveProp->clicked();
        }

        std::vector<OctreeSave*> octreeSaves = network->getProcessorsByType<OctreeSave>();
        for (size_t i=0; i<octreeSaves.size(); i++) {
            ButtonProperty* saveProp = dynamic_cast<ButtonProperty*>(octreeSaves.at(i)->getProperty("save"));
            if (!saveProp) {
                LERROR("'save' property of OctreeSave '" << octreeSaves.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of OctreeSave: " << octreeSaves.at(i)->getID());
            saveProp->clicked();
        }

        /*
        std::vector<OctreeToVVDStorage*> octreeToVVDStorages = network->getProcessorsByType<OctreeToVVDStorage>();
        for (size_t i=0; i<octreeToVVDStorages.size(); i++) {
            ButtonProperty* saveProp = dynamic_cast<ButtonProperty*>(octreeToVVDStorages.at(i)->getProperty("save"));
            if (!saveProp) {
                LERROR("'save' property of OctreeToVVDStorage '" << octreeToVVDStorages.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of OctreeToVVDStorage: " << octreeToVVDStorages.at(i)->getID());
            saveProp->clicked();
        }
        */
    }

    // trigger GeometrySave and TextSave processors
    if (triggerGeometrySaves) {
        std::vector<GeometrySave*> geometrySaves = network->getProcessorsByType<GeometrySave>();
        for (size_t i=0; i<geometrySaves.size(); i++) {
            ButtonProperty* saveProp = dynamic_cast<ButtonProperty*>(geometrySaves.at(i)->getProperty("save"));
            if (!saveProp) {
                LERROR("'save' property of GeometrySave '" << geometrySaves.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of GeometrySave: " << geometrySaves.at(i)->getID());
            saveProp->clicked();
        }

        //TODO: give text save a own key word
        std::vector<TextSave*> textSaves = network->getProcessorsByType<TextSave>();
        for (size_t i=0; i<textSaves.size(); i++) {
            ButtonProperty* saveProp = dynamic_cast<ButtonProperty*>(textSaves.at(i)->getProperty("save"));
            if (!saveProp) {
                LERROR("'save' property of TextSave '" << textSaves.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of TextSave: " << textSaves.at(i)->getID());
            saveProp->clicked();
        }

#ifdef VRN_MODULE_FLOWREEN
        //TODO: give streamline save a own key word
        std::vector<StreamlineSave*> streamlineSaves = network->getProcessorsByType<StreamlineSave>();
        for (size_t i=0; i<streamlineSaves.size(); i++) {
            ButtonProperty* saveProp = dynamic_cast<ButtonProperty*>(streamlineSaves.at(i)->getProperty("saveButton"));
            if (!saveProp) {
                LERROR("'save' property of StreamlineSave '" << streamlineSaves.at(i)->getID() << "' not found");
                continue;
            }
            LDEBUG("triggering 'save' property of StreamlineSave: " << streamlineSaves.at(i)->getID());
            saveProp->clicked();
        }
#endif
    }
}

//-------------------------------------------------------------------------------------------------

void exitFailure(const std::string& errorMsg) {
    LFATAL(errorMsg);

    if (networkEvaluator_)
        networkEvaluator_->deinitializeNetwork();
    delete networkEvaluator_;

    if (VoreenApplication::app()->isInitializedGL()) {
        LDEBUG("Deinitializing OpenGL ...");
        VoreenApplication::app()->deinitializeGL();
    }

    VoreenApplication::app()->deinitialize();
    delete qtApp_;
    qtApp_ = 0;
    std::cerr << std::endl << "FAILURE: " << errorMsg << std::endl;
    exit(EXIT_FAILURE);
}
