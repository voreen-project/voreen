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

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/version.h"
#include "voreen/core/utils/commandlineparser.h"
#include "voreen/core/utils/commandqueue.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/utils/memoryinfo.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/propertywidget.h"
#include "voreen/core/properties/propertywidgetfactory.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/processors/cache.h"

// core module is always available
#include "modules/core/coremodule.h"
#include "modules/core/processors/input/octreecreator.h"

#include "tgt/init.h"
#include "tgt/filesystem.h"
#include "tgt/glcanvas.h"
#include "tgt/timer.h"
#include "tgt/gpucapabilities.h"
#include "tgt/shadermanager.h"
#include "tgt/glcontextmanager.h"

// volume memory management
#include "voreen/core/memorymanager/volumememorymanager.h"

// file system watch
#include "voreen/core/utils/voreenfilewatcher.h"

#include "gen_moduleregistration.h"

#include <string>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

#ifdef WIN32
    #include <shlobj.h>
    #include <windows.h>
#else
    #include <stdlib.h>
#endif

#ifdef __APPLE__
    #include "CoreFoundation/CFBundle.h"
#endif

using std::string;

namespace {

string findWithSubDir(const string& path, const string& subdir, int iterations = 0) {
    string p = path;

    // try in start directory
    if (tgt::FileSystem::dirExists(p + "/" + subdir))
        return p;

    // now try parent dirs
    for (int i = 0; i < iterations; i++) {
        p += "/..";
        if (tgt::FileSystem::dirExists(p + "/" + subdir))
            return p;
    }

    return "";
}

// Convert define into a string, compare
// http://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define VRN_XSTRINGIFY(s) VRN_STRINGIFY(s)
#define VRN_STRINGIFY(s) #s

string findBasePath(const string& path) {
    return findWithSubDir(path, "resource/voreencore", 7);
}

#ifdef __APPLE__
string findAppBundleResourcesPath() {
    CFBundleRef bundle;

    bundle = CFBundleGetMainBundle();
    if(!bundle)
        return "";

    CFURLRef resourceURL = CFBundleCopyResourcesDirectoryURL(bundle);
    char* path = new char[200];
    if(!CFURLGetFileSystemRepresentation(resourceURL, true, (UInt8*)path, 200))
        return "";
    CFRelease(resourceURL);

    string pathStr;
    if (path)
        pathStr = string(path);

    delete[] path;
    return pathStr;
}

#endif

} // namespace anonymous

namespace tgt {

/// istream operator for converting string to LogLevel (used by CommandLineParser)
std::istream& operator>>(std::istream& in, tgt::LogLevel& level) {
    std::string token;
    in >> token;
    std::string tokenLower = voreen::toLower(token);
    if (tokenLower == "debug")
        level = tgt::Debug;
    else if (tokenLower == "info")
        level = tgt::Info;
    else if (tokenLower == "warning")
        level = tgt::Warning;
    else if (tokenLower == "error")
        level = tgt::Error;
    else if (tokenLower == "fatal")
        level = tgt::Fatal;
    else
#ifndef VRN_OLD_BOOST
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value, token);
#else
    ;
#endif
    return in;
}

std::string loglevelToString(tgt::LogLevel level) {
    if (level == tgt::Debug)
        return "debug";
    else if (level == tgt::Info)
        return "info";
    else if (level == tgt::Warning)
        return "warning";
    else if (level == tgt::Error)
        return "error";
    else if (level == tgt::Fatal)
        return "fatal";
    else
        return "unknown";
}

} // namespace tgt

namespace voreen {

// define static members
VoreenApplication* VoreenApplication::app_ = 0;
const std::string VoreenApplication::loggerCat_ = "voreen.VoreenApplication";

static const std::string VOREEN_LOCK_NAME = "voreen_lock";

void VoreenApplication::printAsciiLogo() const {
#ifdef WIN32
        HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
        if (hStdOut != INVALID_HANDLE_VALUE)
          SetConsoleTextAttribute(hStdOut, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
#else
        std::cout << std::string("\033[1;33m"); // switch to yellow
#endif
        std::cout << std::string(
                                #include "resource/voreencore/fonts/logo.h"
                                ) << std::endl;
#ifdef WIN32
        if (hStdOut != INVALID_HANDLE_VALUE)
          SetConsoleTextAttribute(hStdOut, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#else
         std::cout << std::string("\033[00m"); // switch back to default
#endif
}

VoreenApplication::VoreenApplication(const std::string& binaryName, const std::string& guiName, const std::string& description,
                                     int argc, char** argv, ApplicationFeatures appType)
    : PropertyOwner(binaryName, guiName)
    , appFeatures_(appType)
    , binaryName_(binaryName)
    , cmdParser_(new CommandLineParser(binaryName, description, po::command_line_style::default_style))
    , cmdQueue_(new CommandQueue())
    , enableLogging_("enableLogging", "Enable Logging", true, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    , logLevel_("logLevel", "Log Level", Processor::INVALID_RESULT, false, Property::LOD_DEVELOPMENT)
    , enableHTMLLogging_("htmlLogging", "Enable HTML File Logging", true, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    , htmlLogFile_("htmlLogFile", "HTML Log File", "Select HTML Log File", "", ".html", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    , overrideGLSLVersion_("")
    , loadModules_(true)
#ifdef VRN_DEPLOYMENT
    , deploymentMode_(true)
#else
    , deploymentMode_(false)
#endif
    , networkEvaluators_()
    , schedulingTimer_(0)
    , eventHandler_()
    , useCaching_("useCaching", "Use Caching", true, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    , volumeCacheLimit_("cacheLimit", "Volume Cache Size (GB)", 10, 0, 1000 /*1 TB*/, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_APPLICATION)
    , octreeCacheLimit_("octreeCacheLimit", "Octree Cache Size (GB)", 100, 0, 1000 /*1 TB*/, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_APPLICATION)
    , cachePath_("cachePath", "Cache Directory", "Select Cache Directory...", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_APPLICATION, VoreenFileWatchListener::ALWAYS_OFF)
    , resetCachePath_("resetCachePath", "Reset Cache Path" , Processor::INVALID_RESULT, Property::LOD_APPLICATION)
    , deleteCache_("deleteCache", "Delete Cache", Processor::INVALID_RESULT, Property::LOD_APPLICATION)
    , determineCpuRamLimitAutomatically_("automaticRamLimit", "Determine Volume CPU RAM Limit automatically (recommended)", true)
    // set maximum to 128 GB minus some overhead (operating system, other Voreen stuff etc.)
    , cpuRamLimit_("cpuRamLimit", "Volume CPU RAM Limit (MB)", 4000, std::min(1000, static_cast<int>(MemoryInfo::getTotalPhysicalMemory() / (1024 * 1024))),
            MemoryInfo::getTotalPhysicalMemory() ? (static_cast<int>(MemoryInfo::getTotalPhysicalMemory() / (1024 * 1024))) : 120000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_APPLICATION)
    , determineGpuMemoryLimitAutomatically_("automaticGpuMemoryLimit", "Determine Volume GPU Texture Memory Limit automatically (recommended)", true)
    , gpuMemoryLimit_("gpuMemoryLimit", "Volume GPU Texture Memory Limit (MB)", 1000, 100, 4000,
            Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC, Property::LOD_APPLICATION)
    , availableGraphicsMemory_("availableGraphicsMemory", "Available Graphics Memory (MB)", -1, -1, 10000, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC, Property::LOD_DEVELOPMENT)
    , refreshAvailableGraphicsMemory_("refreshAvailableGraphicsMemory", "Refresh", Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    , testDataPath_("testDataPath", "Test Data Directory", "Select Test Data Directory...",
        "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT, VoreenFileWatchListener::ALWAYS_OFF)
    , tempDataPath_("tempDataPath", "Temporary Data Directory", "Select Temporary Data Directory...",
        "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT, VoreenFileWatchListener::ALWAYS_OFF)
    , showSplashScreen_("showSplashScreen", "Show Splash Screen", true, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    , askForSave_("askForSave", "Notify Unsaved Changes", true, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    //, loadLastWorkspaceOnStartup_("loadLastWorkspaceOnStartup", "Load Last Workspace on Startup", true, Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
    , showStartupWizard_("showStartupWizard", "Show Startup Wizard", true, Processor::INVALID_RESULT)
    , useDoubleBuffering_(true)
    , initialized_(false)
    , initializedGL_(false)
    , networkEvaluationRequired_(false)
    , mersenneTwister_()
    , uuidGenerator_(mersenneTwister_)
{
    id_ = guiName;
    guiName_ = guiName;
    app_ = this;
    cmdParser_->setCommandLine(argc, argv);

    // command line options
    cmdParser_->addFlagOption("help,h", CommandLineParser::AdditionalOption, "Print help message");

    cmdParser_->addFlagOption("revision", CommandLineParser::AdditionalOption, "Print revision of binary and exit.");

    cmdParser_->addOption<bool>("logging", CommandLineParser::AdditionalOption,
        "If set to false, logging is disabled entirely (not recommended)"
        /*, enableLogging_->get(), enableLogging_->get() ? "true" : "false" */);

    cmdParser_->addOption<tgt::LogLevel>("logLevel", CommandLineParser::AdditionalOption,
        "Sets the verbosity of the logger \n(debug|info|warning|error|fatal)"
        /*, tgt::Info, loglevelToString(tgt::Info)*/);

    cmdParser_->addOption<bool>("fileLogging", CommandLineParser::AdditionalOption,
        "Enables HTML file logging \n(ignored, if logging is disabled)"
        /*, enableHTMLLogging_->get(), enableHTMLLogging_->get() ? "true" : "false" */);

    cmdParser_->addOption<std::string>("logFile", CommandLineParser::AdditionalOption,
        "Specifies the HTML log file"
        /*, htmlLogFile_->get() */);

    cmdParser_->addOption<bool>("useCaching", CommandLineParser::AdditionalOption,
        "Enables or disables data caching. Overrides the setting stored in the application settings.");

    cmdParser_->addOption("glslVersion", overrideGLSLVersion_, CommandLineParser::AdditionalOption,
        "Overrides the detected GLSL version (1.10|1.20|1.30|1.40|1.50|3.30|4.00|..)",
        overrideGLSLVersion_);

    cmdParser_->addOption<bool>("useDoubleBuffering", CommandLineParser::AdditionalOption,
        "Enables or disables double buffering.");

    // caching properties
    resetCachePath_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::resetCachePath));
    deleteCache_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::deleteCache));
    addProperty(useCaching_);
    addProperty(volumeCacheLimit_);
    addProperty(octreeCacheLimit_);
    addProperty(cachePath_);
    addProperty(resetCachePath_);
    addProperty(deleteCache_);
    useCaching_.setGroupID("caching");
    volumeCacheLimit_.setGroupID("caching");
    octreeCacheLimit_.setGroupID("caching");
    cachePath_.setGroupID("caching");
    resetCachePath_.setGroupID("caching");
    deleteCache_.setGroupID("caching");
    setPropertyGroupGuiName("caching", "Data Caching");

    // CPU RAM properties
    uint64_t memSize = MemoryInfo::getTotalPhysicalMemory() / (1024 * 1024);
    if (!memSize) {
        // could not determine RAM size -> cannot be set automatically
        determineCpuRamLimitAutomatically_.set(false);
        determineCpuRamLimitAutomatically_.setReadOnlyFlag(true);
    }
    cpuRamLimit_.setGroupID("cpu-ram");
    cpuRamLimit_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::cpuRamLimitChanged));
    determineCpuRamLimitAutomatically_.setGroupID("cpu-ram");
    determineCpuRamLimitAutomatically_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::determineRamLimitAutomaticallyChanged));
    addProperty(determineCpuRamLimitAutomatically_);
    addProperty(cpuRamLimit_);
    setPropertyGroupGuiName("cpu-ram", "CPU Memory");

    // gpu memory properties
    determineGpuMemoryLimitAutomatically_.setGroupID("graphicsMemory");
    gpuMemoryLimit_.setGroupID("graphicsMemory");
    determineGpuMemoryLimitAutomatically_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::determineGpuMemoryLimitAutomaticallyChanged));
    gpuMemoryLimit_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::gpuMemoryLimitChanged));
    addProperty(determineGpuMemoryLimitAutomatically_);
    addProperty(gpuMemoryLimit_);
    availableGraphicsMemory_.setReadOnlyFlag(true);
    addProperty(availableGraphicsMemory_);
    refreshAvailableGraphicsMemory_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::queryAvailableGraphicsMemory));
    addProperty(refreshAvailableGraphicsMemory_);
    availableGraphicsMemory_.setGroupID("graphicsMemory");
    refreshAvailableGraphicsMemory_.setGroupID("graphicsMemory");
    setPropertyGroupGuiName("graphicsMemory", "Graphics Memory");

    // logging properties
    addProperty(enableLogging_);
    enableLogging_.setVisibleFlag(false);
    logLevel_.addOption("debug",   "Debug",    tgt::Debug);
    logLevel_.addOption("info",    "Info",     tgt::Info);
    logLevel_.addOption("warning", "Warning",  tgt::Warning);
    logLevel_.addOption("error",   "Error",    tgt::Error);
    logLevel_.addOption("fatal",   "Fatal",    tgt::Fatal);
    logLevel_.selectByKey("info");
    logLevel_.setDefaultValue("info");
    addProperty(logLevel_);
    addProperty(enableHTMLLogging_);
    addProperty(htmlLogFile_);

    logLevel_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::logLevelChanged));

    enableLogging_.setGroupID("logging");
    logLevel_.setGroupID("logging");
    enableHTMLLogging_.setGroupID("logging");
    htmlLogFile_.setGroupID("logging");
    setPropertyGroupGuiName("logging", "Logging");

    // regression test properties
    addProperty(testDataPath_);
    testDataPath_.setGroupID("regressionTesting");
    setPropertyGroupGuiName("regressionTesting", "Regression Testing");

    //  temporary data properties
    addProperty(tempDataPath_);
    tempDataPath_.onChange(MemberFunctionCallback<VoreenApplication>(this, &VoreenApplication::tempDataPathChanged));
    tempDataPath_.setGroupID("tempData");
    setPropertyGroupGuiName("tempData", "Temporary Data");

    // user interface properties
    addProperty(showSplashScreen_);
    showSplashScreen_.setGroupID("user-interface");
    addProperty(askForSave_);
    askForSave_.setGroupID("user-interface");
    //addProperty(loadLastWorkspaceOnStartup_);
    //loadLastWorkspaceOnStartup_.setGroupID("user-interface");
    addProperty(showStartupWizard_);
    showStartupWizard_.setGroupID("user-interface");
    setPropertyGroupGuiName("user-interface", "User Interface");
}

VoreenApplication::~VoreenApplication() {
    if (initializedGL_) {
        std::cerr << "~VoreenApplication(): OpenGL deinitialization has not been performed. Call deinitializeGL() before destruction.\n";
        return;
    }

    if (initialized_) {
        std::cerr << "~VoreenApplication(): application has not been deinitialized. Call deinitialize() before destruction.\n";
        return;
    }
}

VoreenApplication* VoreenApplication::app() {
    return app_;
}

const std::string& VoreenApplication::getBinaryName() const {
    return binaryName_;
}

VoreenApplication::ApplicationFeatures VoreenApplication::getApplicationType() const {
    return appFeatures_;
}

void VoreenApplication::setLoggingEnabled(bool enabled) {
    enableLogging_.set(enabled);
}

bool VoreenApplication::isLoggingEnabled() const {
    return enableLogging_.get();
}

void VoreenApplication::setLogLevel(tgt::LogLevel logLevel) {
    logLevel_.selectByValue(logLevel);
}

tgt::LogLevel VoreenApplication::getLogLevel() const {
    return logLevel_.getValue();
}

void VoreenApplication::setFileLoggingEnabled(bool enabled) {
    enableHTMLLogging_.set(enabled);
}

bool VoreenApplication::isFileLoggingEnabled() const {
    return enableHTMLLogging_.get();
}

void VoreenApplication::setLogFile(const std::string& logFile) {
    htmlLogFile_.set(logFile);
}

const std::string& VoreenApplication::getLogFile() const {
    return htmlLogFile_.get();
}

/*bool VoreenApplication::getLoadLastWorkspaceOnStartup() const {
    return loadLastWorkspaceOnStartup_.get();
}*/

bool VoreenApplication::getShowStartupWizard() const {
    return showStartupWizard_.get();
}


bool VoreenApplication::getShowSplashScreen() const {
    return showSplashScreen_.get();
}

bool VoreenApplication::getAskForSave() const {
    return askForSave_.get();
}

void VoreenApplication::setOverrideGLSLVersion(const std::string& version) {
    if (isInitialized()) {
        LERROR("Trying to override GLSL version after application initialization");
    }
    else
        overrideGLSLVersion_ = version;
}

const std::string& VoreenApplication::getOverrideGLSLVersion() const {
    return overrideGLSLVersion_;
}

void VoreenApplication::setModuleLoadingEnabled(bool enabled) {
    if (isInitialized()) {
        LERROR("Trying to change module loading after application initialization");
    }
    else
        loadModules_ = enabled;
}

bool VoreenApplication::isModuleLoadingEnabled() const {
    return loadModules_;
}

void VoreenApplication::setDeploymentMode(bool dm) {
    if (isInitialized()) {
        LERROR("Trying to change deployment mode after application initialization");
    }
    else
        deploymentMode_ = dm;
}

bool VoreenApplication::getDeploymentMode() const {
    return deploymentMode_;
}

CommandLineParser* VoreenApplication::getCommandLineParser() const {
    return cmdParser_.get();
}

CommandQueue* VoreenApplication::getCommandQueue() const {
    return cmdQueue_.get();
}

void VoreenApplication::loadModules() {
    if (isModuleLoadingEnabled()) {
        LDEBUG("Loading modules from module registration header");
        registerAllModules(this); //< included from gen_moduleregistration.h
    }
    else {
        LDEBUG("Module auto loading disabled");
        registerModule(new CoreModule(getBasePath("modules/core")));  //< core module is always included
    }
}

void VoreenApplication::initialize() {

    if (initialized_)
        throw VoreenException("Application already initialized");

    //
    // Path detection
    //

    // detect documents path first, needed for log file
    std::string documentsPath;
#ifdef WIN32
    TCHAR szPath[MAX_PATH];
    // get "my documents" directory
    if (SHGetFolderPath(NULL, CSIDL_PERSONAL, NULL, 0, szPath) == S_OK)
        documentsPath = szPath;
#else
    if (getenv("HOME") != 0) {
        documentsPath = getenv("HOME");
    }
    else {
        std::cerr << "Failed to detect home.";
    }
#endif

    std::string programPath = getProgramPath();

#if defined(__APPLE__) && defined(VRN_DEPLOYMENT)
    basePath_ = findAppBundleResourcesPath() + "/voreenRoot";
#else
    #if defined(VRN_BASE_PATH) // use base path passed by CMAKE, if present

// See https://gcc.gnu.org/onlinedocs/cpp/Stringizing.html:
#define VRN_BASE_PATH_STRINGIFY(x) VRN_BASE_PATH_STRINGIFY_CONTENT(x)
#define VRN_BASE_PATH_STRINGIFY_CONTENT(x) #x
        basePath_ = VRN_BASE_PATH_STRINGIFY(VRN_BASE_PATH);
#undef VRN_BASE_PATH_STRINGIFY
#undef VRN_BASE_PATH_STRINGIFY_CONTENT

        if (!tgt::FileSystem::dirExists(basePath_)) {
            std::cerr << "WARNING: Passed base path does not exist: " << basePath_ << ". Using current directory instead.\n";
            basePath_ = ".";
        }
        else {
            basePath_ = tgt::FileSystem::cleanupPath(basePath_);
        }
    #else // else try to find base path starting at program path

        // detect base path based on program location
        // (program path is available before command line parser execution)
        basePath_ = programPath;
        // If we failed to find the program path for whatever reason we fall back to the current directory
        if(basePath_.empty()) {
            basePath_ = ".";
        }

        // try to find base path starting at program path
        basePath_ = findBasePath(basePath_);
        if (basePath_.empty()) {
            std::cerr << "WARNING: Base path not found. Using current directory instead.\n";
            basePath_ = ".";
        }
        basePath_ = tgt::FileSystem::absolutePath(basePath_);
    #endif
#endif

    // use user directory for custom data, if in deployment mode
//#ifdef WIN32
    userDataPath_ = tgt::FileSystem::cleanupPath(getBasePath("data"));
//#else
    if (getDeploymentMode() && !documentsPath.empty()) {
        char lastChar = *documentsPath.rbegin();
        if (lastChar != '/' && lastChar != '\\')
            documentsPath += "/";
        userDataPath_ = tgt::FileSystem::cleanupPath(documentsPath + "Voreen");
    }
    else // use VRN_HOME/data as data directory
        userDataPath_ = tgt::FileSystem::cleanupPath(getBasePath("data"));
//#endif

    //
    // Execute command line parser
    //
    if (appFeatures_ & APP_COMMAND_LINE) {
        tgtAssert(cmdParser_, "no command line parser");

        if (appFeatures_ & APP_CONFIG_FILE)
            cmdParser_->setConfigFile(getUserDataPath(tgt::FileSystem::baseName(binaryName_) + ".cfg"));

        try {
            cmdParser_->execute();
        }
        catch (VoreenException& e) {
            // a missing required argument causes an exception even if the --help flag has been passed,
            // so check for it here
            bool helpFlag = false;
            cmdParser_->getOptionValue("help", helpFlag);
            if (helpFlag) {
                std::cout << cmdParser_->getUsageString() << std::endl;
                exit(EXIT_SUCCESS);
            }
            else {
                std::cerr << "\nError: " << e.what() << "\n ";
                //std::cerr << cmdParser_->getUsageString(CommandLineParser::AllTypes, false) << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        // if --help has been passed, print usage and exit program
        bool helpFlag = false;
        cmdParser_->getOptionValue("help", helpFlag);
        if (helpFlag) {
            std::cout << cmdParser_->getUsageString() << std::endl;
            exit(EXIT_SUCCESS);
        }

        // if --revision has been passed, print revision and exit
        bool revisionFlag = false;
        cmdParser_->getOptionValue("revision", revisionFlag);
        if (revisionFlag) {
            std::cout << VoreenVersion::getRevision() << std::endl;
            exit(EXIT_SUCCESS);
        }
    }
    else {
        std::cout << "Command line disabled." << std::endl;
    }

    //
    // tgt initialization
    //
    tgt::init(tgt::InitFeature::ALL, tgt::Info);

    // Init logging.
    // Setting properties needs to be done here, since the base path must already be set.
    htmlLogFile_.set(binaryName_ + "-log.html");
    htmlLogFile_.setDefaultValue(binaryName_ + "-log.html");

    std::string absLogPath;
    initLogging(absLogPath);

    // output ascii logo to stdout
    printAsciiLogo();

    // log basic information
    LINFO("Program: " << getBinaryName());
    VoreenVersion::logAll("voreen.VoreenApplication");
    LINFO("Base path:   " << basePath_);
    LINFO("App Resource Path:    " << getApplicationResourcePath());
    LINFO("Program path:    " << tgt::FileSystem::dirName(programPath));
    LINFO("User data path:  " << getUserDataPath() << (getDeploymentMode() ? " (Deployment mode)" : " (Developer mode)"));

    // Show total system main memory
    LINFO("Total physical CPU RAM: " << MemoryInfo::getTotalPhysicalMemory()  / (1024 * 1024) << " MB");

    // load modules
    try {
        LDEBUG("Loading modules");
        loadModules();
    }
    catch (VoreenException& e) {
        LERROR("Failed to load modules: " << e.what());
        throw;
    }
    if (modules_.empty())
        LWARNING("No modules loaded");
    else {
        std::vector<std::string> moduleNames;
        for (size_t i=0; i<modules_.size(); i++)
            moduleNames.push_back(modules_[i]->getID());
        LINFO("Modules: " << strJoin(moduleNames, ", "));
    }

    // load settings
    initApplicationSettings();
    loadApplicationSettings();

    // We need to set a proper initial tmp path.
    if(tempDataPath_.get().empty())
        tempDataPath_.set(getUserDataPath("tmp"));

    // Clear old temporary data from other instances.
    cleanOrphanedTemporaryData();

    // initialize the file watch system
    VoreenFileWatcher::init();

    // initialize the memory manager now that memory settings have been set in the application
    VolumeMemoryManager::init();

    // apply command-line logging parameters, if specified
    if (cmdParser_->isOptionSet("logging")) {
        bool logging;
        cmdParser_->getOptionValue("logging", logging);
        setLoggingEnabled(logging);
    }
    else {
        setLoggingEnabled(true);
    }
    if (cmdParser_->isOptionSet("logLevel")) {
        tgt::LogLevel logLevel = tgt::Info;
        cmdParser_->getOptionValue("logLevel", logLevel);
        setLogLevel(logLevel);
    }
    if (cmdParser_->isOptionSet("fileLogging")) {
        bool fileLogging;
        cmdParser_->getOptionValue("fileLogging", fileLogging);
        setFileLoggingEnabled(fileLogging);
    }
    if (cmdParser_->isOptionSet("logFile")) {
        std::string logFile;
        cmdParser_->getOptionValue("logFile", logFile);
        if (!logFile.empty()) {
            // Resolve non-absolute paths relative to cwd:
            if(!tgt::FileSystem::isAbsolutePath(logFile)) {
                logFile = tgt::FileSystem::cleanupPath(tgt::FileSystem::currentDirectory() + "/" + logFile);
            }
            setLogFile(logFile);
        }
    }

    // reinit logging, if default settings have been overwritten
    if (enableLogging_.get() != enableLogging_.getDefault()         || logLevel_.get() != logLevel_.getDefault()        ||
        enableHTMLLogging_.get() != enableHTMLLogging_.getDefault() || htmlLogFile_.get() != htmlLogFile_.getDefault()  )
    {
        initLogging(absLogPath);
    }
    if (!absLogPath.empty())
        LINFO("HTML log file:  " << absLogPath);

    // override caching setting, if specified on command line
    if (cmdParser_->isOptionSet("useCaching")) {
        bool caching = true;
        cmdParser_->getOptionValue<bool>("useCaching", caching);
        LINFO((caching ? "Enabled" : "Disabled") << " caching via command line");
        setUseCaching(caching);
    }

    if (cmdParser_->isOptionSet("useDoubleBuffering")) {
        bool db = true;
        cmdParser_->getOptionValue<bool>("useDoubleBuffering", db);
        LINFO((db ? "Enabled" : "Disabled") << " double buffering via command line");
        useDoubleBuffering_ = db;
    }

    // initialize modules
    LINFO("Initializing modules");
    for (size_t i=0; i<modules_.size(); i++) {
        try {
            LDEBUG("Initializing module '" << modules_.at(i)->getID() << "'");
            modules_.at(i)->initialize();
            modules_.at(i)->initialized_ = true;
        }
        catch (const VoreenException& e) {
            LERROR("VoreenException during initialization of module '" << modules_.at(i)->getID() << "': " << e.what());
            modules_.at(i)->initialized_ = false;
        }
        catch (const std::exception& e) {
            LERROR("std::exception during initialization of module '" << modules_.at(i)->getID() << "': " << e.what());
            modules_.at(i)->initialized_ = false;
        }
        catch (...) {
            LERROR("Unknown exception during initialization of module '" << modules_.at(i)->getID() << "'");
            modules_.at(i)->initialized_ = false;
        }
    }

    // init timer
    schedulingTimer_ = createTimer(&eventHandler_);
    eventHandler_.addListenerToFront(this);

    // init quality mode
    VoreenQualityMode::init();

    initialized_ = true;
}

void VoreenApplication::deinitialize() {

    if (!initialized_)
        throw VoreenException("Application not initialized");

    if (initializedGL_)
        throw VoreenException("OpenGL deinitialization not performed. Call deinitializeGL() before deinitialization!");

    // Clean cache and temp files.
    cleanCache();
    cleanTemporaryData();

    delete schedulingTimer_;
    schedulingTimer_ = 0;

    saveApplicationSettings();

    // deinitialize modules
    LINFO("Deinitializing modules");
    for (int i=(int)modules_.size()-1; i>=0; i--) {
        try {
            LDEBUG("Deinitializing module '" << modules_.at(i)->getID() << "'");
            modules_.at(i)->deinitialize();
            modules_.at(i)->initialized_ = false;
        }
        catch (const VoreenException& e) {
            LERROR("VoreenException during deinitialization of module '" << modules_.at(i)->getID() << "': " << e.what());
        }
        catch (const std::exception& e) {
            LERROR("std::exception during deinitialization of module '" << modules_.at(i)->getID() << "': " << e.what());
        }
        catch (...) {
            LERROR("Unknown exception during deinitialization of module '" << modules_.at(i)->getID() << "'");
        }
    }

    // clear modules
    LDEBUG("Deleting modules");
    for (int i=(int)modules_.size()-1; i>=0; i--) {
        LDEBUG("Deleting module '" << modules_.at(i)->getID() << "'");
        delete modules_.at(i);
    }
    modules_.clear();

    // deinitialize volume memory manager
    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::deinit();

    // deinitialize the file watch system
    if (VoreenFileWatcher::isInited())
        VoreenFileWatcher::deinit();

    LDEBUG("tgt::deinit()");
    tgt::deinit();

    // deinit quality mode
    VoreenQualityMode::deinit();

    initialized_ = false;
}

void VoreenApplication::initializeGL() {
    if (!initialized_)
        throw VoreenException("Application not initialized. Call initialize() before OpenGL initialization!");

    if (initializedGL_)
        throw VoreenException("OpenGL initialization already done");

    LINFO("Initializing OpenGL");

    LDEBUG("tgt::initGL");
    tgt::initGL();

    if (!overrideGLSLVersion_.empty()) {
        LWARNING("Overriding detected GLSL version " << GpuCaps.getShaderVersion()
            << " with version: " << overrideGLSLVersion_);
        GpuCaps.overrideGLSLVersion(overrideGLSLVersion_);
    }

    GpuCaps.logCapabilities(false, true);

    // now that OpenGL is initialized, we can check the GPU memory settings
    int texMemSize = GpuCaps.retrieveTotalTextureMemory() / 1024; // kilobytes to MB
    if (!(texMemSize > 0)) {
        // could not determine GPU memory -> cannot be set automatically
        determineGpuMemoryLimitAutomatically_.set(false);
        determineGpuMemoryLimitAutomatically_.setReadOnlyFlag(true);
    }
    else {
        // set the maximum, even if the actual value is not set automatically
        gpuMemoryLimit_.setMaxValue(texMemSize);

        if (determineGpuMemoryLimitAutomatically_.get()) {
            // set a "good" default value: take 3/4 of the available GPU memory
            gpuMemoryLimit_.set((texMemSize / 4) * 3);
        }

    }
    gpuMemoryLimit_.setReadOnlyFlag(determineGpuMemoryLimitAutomatically_.get());

    // to load the immediate mode shader
    ShdrMgr.addPath(getBasePath() + "/ext/tgt/glsl");

    //set depth test enable to satisfy voreen guidelines
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // OpenGL initialize modules
    LINFO("OpenGL initializing modules");
    for (size_t i=0; i<modules_.size(); i++) {
        if (!modules_.at(i)->isInitialized()) {
            LERROR("Module '" << modules_.at(i)->getID() << "' has not been initialized before OpenGL initialization");
            modules_.at(i)->initializedGL_ = false;
            continue;
        }
        try {
            LDEBUG("OpenGL initialization of module '" << modules_.at(i)->getID() << "'");
            modules_.at(i)->initializeGL();
            modules_.at(i)->initializedGL_ = true;
        }
        catch (const VoreenException& e) {
            LERROR("VoreenException during OpenGL initialization of module '" << modules_.at(i)->getID() << "': " << e.what());
            modules_.at(i)->initializedGL_ = false;
        }
        catch (const std::exception& e) {
            LERROR("std::exception during OpenGL initialization of module '" << modules_.at(i)->getID() << "': " << e.what());
            modules_.at(i)->initializedGL_ = false;
        }
        catch (...) {
            LERROR("Unknown exception during OpenGL initialization of module '" << modules_.at(i)->getID() << "'");
            modules_.at(i)->initializedGL_ = false;
        }
    }

    queryAvailableGraphicsMemory();

    initializedGL_ = true;
}

void VoreenApplication::deinitializeGL() {

    if (!initializedGL_)
        throw VoreenException("OpenGL not initialized");

    LINFO("OpenGL deinitializing modules");
    for (int i=(int)modules_.size()-1; i>=0; i--) {

        if (modules_.at(i)->isInitializedGL()) {
            try {
                LDEBUG("OpenGL deinitialization of module '" << modules_.at(i)->getID() << "'");
                modules_.at(i)->deinitializeGL();
                modules_.at(i)->initializedGL_ = false;
            }
            catch (const VoreenException& e) {
                LERROR("VoreenException during OpenGL deinitialization of module '" << modules_.at(i)->getID() << "': " << e.what());
            }
            catch (const std::exception& e) {
                LERROR("std::exception during OpenGL deinitialization of module '" << modules_.at(i)->getID() << "': " << e.what());
            }
            catch (...) {
                LERROR("unknown exception during OpenGL deinitialization of module '" << modules_.at(i)->getID() << "'");
            }
        }
        else {
            LWARNING("Skipping OpenGL deinitialization of module '" << modules_.at(i)->getID() << "': not OpenGL initialized");
        }
    }

    LDEBUG("tgt::deinitGL");
    tgt::deinitGL();

    initializedGL_ = false;
}

bool VoreenApplication::isInitialized() const {
    return initialized_;
}

bool VoreenApplication::isInitializedGL() const {
    return initializedGL_;
}

bool VoreenApplication::useCaching() const {
    return useCaching_.get();
}

void VoreenApplication::setUseCaching(bool useCaching) {
    useCaching_.set(useCaching);
}

bool VoreenApplication::useDoubleBuffering() const {
    return useDoubleBuffering_;
}

void VoreenApplication::setUseDoubleBuffering(bool useDoubleBuffering) {
    useDoubleBuffering_ = useDoubleBuffering;
}

size_t VoreenApplication::getCpuRamLimit() const {
    if (cpuRamLimit_.get() > 0)
        return static_cast<size_t>(static_cast<uint64_t>(cpuRamLimit_.get()) << 20); //< property specifies the limit in MB
    else
        return static_cast<size_t>(static_cast<uint64_t>(cpuRamLimit_.getMaxValue()) << 20); //< use max value
}

size_t VoreenApplication::getGpuMemoryLimit() const {
    return static_cast<size_t>(static_cast<uint64_t>(gpuMemoryLimit_.get()) << 20); //< property specifies limit in MB
}

void VoreenApplication::setCpuRamLimit(size_t ramLimit) {
    cpuRamLimit_.set(static_cast<int>(ramLimit >> 20)); //< property specifies the limit in MB
}

void VoreenApplication::registerModule(VoreenModule* module) {
    tgtAssert(module, "null pointer passed");

    // check if module's name and dirName have been set
    if (module->getID().empty() || module->getID() == "undefined") {
        tgtAssert(false, "module has no name (set in constructor!)");
        LERROR("Module has no name (set in constructor!). Skipping.");
        return;
    }

    // check if module directory exists
    /*std::string moduleDir = module->getModulePath();
    if (!tgt::FileSystem::dirExists(moduleDir))
        LWARNING("Module '" << module->getID() << "': module directory '" << moduleDir << "' does not exist"); */

    // register module
    if (std::find(modules_.begin(), modules_.end(), module) == modules_.end())
        modules_.push_back(module);
    else
        LWARNING("Module '" << module->getID() << "' has already been registered. Skipping.");
}

const std::vector<VoreenModule*>& VoreenApplication::getModules() const {
    return modules_;
}

VoreenModule* VoreenApplication::getModule(const std::string& moduleName) const {
    // search by module name first
    for (size_t i = 0 ; i < modules_.size() ; ++i) {
        VoreenModule* module = modules_.at(i);
        if (module->getID() == moduleName)
            return module;
    }

    // then search by module directory name
    for (size_t i = 0 ; i < modules_.size() ; ++i) {
        VoreenModule* module = modules_.at(i);
        if (module->getDirName() == moduleName)
            return module;
    }

    return 0;
}

void VoreenApplication::registerSerializerFactory(SerializableFactory* factory) {
    tgtAssert(factory, "null pointer passed");
    if (std::find(serializerFactories_.begin(), serializerFactories_.end(), factory) == serializerFactories_.end())
        serializerFactories_.push_back(factory);
    else
        LWARNING("SerializerFactory already registered. Skipping.");
}

const std::vector<SerializableFactory*>& VoreenApplication::getSerializerFactories() const {
    return serializerFactories_;
}

ProcessorWidget* VoreenApplication::createProcessorWidget(Processor* processor) const {
    tgtAssert(processor, "null pointer passed");
    if ((appFeatures_ & APP_PROCESSOR_WIDGETS) == 0)
        return 0;

    const std::vector<VoreenModule*>& modules = getModules();
    for (size_t m=0; m<modules.size(); m++) {
        const std::vector<ProcessorWidgetFactory*>& factories = modules_.at(m)->getRegisteredProcessorWidgetFactories();
        for (size_t f=0; f<factories.size(); f++) {
            ProcessorWidget* processorWidget = factories.at(f)->createWidget(processor);
            if (processorWidget)
                return processorWidget;
        }
    }
    return 0;
}

PropertyWidget* VoreenApplication::createPropertyWidget(Property* property) const {
    tgtAssert(property, "null pointer passed");
    if ((appFeatures_ & APP_PROPERTY_WIDGETS) == 0)
        return 0;

    const std::vector<VoreenModule*>& modules = getModules();
    for (size_t m=0; m<modules.size(); m++) {
        const std::vector<PropertyWidgetFactory*>& factories = modules_.at(m)->getRegisteredPropertyWidgetFactories();
        for (size_t f=0; f<factories.size(); f++) {
            PropertyWidget* propertyWidget = factories.at(f)->createWidget(property);
            if (propertyWidget)
                return propertyWidget;
        }
    }
    return 0;
}

tgt::Timer* VoreenApplication::createTimer(tgt::EventHandler* /*handler*/) const {
    return 0;
}

ProgressBar* VoreenApplication::createProgressDialog() const {
    return 0;
}

void VoreenApplication::showMessageBox(const std::string& /*title*/, const std::string& message, bool /*error=false*/) const {
    LINFO("showMessageBox() not implemented: message=" << message);
}

std::string VoreenApplication::getBasePath(const std::string& filename) const {
    return tgt::FileSystem::cleanupPath(basePath_ + (filename.empty() ? "" : "/" + filename));
}

boost::uuids::uuid VoreenApplication::generateUUID() const {
    return uuidGenerator_();
}

std::string VoreenApplication::getUniqueFilePath(const std::string& root, const std::string& suffix) const {
    std::string outputPath = "";
    do {
        std::string name = boost::lexical_cast<std::string>(generateUUID()) + suffix;
        outputPath = tgt::FileSystem::cleanupPath(root + "/" + name);
    } while (tgt::FileSystem::fileExists(outputPath) || tgt::FileSystem::dirExists(outputPath));
    return outputPath;
}

std::string VoreenApplication::getCachePath(const std::string& filename) const {
    std::string cacheBasePath;

    // use property value, if set and directory does exist
    if (!cachePath_.get().empty()) {
        cacheBasePath = cachePath_.get();
        if (!tgt::FileSystem::dirExists(cacheBasePath)) {
            LWARNING("Cache path does not exist: " << cacheBasePath << ". Switching to user data path: " << getUserDataPath("cache"));
            cacheBasePath = "";
        }
    }
    // otherwise use user data path
    if (cacheBasePath.empty())
        cacheBasePath = getUserDataPath("cache");

    return tgt::FileSystem::cleanupPath(cacheBasePath + (filename.empty() ? "" : "/" + filename));
}

std::string VoreenApplication::getProgramPath() const {
    tgtAssert(cmdParser_, "no command line parser");
    // Resolve symlinks and thus (hopefully) get the actual path to the voreen binary
    std::string binaryPath = tgt::FileSystem::absolutePath(cmdParser_->getProgramPath());
    return tgt::FileSystem::dirName(binaryPath);
}

void VoreenApplication::cleanCache() const {
    // clean volume cache
    int volumeCacheLimitMB = volumeCacheLimit_.get() >> 10; //< property specifies GB
    CacheCleaner cleaner;
    cleaner.initialize(getCachePath());
    cleaner.deleteUnused();
    cleaner.limitCache(volumeCacheLimitMB);

    // clean octree cache
    uint64_t octreeCacheLimitBytes = (uint64_t)octreeCacheLimit_.get() << 30;//< property specifies GB
    OctreeCreator::limitCacheSize(getCachePath(), octreeCacheLimitBytes, false);
}

void VoreenApplication::deleteCache() {
    std::string cacheDir = getCachePath();

    std::vector<std::string> files = FileSys.listFiles(cacheDir);
    for(size_t i=0; i<files.size(); i++)
        FileSys.deleteFile(cacheDir+"/"+files[i]);

    std::vector<std::string> dirs = FileSys.listSubDirectories(cacheDir);
    for(size_t i=0; i<dirs.size(); i++)
        FileSys.deleteDirectoryRecursive(cacheDir+"/"+dirs[i]);
}

void VoreenApplication::resetCachePath() {
    cachePath_.set("");
}

std::string VoreenApplication::getCoreResourcePath(const std::string& filename) const {
    return tgt::FileSystem::cleanupPath(getBasePath() + "/resource/voreencore" + (filename.empty() ? "" : "/" + filename));
}

std::string VoreenApplication::getApplicationResourcePath(const std::string& filename) const {
    return getCoreResourcePath(filename);
}

std::string VoreenApplication::getUserDataPath(const std::string& filename) const {
    return tgt::FileSystem::cleanupPath(userDataPath_ + (filename.empty() ? "" : "/" + filename));
}

std::string VoreenApplication::getFontPath(const std::string& filename) const {
    return tgt::FileSystem::cleanupPath(getCoreResourcePath("fonts") + (filename.empty() ? "" : "/" + filename));
}

std::string VoreenApplication::getModulePath(const std::string& moduleName) const {
    tgtAssert(moduleName != "", "empty module name passed");
    if (getModule(moduleName))
        return getModule(moduleName)->getModulePath();
    else
        return "";
}

std::string VoreenApplication::getTemporaryPath(const std::string& filename) const {
    return tgt::FileSystem::cleanupPath(tempDataPathInstance_ + (filename.empty() ? "" : "/" + filename));
}

std::string VoreenApplication::getUniqueTmpFilePath(const std::string& suffix) const {
    return getUniqueFilePath(tempDataPathInstance_, suffix);
}

void VoreenApplication::cleanTemporaryData() {
    if(tgt::FileSystem::dirExists(tempDataPathInstance_)) {
        tempDataPathLock_.reset(nullptr);
        if (tgt::FileSystem::deleteDirectoryRecursive(tempDataPathInstance_)) {
            LINFO("Successfully removed temporary data");
        }
        else {
            LWARNING("Failed to delete temporary data");
        }
    }
}

void VoreenApplication::cleanOrphanedTemporaryData() {
    std::vector<std::string> directories = tgt::FileSystem::listSubDirectories(tempDataPath_.get());
    for (const std::string& dir : directories) {
        std::string dirPath = tgt::FileSystem::cleanupPath(tempDataPath_.get() + "/" + dir);

        // Skip current instance.
        if (dirPath == tempDataPathInstance_)
            continue;

        std::vector<std::string> files = tgt::FileSystem::listFiles(dirPath);
        for (const std::string& file : files) {
            if (file == VOREEN_LOCK_NAME) {
                // Figure out if the instance is still running.
                std::string filePath = tgt::FileSystem::cleanupPath(dirPath + "/" + file);
                try {
                    if(boost::interprocess::file_lock(filePath.c_str()).try_lock()) {
                        if (tgt::FileSystem::deleteDirectoryRecursive(dirPath)) {
                            LINFO("Successfully removed temporary data of instance: " << dir);
                        }
                        else {
                            LWARNING("Failed to delete temporary data of instance " << dir);
                        }
                    }
                }
                catch (const boost::interprocess::interprocess_exception&) {
                    // Ignore silently.
                }

                // Go to next directory.
                break;
            }
        }
    }
}

void VoreenApplication::tempDataPathChanged() {

    // Do NOT delete old temporary data, since it might be still in use.
    //cleanTemporaryData();
    // But DO release the lock.
    tempDataPathLock_.reset(nullptr);

    // Use property value, if possible.
    std::string tempBasePath;
    if (!tempDataPath_.get().empty()) {
        tempBasePath = tempDataPath_.get();
        if (!tgt::FileSystem::dirExists(tempBasePath)) {
            LWARNING("Temporary path does not exist: " << tempBasePath << ". Switching to user data path: " << getUserDataPath("tmp"));
            tempBasePath = "";
        }
    }
    // Otherwise use user data path.
    if (tempBasePath.empty())
        tempBasePath = getUserDataPath("tmp");

    // Update property.
    tempDataPath_.set(tempBasePath);

    // Set temporary data path for this instance.
    tempDataPathInstance_ = getUniqueFilePath(tempBasePath);

    // Check for temporary directory.
    if (!tgt::FileSystem::dirExists(tempDataPathInstance_)) {
        if (!tgt::FileSystem::createDirectoryRecursive(tempDataPathInstance_)) {
            throw VoreenException("Failed to create temporary directory");
        }
        else
            LINFO("Created temporary directory: " << tempDataPathInstance_);
    }

    // try to write to temporary directory (we can assume it exists)
    std::string file_lock = tgt::FileSystem::cleanupPath(tempDataPathInstance_ + "/" + VOREEN_LOCK_NAME);
    std::ifstream fs;
    fs.open(file_lock, std::fstream::in | std::fstream::out | std::fstream::trunc);
    if (!fs.is_open()) {
        throw VoreenException("Failed to write to temporary directory");
    }
    fs.close();

    // We can assume the directory exists, so create a file lock inside.
    try {
        tempDataPathLock_.reset(new boost::interprocess::file_lock(file_lock.c_str()));
        tempDataPathLock_->lock();
    }
    catch (const boost::interprocess::interprocess_exception& exception) {
        throw VoreenException(std::string("Failed to create file lock: ") + exception.what());
    }
}

void VoreenApplication::setTempDataPath(const std::string& path) {
    tempDataPath_.set(path);
}

std::string VoreenApplication::getTestDataPath() const {
    return testDataPath_.get();
}

void VoreenApplication::setTestDataPath(const std::string& path) {
    testDataPath_.set(path);
}

void VoreenApplication::scheduleNetworkProcessing() {
    if (schedulingTimer_ && !networkEvaluators_.empty() /*&& schedulingTimer_->isStopped()*/) {
        // schedule network for immediate re-evaluation
        networkEvaluationRequired_ = true;
        if (!schedulingTimer_->isStopped())
            schedulingTimer_->stop();

        // This line used to be indented in line with the above if statement
        // which obviously had no effect, but could be what we want...
        schedulingTimer_->start(0, 1);
    }
}

void VoreenApplication::timerEvent(tgt::TimeEvent* /*e*/) {
    if (!isInitialized())
        return;

    // Execute Command queue.
    cmdQueue_->executeAll();

    if (networkEvaluationRequired_) {
        networkEvaluationRequired_ = false;
        for (std::set<NetworkEvaluator*>::iterator it = networkEvaluators_.begin(); it != networkEvaluators_.end(); it++) {
            if (!(*it)->isLocked()) 
                (*it)->process();
        }
    }

    // check every 100 ms if the network has to be re-evaluated,
    // unless a re-evaluation has been explicitly scheduled via scheduleNetworkProcessing
    if (schedulingTimer_->isStopped() || schedulingTimer_->getTickTime() != 100) {
        schedulingTimer_->stop();
        schedulingTimer_->start(100, 0);
    }
}

void VoreenApplication::queryAvailableGraphicsMemory() {
    if (tgt::GpuCapabilities::isInited()) {
        int mem = (GpuCaps.retrieveAvailableTextureMemory() != -1 ? (GpuCaps.retrieveAvailableTextureMemory() >> 10) : -1);
        // do a check to avoid a warning when the property is set
        if (mem < availableGraphicsMemory_.getMinValue()) {
            availableGraphicsMemory_.setMinValue(mem);
            availableGraphicsMemory_.setMaxValue(mem);
        }
        else {
            availableGraphicsMemory_.setMaxValue(mem);
            availableGraphicsMemory_.setMinValue(mem);
        }
        availableGraphicsMemory_.setDefaultValue(mem);
        availableGraphicsMemory_.set(mem);
    }
    else {
        LWARNING("queryAvailableGraphicsMemory(): GpuCapabilities not instantiated");
    }
}

void VoreenApplication::registerNetworkEvaluator(NetworkEvaluator* evaluator) {
    networkEvaluators_.insert(evaluator);
}

void VoreenApplication::deregisterNetworkEvaluator(NetworkEvaluator* evaluator) {
    networkEvaluators_.erase(evaluator);
}

NetworkEvaluator* VoreenApplication::getNetworkEvaluator() const {
    if(networkEvaluators_.empty())
        return 0;
    else
        return *(networkEvaluators_.begin());
}

NetworkEvaluator* VoreenApplication::getNetworkEvaluator(Processor* p) const {
    tgtAssert(p, "null pointer passed");
    for(std::set<NetworkEvaluator*>::iterator it = networkEvaluators_.begin(); it != networkEvaluators_.end(); it++) {
        std::vector<Processor*> pv = (*it)->getProcessorNetwork()->getProcessors();
        for(size_t i = 0; i < pv.size(); i++) {
            if(pv.at(i) == p)
                return *it;
        }
    }
    return 0;
}

NetworkEvaluator* VoreenApplication::getNetworkEvaluator(ProcessorNetwork* network) const {
    tgtAssert(network, "null pointer passed");
    for(std::set<NetworkEvaluator*>::iterator it = networkEvaluators_.begin(); it != networkEvaluators_.end(); it++) {
        if ((*it)->getProcessorNetwork() == network)
            return (*it);
    }
    return 0;
}

void VoreenApplication::initApplicationSettings() {
}

void VoreenApplication::determineRamLimitAutomaticallyChanged() {
    // determine RAM size
    uint64_t memSize = MemoryInfo::getTotalPhysicalMemory() / (1024 * 1024);
    if (!memSize) {
        // could not determine RAM size -> cannot be set automatically
        determineCpuRamLimitAutomatically_.set(false);
        determineCpuRamLimitAutomatically_.setReadOnlyFlag(true);
        return;
    }

    cpuRamLimit_.setReadOnlyFlag(determineCpuRamLimitAutomatically_.get());

    if (determineCpuRamLimitAutomatically_.get()) {
        // set a "good" default value... take 3/4 of the available RAM
        cpuRamLimit_.set(static_cast<int>((memSize / 4) * 3));
    }
}

void VoreenApplication::determineGpuMemoryLimitAutomaticallyChanged() {
    // do nothing if OpenGL has not been properly initialized
    if (!tgt::GpuCapabilities::isInited())
        return;

    // try to determine the GPU memory automatically
    int memSize = GpuCaps.retrieveTotalTextureMemory() / 1024; // kilobytes to MB
    if (!(memSize > 0)) {
        // could not determine GPU memory -> cannot be set automatically
        determineGpuMemoryLimitAutomatically_.set(false);
        determineGpuMemoryLimitAutomatically_.setReadOnlyFlag(true);
        return;
    }

    gpuMemoryLimit_.setReadOnlyFlag(determineGpuMemoryLimitAutomatically_.get());

    if (determineGpuMemoryLimitAutomatically_.get()) {
        // set a "good" default value: take 3/4 of the available GPU memory
        gpuMemoryLimit_.set((memSize / 4) * 3);
    }
}

void VoreenApplication::cpuRamLimitChanged() {
    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().updateMainMemory();
}

void VoreenApplication::gpuMemoryLimitChanged() {
    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().updateGraphicsMemory();
}

void VoreenApplication::loadApplicationSettings() {
    std::string filename = getUserDataPath(toLower(binaryName_) + "_settings.xml");
    if (!deserializeSettings(this, filename)) {
        // try old voreensettings.xml
        if (!deserializeSettings(this, getUserDataPath("voreensettings.xml")))
            LWARNING("Failed to deserialize application settings");
    }

    // check CPU RAM properties, which might have changed (e.g., when upgrading RAM)
    uint64_t mainMemSize = MemoryInfo::getTotalPhysicalMemory() / (1024 * 1024); // bytes to MB
    if (!mainMemSize) {
        // could not determine RAM size -> cannot be set automatically
        determineCpuRamLimitAutomatically_.set(false);
        determineCpuRamLimitAutomatically_.setReadOnlyFlag(true);
    }
    else if (determineCpuRamLimitAutomatically_.get()) {
        // set a "good" default value... take 3/4 of the available RAM
        cpuRamLimit_.set(static_cast<int>((mainMemSize / 4) * 3));
    }
    cpuRamLimit_.setReadOnlyFlag(determineCpuRamLimitAutomatically_.get());

    // do not set the GPU properties here, since the settings are loaded before initializeGL() is called - settings are checked there!
    // check GPU memory properties, which might have changed (e.g., after installing a new graphics card)
    /*int texMemSize = GpuCaps.retrieveTotalTextureMemory() / 1024; // kilobytes to MB
    if (!(texMemSize > 0)) {
        // could not determine GPU memory -> cannot be set automatically
        determineGpuMemoryLimitAutomatically_.set(false);
        determineGpuMemoryLimitAutomatically_.setReadOnlyFlag(true);
    }
    else {
        // set the maximum, even if the actual value is not set automatically
        gpuMemoryLimit_.setMaxValue(texMemSize);

        if (determineGpuMemoryLimitAutomatically_.get()) {
            // set a "good" default value: take 3/4 of the available GPU memory
            gpuMemoryLimit_.set((texMemSize / 4) * 3);
        }

    }
    gpuMemoryLimit_.setReadOnlyFlag(determineGpuMemoryLimitAutomatically_.get());*/

    const std::vector<VoreenModule*>& modules = getModules();
    for(size_t i=0; i<modules.size(); i++) {
        if(!modules[i]->getProperties().empty()) {
            deserializeSettings(modules[i], getUserDataPath(toLower(modules[i]->getID()) + "_settings.xml"));
        }
    }
}

void VoreenApplication::saveApplicationSettings() {
    std::string filename = getUserDataPath(toLower(binaryName_) + "_settings.xml");

    if(!serializeSettings(this, filename))
        LWARNING("Failed to save application settings");

    const std::vector<VoreenModule*>& modules = getModules();
    for(size_t i=0; i<modules.size(); i++) {
        if(!modules[i]->getProperties().empty()) {
            serializeSettings(modules[i], getUserDataPath(toLower(modules[i]->getID()) + "_settings.xml"));
        }
    }
}

void VoreenApplication::serialize(Serializer& s) const {
    PropertyOwner::serialize(s);
}

void VoreenApplication::deserialize(Deserializer& s) {
    PropertyOwner::deserialize(s);
}

std::string VoreenApplication::getSerializableTypeString(const std::type_info& type) const {
    for (size_t i=0; i<modules_.size(); i++) {
        std::string typeString = modules_[i]->getSerializableTypeString(type);
        if (!typeString.empty())
            return typeString;
    }
    return "";
}

VoreenSerializableObject* VoreenApplication::createSerializableType(const std::string& typeString) const {
    const VoreenSerializableObject* type = getSerializableType(typeString);
    if (type)
        return type->create();
    else
        return 0;
}

const VoreenSerializableObject* VoreenApplication::getSerializableType(const std::string& typeString) const {
    for (size_t i=0; i<modules_.size(); i++) {
        const VoreenSerializableObject* instance = modules_[i]->getSerializableType(typeString);
        if (instance)
            return instance;
    }
    return 0;
}

void VoreenApplication::initLogging(std::string& htmlLogFile) {
    if (!tgt::Singleton<tgt::LogManager>::isInited()) {
        std::cerr << "LogManager not initialized";
        return;
    }

    // extract current HTML file path
    std::string previousLogFile;
    std::vector<tgt::Log*> previousLogs = LogMgr.getLogs();
    for (size_t i=0; i<previousLogs.size(); i++)
        if (dynamic_cast<const tgt::HtmlLog*>(previousLogs.at(i)))
            previousLogFile = static_cast<const tgt::HtmlLog*>(previousLogs.at(i))->getAbsFilename();

    LogMgr.clear();

    if (!enableLogging_.get()) {
        std::cout << "Logging disabled!" << std::endl;
        return;
    }

    // Console log
    tgt::ConsoleLog* log = new tgt::ConsoleLog();
    log->addCat("", true, logLevel_.getValue());
    LogMgr.addLog(log);

    // HTML file logging
    if (!htmlLogFile_.get().empty()) {
        std::string absLogPath;

        // add HTML file logger
        tgt::Log* log = 0;
        htmlLogFile = htmlLogFile_.get();
        if (tgt::FileSystem::isAbsolutePath(htmlLogFile)) {
            absLogPath = htmlLogFile;
        }
        else {
            LogMgr.reinit(getUserDataPath()); //< write log file to user data dir by default
            absLogPath = tgt::FileSystem::absolutePath(LogMgr.getLogDir() + "/" + htmlLogFile);
        }
        log = new tgt::HtmlLog(htmlLogFile, false, true, true, true, tgt::FileSystem::cleanupPath(absLogPath) == tgt::FileSystem::cleanupPath(previousLogFile));
        tgtAssert(log, "no log");

        log->addCat("", true, logLevel_.getValue());
        LogMgr.addLog(log);

        htmlLogFile = absLogPath;

        //std::cout << "HTML Log File: " << htmlLogFile << "\n";
    }
    else {
        // should be never be empty (neither default value nor cmd line param)
        LERROR("HTML log file path is empty.");
    }

}

void VoreenApplication::logLevelChanged() {
    if (!initialized_)
        return;

    if (!tgt::Singleton<tgt::LogManager>::isInited()) {
        std::cerr << "LogManager not initialized";
        return;
    }

    LogMgr.setLogLevel(logLevel_.getValue());

}

} // namespace
