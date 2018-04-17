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

#ifndef VRN_APPLICATION_H
#define VRN_APPLICATION_H

#include "voreen/core/voreenmodule.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/propertyowner.h"
#include "voreen/core/io/serialization/serializablefactory.h"

#include <string>
#include <boost/uuid/uuid_generators.hpp>
#include "tgt/logmanager.h"
#include "tgt/event/eventlistener.h"
#include "tgt/event/eventhandler.h"
#include "tgt/timer.h"


namespace tgt {
    class EventHandler;
    class GLCanvas;
}

namespace voreen {

class Processor;
class ProcessorWidget;
class PropertyWidget;
class VoreenModule;
class NetworkEvaluator;
class ProcessorNetwork;
class ProgressBar;
class CommandLineParser;
class CommandQueue;

/**
 * Represents basic properties of a Voreen application. There should only be one instance of
 * this class, which can be access via the static method app().
 */
class VRN_CORE_API VoreenApplication : private tgt::EventListener, public PropertyOwner, public SerializableFactory {

    friend class Processor;
    friend class NetworkEvaluator;

public:
    /// Features used in this application
    enum ApplicationFeatures {
        APP_NONE                =  0,       ///< nothing
        APP_COMMAND_LINE        =  1,       ///< enable command line parsing
        APP_CONFIG_FILE         =  2,       ///< read program parameters from ini file <binaryName>.cfg in user data path
        APP_PROCESSOR_WIDGETS   =  4,       ///< activate creation of processor widgets
        APP_PROPERTY_WIDGETS    =  8,       ///< activate creation of property widgets
        APP_WIDGETS             =  APP_PROCESSOR_WIDGETS | APP_PROPERTY_WIDGETS, ///< activate creation of widgets
        APP_ALL                 =  0xFFFF,  ///< all features
        APP_DEFAULT = APP_ALL               ///< default: all features
    };

    /**
     * @param binaryName Name of the application binary ("voreenve")
     * @param displayName Nice-looking name of the application ("VoreenVE")
     * @param description Description of the application to be printed in the help message
     * @param argc Number of arguments as retrieved from main()
     * @param argv Argument vector as retrieved from main()
     * @param appType Features to activate
     */
    VoreenApplication(const std::string& binaryName, const std::string& guiName, const std::string& description,
                      int argc, char** argv, ApplicationFeatures appType = APP_DEFAULT);

    virtual ~VoreenApplication();

    /// Allows access to the global instance of this class.
    static VoreenApplication* app();

    virtual std::string getClassName() const { return "VoreenApplication"; }

    /// Do not use!
    virtual VoreenApplication* create() const { return nullptr; }

    /// Returns the name of the application binary.
    const std::string& getBinaryName() const;

    /// Returns the application features as passed to the constructor.
    ApplicationFeatures getApplicationType() const;

    //
    // Command line options
    //

    /**
     * If set to false, console and HTML file logging are disabled entirely.
     * May be overridden by command line option.
     * Default: true
     *
     * @note Can only be set before application initialization.
     */
    void setLoggingEnabled(bool enabled);

    /// Returns whether logging is enabled.
    bool isLoggingEnabled() const;

    /// returns the state of the property with the same name
    //bool getLoadLastWorkspaceOnStartup() const;

    /// returns the state of the property with the same name
    bool getShowStartupWizard() const;

    /// returns the state of the property with the same name
    bool getShowSplashScreen() const;

    /// returns the state of the property with the same name
    bool getAskForSave() const;


    /**
     * Sets the log level to be used for all loggers.
     * May be overridden by command line option.
     * Default: Info
     *
     * @note Can only be set before application initialization.
     */
    void setLogLevel(tgt::LogLevel logLevel);

    /// Returns the log level used by all loggers.
    tgt::LogLevel getLogLevel() const;

    /**
     * Specifies whether logging to an HTML file is enabled.
     * May be overridden by command line option.
     * If logging is disabled entirely, this setting is ignored.
     * Default: true
     *
     * @see setLogFile
     *
     * @note Can only be set before application initialization.
     */
    void setFileLoggingEnabled(bool enabled);

    /// Returns whether logging to an HTML file is enabled.
    bool isFileLoggingEnabled() const;

    /**
     * Sets the path to the HTML log file. Setting the log file
     * to empty disables HTML logging.  May be overridden by command line option.
     * Default: <binaryName>-log.html
     *
     * @see setFileLoggingEnabled
     *
     * @note Can only be set before application initialization.
     */
    void setLogFile(const std::string& logFile);

    /// Returns the path to the HTML log file.
    const std::string& getLogFile() const;

    /**
     * Overrides the detected GLSL version with the passed one.
     * May be overridden by command line option.
     * Default: not set
     *
     * @note Can only be set before application initialization.
     */
    void setOverrideGLSLVersion(const std::string& version);

    /// Returns the overriding GLSL version if set, or an empty string otherwise.
    const std::string& getOverrideGLSLVersion() const;

    /**
     * Determines the use of double buffering for all canvases.
     * May be overridden by command line option.
     * Default: enabled
     *
     * @note Can only be set before application initialization.
     */
    void setUseDoubleBuffering(bool useDoubleBuffering);

    /// Returns whether or not double buffering is used.
    bool useDoubleBuffering() const;

    //
    // Further settings
    //

    /**
     * Specifies whether the compiled-in modules are automatically loaded from
     * the module registration headers during initialization. Otherwise, desired modules need
     * to be instantiated "manually" and passed to registerModule() \e before initialization.
     * Default: true
     *
     * @note Can only be set before application initialization.
     */
    void setModuleLoadingEnabled(bool enabled);

    /// Returns whether modules are loaded automatically during initialization.
    bool isModuleLoadingEnabled() const;

    /**
     * In deployment mode the user-data directory is located in the home instead of voreen/data.
     * Default: false
     *
     * @note Can only be set before application initialization.
     */
    void setDeploymentMode(bool dm);

    /// Returns whether the application is in deployment mode.
    bool getDeploymentMode() const;

    //
    // Initialization
    //

    /// Returns the command line parser used by this application.
    CommandLineParser* getCommandLineParser() const;

    /// Returns the command queue used by this application
    CommandQueue* getCommandQueue() const;

    /**
     * Performs basic initializations as controlled by appType_, which do not require OpenGL access:
     *  - initialize tgt
     *  - execute command parser
     *  - start logging
     *  - detect paths
     *  - instantiates and initializes the module classes
     *
     * @note This function should be called right after object construction,
     *  and must especially be called before initializeGL().
     *
     * @throw VoreenException if initialization failed
     */
    virtual void initialize();

    /**
     * Deinitializes the application and deinitializes and deletes the module objects,
     * to be called right before object destruction.
     *
     * @note If OpenGL has been initialized, deinitializeGL() must be called first.
     *
     * @throw VoreenException if deinitialization failed
     */
    virtual void deinitialize();

    /**
     * Does OpenGL-specific initializations and calls initializeGL() of registered modules.
     *
     * @note initialize() must be called first.
     *
     * @throws VoreenException if OpenGL initialization failed
     */
    virtual void initializeGL();

    /**
     * Does OpenGL-specific deinitializations and calls deinitializeGL() of registered modules.
     *
     * @throws VoreenException if OpenGL deinitialization failed.
     */
    virtual void deinitializeGL();

    /**
     * Returns whether the application has been successfully initialized
     * (and not yet deinitialized).
     *
     * @see initialize
     */
    bool isInitialized() const;

    /**
     * Returns whether OpenGL-specific initializations have been performed.
     *
     * @see initializeGL
     */
    bool isInitializedGL() const;

    /**
     * Registers a network evaluator.
     *
     * Is internally used for scheduled network processing.
     * @see scheduleNetworkProcessing
     */
    void registerNetworkEvaluator(NetworkEvaluator* evaluator);

    /**
     * Deregisters a network evaluator.
     *
     * Is internally used for scheduled network processing.
     * @see scheduleNetworkProcessing
     */
    void deregisterNetworkEvaluator(NetworkEvaluator* evaluator);

    /**
     * Returns the first network evaluator in the set.
     */
    NetworkEvaluator* getNetworkEvaluator() const;

    /**
     * Returns the network evaluator belonging to the passed processor,
     * or null if no evaluator has a network assigned that contains the processor.
     */
    NetworkEvaluator* getNetworkEvaluator(Processor* p) const;

    /**
     * Returns the network evaluator that has the passed ProcessorNetwork assigned,
     * or null if no such evaluator exists.
     */
    NetworkEvaluator* getNetworkEvaluator(ProcessorNetwork* network) const;

    //
    // Settings (saved as properties)
    //
    bool useCaching() const;
    void setUseCaching(bool useCaching);

    /// Returns the CPU RAM limit in bytes set for the entire application, which is used for volume and octree memory management.
    size_t getCpuRamLimit() const;

    /// Returns the GPU memory limit in bytes which is used for managing the GPU texture memory available for volumes in memory management
    size_t getGpuMemoryLimit() const;

    //
    // Modules
    //

    /**
     * Registers a module.
     */
    void registerModule(VoreenModule* module);

    /**
     * Returns all registered modules.
     */
    const std::vector<VoreenModule*>& getModules() const;

    /**
     * Returns the VoreenModule with the passed name, or 0 if no such module exists.
     */
    VoreenModule* getModule(const std::string& moduleName) const;

    /**
     * Returns the absolute directory of the module with the passed name,
     * or an empty string, if no such module exists.
     */
    std::string getModulePath(const std::string& moduleName) const;


    //
    // Factory methods
    //

    /**
     * Implementation of SerializableFactory::createType:
     * Returns a new instance of the serializable class corresponding to the given @c typeString,
     * if such has been registered at one of the application's modules.
     *
     * Internally, this function iterates over all registered modules and
     * delegates the call to them.
     *
     * @see SerializableFactory::createType
     *
     * @returns either the new instance or @c 0, if no class with this type string
     *  has been registered.
     */
    virtual VoreenSerializableObject* createSerializableType(const std::string& typeString) const;

    /**
     * Returns the registered instance of the serializable class corresponding to the given @c typeString,
     * if such has been registered at one of the application's modules.
     *
     * Internally, this function iterates over all registered modules and
     * delegates the call to them.
     *
     * @returns either the registered instance or @c 0, if no class with this type string
     *  has been registered.
     */
    const VoreenSerializableObject* getSerializableType(const std::string& typeString) const;

    /**
     * Iterates over all registered modules and collects the registered
     * VoreenSerializableObject instances of the specified type.
     */
    template<typename T>
    std::vector<const T*> getSerializableTypes() const;

    /**
     * Implementation of SerializableFactory::getTypeString:
     * Returns a type string corresponding to the given @c type_info object,
     * in case this type has been registered at any of the application's modules.
     *
     * Internally, this function iterates over all registered modules and
     * delegates the call to them.
     *
     * @see SerializableFactory::getTypeString
     *
     * @returns either the string corresponding to the given @c type_info object,
     *          or an empty string, if the type has not been registered at any of the application's modules.
     */
    virtual std::string getSerializableTypeString(const std::type_info& type) const;

    /**
     * Registers a serializer factory.
     *
     * @deprecated Register types directly at your module as VoreenSerializableObject!
     */
    void registerSerializerFactory(SerializableFactory* factory);

    /**
     * Returns the application's serializer factories.
     *
     * @deprecated Register types directly at your module as VoreenSerializableObject!
     */
    const std::vector<SerializableFactory*>& getSerializerFactories() const;


    /**
     * Creates and returns a ProcessorWidget for the passed processor by delegation to
     * the ProcessorWidgetFactories of the registered modules. If no factory is able
     * to create a widget for the passed processor, the null pointer is returned.
     *
     * @note The creation of ProcessorWidgets can be disabled via the ApplicationFeatures.
     *  In this case, the method always returns null.
     */
    virtual ProcessorWidget* createProcessorWidget(Processor* processor) const;

    /**
     * Creates and returns a PropertyWidget for the passed property by delegation to
     * the PropertyWidgetFactories of the registered modules. If no factory is able
     * to create a widget for the passed property, the null pointer is returned.
     *
     * @note The creation of PropertyWidgets can be disabled via the ApplicationFeatures.
     *  In this case, the method always returns null.
     */
    virtual PropertyWidget* createPropertyWidget(Property* property) const;

    /**
     * Factory method for timers.
     *
     * @param handler The event handler that will be used
     *  for broadcasting the timer events. Must not be null.
     *
     * @note You have to override this function in a toolkit-specific subclass
     *  in order to actual create a timer. The standard implementation returns
     *  the null pointer.
     */
    virtual tgt::Timer* createTimer(tgt::EventHandler* handler) const;

    /**
     * Factory method for progress dialogs.
     *
     * @note You have to override this function in a toolkit-specific subclass
     *  in order to actual create a progress dialog. The standard implementation returns
     *  the null pointer.
     */
    virtual ProgressBar* createProgressDialog() const;

    /**
     * Displays a message box.
     *
     * @param error If true, an error message is shown, otherwise a standard info message.
     *
     * @note You have to override this function in a toolkit-specific subclass
     *  in order to actually display a message box. The standard implementation
     *  does nothing.
     */
    virtual void showMessageBox(const std::string& title, const std::string& message, bool error=false) const;


    //
    // Paths
    //

    /**
     * Returns the application's base path as detected by \sa init().
     */
    std::string getBasePath(const std::string& filename = "") const;

    /**
     * Returns the path, where the program binary is executed.
     */
    std::string getProgramPath() const;

    /**
    * Constructs an unique absolute file path within the specified directory
    *
    * @param root directory the unique path should be generated for
    * @param suffix is guaranteed to be the suffix of the generated file name
    *
    * Note that currently the implementation tries its best to avoid racing conditions,
    * and is guaranteed to generate a unique path for all calls within the same voreen
    * process, but (because it just returns a path and does not create the file) cannot
    * guarantee that (e.g.) the returned file path does not point to a file that has since
    * been created by another process.
    */
    std::string getUniqueFilePath(const std::string& root, const std::string& suffix = "") const;

    /**
     * Constructs an absolute path consisting of the cache directory
     * (getUserDataPath("cache")) and the given filename.
     */
    std::string getCachePath(const std::string& filename = "") const;

    /// Delete cache entries until cache is smaller than maxSize MBs.
    void cleanCache() const;
    /// Delete content of cache.
    void deleteCache();
    /// Clears the cache path property, thereby resetting the cache path to userDataPath/cache.
    void resetCachePath();

    /**
     * Constructs an absolute path consisting of the resource directory (voreen/resource/voreencore) and
     * the given filename.
     */
    std::string getCoreResourcePath(const std::string& filename = "") const;

    /**
     * Constructs an absolute path consisting of the resource directory (voreen/resource/appli) and
     * the given filename.
     */
    virtual std::string getApplicationResourcePath(const std::string& filename = "") const;

    /**
     * Constructs an absolute path consisting of the user-data directory (voreen/data or ~/Voreen) and
     * the given filename.
     */
    std::string getUserDataPath(const std::string& filename = "") const;

    /**
     * Constructs an absolute path consisting of the font directory (typically
     * voreen/data/fonts) and the given filename.
     */
    std::string getFontPath(const std::string& filename = "") const;

    /**
     * Constructs an absolute path consisting of the voreen temporary directory
     * and the given filename.
     */
    std::string getTemporaryPath(const std::string& filename = "") const;

    /**
     * Constructs an unique absolute file path within the voreen temporary directory
     *
     * @param suffix is guaranteed to be the suffix of the generated file name
     * @see getUniqueFilePath
     */
    std::string getUniqueTmpFilePath(const std::string& suffix = "") const;

    /**
     * Removed all temporary data created by the currently running voreen application.
     */
    void cleanTemporaryData();

    /**
     * Returns the test data directory. If not set, an empty string is returned.
     *
     * @see setTestDataPath
     */
    std::string getTestDataPath() const;

    /**
     * Sets the test data directory. It is stored in
     * the application settings using a FileDialogProperty.
     */
    void setTestDataPath(const std::string& path);

#ifdef __APPLE__
    /**
     * Constructs an absolute path consisting of the Mac application bundle's resource
     * directory (path Contents/Resources within the bundle) and the given filename.
     */
    std::string getAppBundleResourcesPath(const std::string& filename = "") const;
#endif

    /**
     * Add initial application properties
     */
    virtual void initApplicationSettings();
    void loadApplicationSettings();
    void saveApplicationSettings();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:
    /**
     * Overwrite this method to load additional modules.
     */
    virtual void loadModules();

    /**
     * This function triggers a non-blocking network processing,
     * i.e., the assigned NetworkEvaluator's process() function is called
     * after the current call stack has been completed (event queue concept).
     *
     * This function is called by each end processor when it receives an invalidation.
     *
     * @note Since internally a tgt::Timer is used for scheduling, this
     *  function does only work in a derived class that overrides createTimer().
     *  Otherwise, it is a no-op.
     */
    virtual void scheduleNetworkProcessing();

    /**
     * Queries the currently available graphics memory via GpuCapabilities
     * and updates the respective property.
     */
    void queryAvailableGraphicsMemory();

    /**
     * When the user activates / deactivates the automatic setting, the cpu ram limit has to be set read-only and / or set to an automatic limit
     */
    virtual void determineRamLimitAutomaticallyChanged();

    /// notifies the memory manager to update the available memory
    virtual void cpuRamLimitChanged();

    /**
     * When the user activates / deactivates the automatic setting, the gpu memory limit has to be set read-only and / or set to an automatic limit
     */
    virtual void determineGpuMemoryLimitAutomaticallyChanged();

    /// notifies the memory manager to update the available gpu texture memory
    virtual void gpuMemoryLimitChanged();

    /// Prints an ascii-art version of the voreen logo to std::cout
    virtual void printAsciiLogo() const;

    /// Sets the application's CPU RAM limit in bytes.
    void setCpuRamLimit(size_t ramLimit);

    static const std::string loggerCat_;

private:
    /**
     * Helper function for scheduled network processing.
     * @see scheduleNetworkProcessing
     */
    virtual void timerEvent(tgt::TimeEvent* e);

    /// Sets up the logging framework according to the current property settings.
    void initLogging(std::string& htmlLogFile);

    void logLevelChanged();
    void tempDataPathChanged();

    static VoreenApplication* app_;

    ApplicationFeatures appFeatures_;  ///< application configuration as passed to the constructor
    std::string binaryName_;    ///< name of the application binary set in the constructor

    // command line options
    std::string overrideGLSLVersion_;   ///< if set, the detected GLSL version is overridden by this

    // further settings
    bool loadModules_;          ///< if true, modules are auto-loaded from the module registration headers
    bool deploymentMode_;       ///< in deployment mode the user-data directory is located in the home instead of voreen/data.

    // paths detected during initialization
    std::string basePath_;
    std::string userDataPath_;

    std::vector<VoreenModule*> modules_;
    std::vector<SerializableFactory*> serializerFactories_;

    std::unique_ptr<CommandLineParser> cmdParser_;
    std::unique_ptr<CommandQueue> cmdQueue_;

    std::set<NetworkEvaluator*> networkEvaluators_;
    tgt::Timer* schedulingTimer_;       ///< Timer for scheduled network processing.
    tgt::EventHandler eventHandler_;    ///< Local event handler for the scheduling events.
    bool networkEvaluationRequired_;    ///< set to true by scheduleNetworkProcessing and checked on timer event

    // settings properties for caching
    BoolProperty useCaching_;
    IntProperty volumeCacheLimit_;
    IntProperty octreeCacheLimit_;
    FileDialogProperty cachePath_;
    ButtonProperty resetCachePath_;
    ButtonProperty deleteCache_;

    // CPU RAM properties
    BoolProperty determineCpuRamLimitAutomatically_;
    IntProperty cpuRamLimit_;

    // GPU memory properties
    BoolProperty determineGpuMemoryLimitAutomatically_;
    IntProperty gpuMemoryLimit_;

    // setting properties regarding GPU memory
    IntProperty availableGraphicsMemory_;
    ButtonProperty refreshAvailableGraphicsMemory_;

    // setting for double buffering
    bool useDoubleBuffering_;

    // logging settings
    BoolProperty enableLogging_;               ///< if false, console and HTML file logging is disabled
    OptionProperty<tgt::LogLevel> logLevel_;   ///< log level to use for console and HTML loggers
    BoolProperty enableHTMLLogging_;           ///< determines whether logging to an HTML file is enabled
    FileDialogProperty htmlLogFile_;           ///< HTML log file

    // settings properties for regression testing
    FileDialogProperty testDataPath_;

    // settings properties for temporary data
    FileDialogProperty tempDataPath_;

    // not used directly by the VoreenApplication, but may be quried by an actual application (e.g. VoreenVE)
    BoolProperty showSplashScreen_;             //< if true, a splash screen should be shown on start up
    BoolProperty askForSave_;                   //< if true, the user should be asked, if he wants to save the modified workspace
    //BoolProperty loadLastWorkspaceOnStartup_;   //< if true, the last loaded workspace should be restored
    BoolProperty showStartupWizard_; //< if true, shows the startup wizard, else creates an empty workspace

    // Generator for uuids (mutable, since ()-operator is non-const)
    mutable boost::uuids::random_generator uuidGenerator_;

    // The tempDataPath for this particular instance of Voreen
    std::string tempDataPathInstance_;

    bool initialized_;
    bool initializedGL_;
};

//-----------------------------------------------------------------------------
// template definitions

template<typename T>
std::vector<const T*> voreen::VoreenApplication::getSerializableTypes() const {
    std::vector<const T*> result;
    for (size_t i=0; i<modules_.size(); i++) {
        std::vector<const T*> modResources = modules_[i]->getSerializableTypes<T>();
        result.insert(result.end(), modResources.begin(), modResources.end());
    }
    return result;
}

} // namespace

#endif
