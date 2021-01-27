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

#ifndef VRN_PYTHONMODULE_H
#define VRN_PYTHONMODULE_H

#include "core/pythonscript.h"
#include "core/pythonoutputlistener.h"
#include "core/pyvoreen.h"

#include "voreen/core/voreenmodule.h"
#include "tgt/manager.h"

#include <string>

namespace voreen {

// Helper to convert Unicode PyObject into c-string
std::string PyUnicodeAsString(PyObject* object);

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API tgt::ResourceManager<PythonScript>;
#endif

extern const char* VOREEN_SCRIPT_ID_IDENTIFIER;

class VRN_CORE_API PythonModule : public VoreenModule, public tgt::ResourceManager<PythonScript> {

public:
    PythonModule(const std::string& modulePath);

    ~PythonModule();

    virtual std::string getDescription() const {
        return "Provides a Python interface to Voreen resources, such as the processor network.";
    }

    /**
     * Loads a Python script from the specified file, considering the current search path.
     *
     * @see ResourceManager
     *
     * @return the created script object, or null if loading failed
     */
    PythonScript* loadScript(const std::string& filename, bool compile = true);

    /**
     * Loads a Python script from the specified path and executes it,
     * considering the current search path,
     *
     * @see loadScript
     *
     * @param logErrors if true, error messages are passed to the logger.
     *
     * @throw VoreenException if script execution failed, contains the error message
     */
    void runScript(const std::string& filename, bool logErrors = true);

    /**
     * Modifies sys.argv in the Python environment.
     *
     * @param argc number of command line arguments
     * @param argv array of arguments
     */
    void setArgv(int argc, wchar_t * argv[]);

    /**
     * Appends the passed path to the Python interpreter's module search path.
     */
    void addModulePath(const std::string& path);

    /**
     * Registers a listener to receive the output of Python scripts.
     */
    void addOutputListener(PythonOutputListener* listener);

    /**
     * Removes an output listener.
     */
    void removeOutputListener(PythonOutputListener* listener);

    /**
     * Checks the state of the Python error indicator.
     * If an error has occurred, false is returned and the
     * type and value variables are written to the passed string reference.
     */
    bool checkForPythonError(std::string& errorMsg);

    /**
     * Use this environment in calls to PyEval_EvalCode or PyRun_String. The caller does not
     * assume ownership of the returned object! Do not decrement reference count!
     *
     * You are free to add items to the dictionary via PyDict_Set* calls.
     */
    PyObject* cleanGlobals();

    /**
     * Returns the global instance of this class.
     *
     * @note Does not create the instance. If the module class has not been
     *       instantiated yet, the null pointer is returned.
     */
    static PythonModule* getInstance();

    /// To accessed from the vrn_print() function.
    static const std::vector<PythonOutputListener*>& getOutputListeners();

protected:
    /**
     * Initializes the Python interpreter and the Voreen bindings.
     */
    virtual void initialize();

    /**
     * Disposes the Python interpreter and the Voreen bindings.
     */
    virtual void deinitialize();

private:
    /// The actual Voreen Python bindings (uses RAII)
    PyVoreen pyVoreen_;

    // Dictionary to use for calls to PyEval_EvalCode
    PyObject* globals_;
    PyObject* initGlobals_;

    /// Output of Python scripts is redirected to these.
    static std::vector<PythonOutputListener*> outputListeners_;

    static PythonModule* instance_;

    static const std::string loggerCat_;
};

} // namespace

#endif
