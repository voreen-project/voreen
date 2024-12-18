/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

// See https://docs.python.org/3/c-api/arg.html#strings-and-buffers
#define PY_SSIZE_T_CLEAN
// include this at very first
#include <Python.h>
#include "modules/python/pythonmodule.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"

#include "processors/dynamicpythonprocessor.h"

namespace voreen {

std::string PyUnicodeAsString(PyObject* object) {

    tgtAssert(object != nullptr, "Object was null");
    tgtAssert(PyUnicode_Check(object), "No Unicode Object");

    std::string result;
    PyObject* bytes = PyUnicode_AsEncodedString(object, "UTF-8", "strict");
    if (bytes != NULL) {
        result = PyBytes_AS_STRING(bytes);
        Py_DECREF(bytes);
    }

    return result;
}

const char* VOREEN_SCRIPT_ID_IDENTIFIER = "__voreen_script_id__";

const std::string PythonModule::loggerCat_("voreen.Python.Module");
PythonModule* PythonModule::instance_ = 0;
std::vector<PythonOutputListener*> PythonModule::outputListeners_;

// -----------------------------------------------------------------
// Python functions

static PyObject* voreen_print(PyObject* /*self*/, PyObject* args) {
    char* msg;
    Py_ssize_t len;
    int isStderr;
    if (!PyArg_ParseTuple(args, "s#i", &msg, &len, &isStderr)) {
        LWARNINGC("voreen.Python.voreen_print", "failed to parse log message");
    } else {
        if (len > 1 || ((len == 1) && (msg[0] != '\0') && (msg[0] != '\r') && (msg[0] != '\n'))) {
            PyObject* scriptIdObj = PyDict_GetItemString(PyEval_GetGlobals(), VOREEN_SCRIPT_ID_IDENTIFIER);
            std::string id = PyUnicodeAsString(scriptIdObj);
            std::string message(msg);
            const std::vector<PythonOutputListener*> listeners = PythonModule::getOutputListeners();
            if (!listeners.empty()) {
                // pass output to listeners
                for (size_t i=0; i<listeners.size(); i++) {
                    if (isStderr)
                        listeners[i]->pyStderr(message, id);
                    else
                        listeners[i]->pyStdout(message, id);
                }
            }
            else {
                // log output if no listener registered
                if (isStderr)
                    LWARNINGC("voreen.Python.voreen_print", message);
                else
                    LINFOC("voreen.Python.voreen_print", message);
            }
        }
    }

    Py_RETURN_NONE;
}

// module 'voreen_internal'
static PyMethodDef internal_methods[] = {
    {
        "vrnPrint",
        voreen_print,
        METH_VARARGS,
        "Internal helper function used for accessing script output"
    },
    { NULL, NULL, 0, NULL} // sentinal
};

static struct PyModuleDef voreenInternalModuleDef =
{
    PyModuleDef_HEAD_INIT,
    "voreen_internal",
    NULL,
    -1,
    internal_methods
};

PyMODINIT_FUNC
PyInit_voreenInternalModule(void)
{
    return PyModule_Create(&voreenInternalModuleDef);
}

// -----------------------------------------------------------------
PythonModule::PythonModule(const std::string& modulePath)
    : VoreenModule(modulePath)
    , tgt::ResourceManager<PythonScript>(false)
    , globals_(nullptr)
    , initGlobals_(nullptr)
{
    setID("Python");
    setGuiName("Python");
    instance_ = this;

    registerProcessor(new DynamicPythonProcessor());
}

PythonModule::~PythonModule() {
    instance_ = nullptr;
}

PythonModule* PythonModule::getInstance() {
    return instance_;
}

void PythonModule::initialize() {
    VoreenModule::initialize();

    //
    // Register module extensions (needs to be done before(!) Py_Initialize()
    //
    if(PyImport_AppendInittab("voreen_internal", PyInit_voreenInternalModule) == -1)
        LWARNING("Failed to init helper module 'voreen_internal'");

    //
    // Initialize Python interpreter
    //
    LINFO("Python version: " << Py_GetVersion());

    // Pass program name to the Python interpreter
    static wchar_t str_pyvoreen[] = L"PyVoreen";
    Py_SetProgramName(str_pyvoreen);

#ifdef WIN32
#ifndef VRN_USE_PYTHON_VERSION
#error VRN_USE_PYTHON_VERSION needs to be defined by module cmake file
#endif

    // Set Python home to embedded build
    std::string libraryPath = getModulePath("ext/" VRN_USE_PYTHON_VERSION);
    static std::wstring str_pyhome; // we need static storage duration for PythonHome.
    str_pyhome = std::wstring(libraryPath.begin(), libraryPath.end());
    Py_SetPythonHome(str_pyhome.c_str());

    // disable import of 'site' module (not available in embedded windows build)
    Py_NoSiteFlag = 1;
#endif

    // Initialize the Python interpreter. Required.
    Py_InitializeEx(0);
    if (!Py_IsInitialized())
        throw VoreenException("Failed to initialize Python interpreter");

    // required in order to use threads.
    PyEval_InitThreads();

    // init ResourceManager search path
    addPath("");
    addPath(VoreenApplication::app()->getCoreResourcePath("scripts"));
    addPath(getModulePath("scripts"));

    // init Python's internal module search path
    addModulePath(VoreenApplication::app()->getCoreResourcePath("scripts"));
    addModulePath(getModulePath("scripts"));
#ifdef WIN32
    addModulePath(libraryPath + "/lib");
#endif

    // Define an initial (empty) global environment
    globals_ = PyDict_New();
    initGlobals_ = PyDict_New();

    //
    // Redirect script output from std::cout to voreen_print function (see above)
    //

    // load output redirector script and run it once
    std::string filename = "outputcatcher.py";
    LDEBUG("Loading Python init script '" << filename << "'");
    PythonScript* initScript = loadScript(filename);
    if (initScript) {
        if (!initScript->run())
            LWARNING("Failed to run init script '" << filename << "': " << initScript->getLog());
        dispose(initScript);
    }
    else {
        LWARNING("Failed to load init script '" << filename << "'");
    }
    Py_CLEAR(initGlobals_);
    initGlobals_ = PyDict_Copy(globals_);
    PyDict_DelItemString(initGlobals_, VOREEN_SCRIPT_ID_IDENTIFIER);
}

void PythonModule::deinitialize() {

    if (!isInitialized())
        throw VoreenException("not initialized");

    Py_CLEAR(globals_);
    Py_CLEAR(initGlobals_);

    // clean up python interpreter
    Py_Finalize();

    VoreenModule::deinitialize();
}

PythonScript* PythonModule::loadScript(const std::string& filename, bool compileDirectly) {

    // do not check isInitialized(), since we call this function from initialize()
    if (!Py_IsInitialized()) {
        LWARNING("load(): not initialized");
        return 0;
    }

    if (isLoaded(filename)) {
        increaseUsage(filename);
        return get(filename);
    }

    PythonScript* script = new PythonScript();
    if (script->load(completePath(filename), compileDirectly)) {
        reg(script, filename);
        return script;
    }
    delete script;

    return 0;
}

void PythonModule::runScript(const std::string& filename, bool logErrors) {

    if (!isInitialized())
        throw VoreenException("PythonModule not initialized");

    PythonScript* script = loadScript(filename, false);
    if (!script)
        throw VoreenException("Failed to load Python script '" + filename + "'");

    std::string errorMsg;
    if (script->compile()) {
        LDEBUG("Running Python script '" << filename << "' ...");
        if (script->run(false))
            LDEBUG("Python script finished.");
        else {
            errorMsg = "Python runtime error:\n" + script->getLog();
            if (logErrors)
                LERROR(errorMsg);
        }
    }
    else {
        errorMsg = "Python compile error:\n" + script->getLog();
        if (logErrors)
            LERROR(errorMsg);
    }
    dispose(script);

    if (!errorMsg.empty())
        throw VoreenException(errorMsg);
}

void PythonModule::setArgv(int argc, wchar_t* argv[]) {
    PySys_SetArgv(argc, argv);
}

void PythonModule::addModulePath(const std::string& path) {

    // do not check isInitialized(), since we call this function from initialize()
    if (!Py_IsInitialized()) {
        LWARNING("addModulePath(): not initialized");
        return;
    }

    // convert windows back slashes to slashes
    std::string pathConv = strReplaceAll(path, "\\", "/");

    LDEBUG("Adding '" << pathConv << "' to Python module search path");
    std::string runString = "import sys\n";
    runString.append(std::string("sys.path.append('") + pathConv + std::string("')"));
    int ret = PyRun_SimpleString(runString.c_str());
    if (ret != 0)
        LWARNING("Failed to add '" << pathConv << "' to Python module search path");
}

bool PythonModule::checkForPythonError(std::string& errorMsg) {

    if (!PyErr_Occurred())
        return true;

    // retrieve error indicator into variables
    PyObject *errtype, *errvalue, *traceback;
    PyErr_Fetch(&errtype, &errvalue, &traceback);
    PyErr_NormalizeException(&errtype, &errvalue, &traceback);

    // convert indicators to Python strings, if possible
    PyObject* errtypeStr = PyObject_Str(errtype);
    PyObject* errvalueStr = PyObject_Str(errvalue);
    if (errtypeStr && PyUnicode_Check(errtypeStr))
        errorMsg = PyUnicodeAsString(errtypeStr);
    if (errvalueStr && PyUnicode_Check(errvalueStr))
        errorMsg += std::string(", ") + PyUnicodeAsString(errvalueStr);

    // free allocated ressources
    Py_XDECREF(errtypeStr);
    Py_XDECREF(errvalueStr);

    Py_XDECREF(errtype);
    Py_XDECREF(errvalue);
    Py_XDECREF(traceback);

    return false;
}

PyObject* PythonModule::cleanGlobals() {
    // Here's the deal: For some reason, the FIRST call to PyEval_EvalCode or PyRun_String "registers" the
    // passed globals dictionary as THE globals dictionary for ever and ever and ever. Subsequent calls will
    // not use the (potentially other) dictionary passed in that call and instead refer to the first passed
    // dictionary. Why? I don't know.
    //
    // However, we want to be able to supply the script with a "clean", but pre setup environment. initGlobals_
    // IS this clean argument with some setup done (i.e., redirections for stdout/stderr, voreen_internal module).
    // We therefore reset the globals_ dictionary (the one that will be passed always (and thus also the first time)
    // to PyEval_EvalCode and PyRun_String) to have the same values as initGlobals_.
    PyDict_Clear(globals_);
    PyDict_Merge(globals_, initGlobals_, true);

    int ret = PyDict_SetItemString(globals_, "__builtins__", PyEval_GetBuiltins());
    tgtAssert(ret == 0, "Setting builtins failed!");

    return globals_;
}

void PythonModule::addOutputListener(PythonOutputListener* listener) {
    tgtAssert(listener, "null pointer passed");
    if (std::find(outputListeners_.begin(), outputListeners_.end(), listener) == outputListeners_.end())
        outputListeners_.push_back(listener);
    else
        LWARNING("Output listener already registered");
}

void PythonModule::removeOutputListener(PythonOutputListener* listener) {
    tgtAssert(listener, "null pointer passed");
    std::vector<PythonOutputListener*>::iterator it = std::find(outputListeners_.begin(), outputListeners_.end(), listener);
    if (it != outputListeners_.end())
        outputListeners_.erase(it);
    else
        LWARNING("Output listener not registered");
}

const std::vector<PythonOutputListener*>& PythonModule::getOutputListeners() {
    return outputListeners_;
}

} // namespace
