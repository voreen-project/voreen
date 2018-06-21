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

// include this at very first
#include <Python.h>
#include "modules/python/pythonmodule.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen {

char* PyUnicodeAsString(PyObject* object) {

    tgtAssert(object != nullptr, "Object was null");
    tgtAssert(PyUnicode_Check(object), "No Unicode Object");

    char* result = NULL;
    PyObject* bytes = PyUnicode_AsEncodedString(object, "UTF-8", "strict");
    if (bytes != NULL) {
        result = PyBytes_AS_STRING(bytes);
        result = strdup(result);
        Py_DECREF(bytes);
    }

    return result;
}

const std::string PythonModule::loggerCat_("voreen.Python.Module");
PythonModule* PythonModule::instance_ = 0;
std::vector<PythonOutputListener*> PythonModule::outputListeners_;

// -----------------------------------------------------------------
// Python functions

static PyObject* voreen_print(PyObject* /*self*/, PyObject* args) {
    char* msg;
    int len;
    int isStderr;
    if (!PyArg_ParseTuple(args, "s#i", &msg, &len, &isStderr)) {
        LWARNINGC("voreen.Python.voreen_print", "failed to parse log message");
    }
    else {
        if (len > 1 || ((len == 1) && (msg[0] != '\0') && (msg[0] != '\r') && (msg[0] != '\n'))) {
            std::string message(msg);
            const std::vector<PythonOutputListener*> listeners = PythonModule::getOutputListeners();
            if (!listeners.empty()) {
                // pass output to listeners
                for (size_t i=0; i<listeners.size(); i++) {
                    if (isStderr)
                        listeners[i]->pyStderr(message);
                    else
                        listeners[i]->pyStdout(message);
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
{
    setID("Python");
    setGuiName("Python");
    instance_ = this;
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
    addModulePath(getModulePath("ext/python27/modules"));
#endif

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
}

void PythonModule::deinitialize() {

    if (!isInitialized())
        throw VoreenException("not initialized");

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
