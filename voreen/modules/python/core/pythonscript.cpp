/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

// include this at very first
#include <Python.h>

#include "../pythonmodule.h"
#include "pythonscript.h"

#include "tgt/filesystem.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <traceback.h>
#include <frameobject.h>

namespace voreen {

const std::string PythonScript::loggerCat_ = "voreen.Python.PythonScript";

PythonScript::PythonScript()
    : source_("")
    , filename_("")
    , byteCode_(nullptr)
    , compiled_(false)
    , errorLine_(-1)
    , errorCol_(-1)
    , printErrorsToStdOut_(false)
{
    if(VoreenApplication::app()) {
        id_ = boost::lexical_cast<std::string>(VoreenApplication::app()->generateUUID());
    }
}

PythonScript::~PythonScript() {
    Py_XDECREF(byteCode_);
}

PythonScript::PythonScript(const PythonScript& other) {
    id_ = other.id_;
    source_ = other.source_;
    filename_ = other.filename_;
    byteCode_ = other.byteCode_;
    Py_XINCREF(byteCode_);
    compiled_ = other.compiled_;
    errorLine_ = other.errorLine_;
    errorCol_ = other.errorCol_;
    printErrorsToStdOut_ = other.printErrorsToStdOut_;
}

PythonScript& PythonScript::operator=(const PythonScript& other) {
    id_ = other.id_;
    source_ = other.source_;
    filename_ = other.filename_;
    Py_XDECREF(byteCode_);
    byteCode_ = other.byteCode_;
    Py_XINCREF(byteCode_);
    compiled_ = other.compiled_;
    errorLine_ = other.errorLine_;
    errorCol_ = other.errorCol_;
    printErrorsToStdOut_ = other.printErrorsToStdOut_;

    return *this;
}

const std::string& PythonScript::getId() const {
    return id_;
}

const std::string& PythonScript::getLog() const {
    return log_;
}

bool PythonScript::load(const std::string& filename, bool pCompile) {
    tgt::File* file = FileSys.open(filename);

    // Check if file is open
    if (!file || !file->isOpen())
        return false;

    size_t len = file->size();

    // Check if file is empty
    if (len == 0)
        return false;

    // allocate memory
    char* sourceChar = nullptr;
    try {
        sourceChar = new char[len+1];
    }
    catch (std::bad_alloc&) {
        return false;
    }

    file->read(sourceChar, len);
    sourceChar[len] = 0; // terminate with 0
    source_ = std::string(sourceChar);
    delete[] sourceChar;

    file->close();
    delete file;
    filename_ = filename;

    // convert windows line breaks to unix, since some Python versions
    // seem to have problems with them
    source_ = strReplaceAll(source_, "\r\n", "\n");

    if (pCompile)
        compile();

    return true;
}

bool PythonScript::compile(bool logErrors) {
    LDEBUG("Compiling script '" << getFilename() << "'");

    Py_XDECREF(byteCode_);
    byteCode_ = Py_CompileString(source_.c_str(), tgt::FileSystem::fileName(filename_).c_str(), Py_file_input);
    compiled_ = checkCompileError(logErrors);

    if (!compiled_) {
        Py_CLEAR(byteCode_);
    }

    return compiled_;
}

bool PythonScript::run(bool logErrors) {
    // use clean environment containing only the built-in symbols,
    // in order to prevent influences of previous script executions
    PyObject* glb = PythonModule::getInstance()->cleanGlobals();

    // Add the script id for correct output redirection.
    PyObject* id = PyUnicode_FromString(id_.c_str());
    PyDict_SetItemString(glb, VOREEN_SCRIPT_ID_IDENTIFIER, id);

    bool success;
    if (compiled_){
        tgtAssert(byteCode_, "No byte code");
        LDEBUG("Running compiled script '" << getFilename() << "' ...");

        // eval compiled byte code and check for errors
        PyObject* dum = PyEval_EvalCode(byteCode_, glb, glb);
        success = checkRuntimeError(logErrors);
        Py_XDECREF(dum);
    }
    else {
        LDEBUG("Running script '" << getFilename() << "' ...");
        PyObject* dum = PyRun_String(source_.c_str(), Py_file_input, glb, glb);
        success = checkRuntimeError(logErrors);
        Py_XDECREF(dum);
    }

    Py_CLEAR(id);

    if (success) {
        LDEBUG("finished.");
    }

    return success;
}

const std::string& PythonScript::getSource() const {
    return source_;
}

void PythonScript::setSource(const std::string& source) {
    source_ = source;
    compiled_ = false;
    Py_CLEAR(byteCode_);
}

void PythonScript::setFilename(const std::string& filename) {
    filename_ = filename;
}

const std::string& PythonScript::getFilename() const {
    return filename_;
}

int PythonScript::getErrorLine() const {
    return errorLine_;
}

int PythonScript::getErrorCol() const {
    return errorCol_;
}

void PythonScript::serialize(voreen::Serializer& s) const {
    s.serialize("source", source_);
    s.serialize("filename", filename_);
}

void PythonScript::deserialize(voreen::Deserializer& s) {
    s.deserialize("source", source_);
    s.deserialize("filename", filename_);
}

bool PythonScript::operator==(const voreen::PythonScript& other) const {
    return (filename_ == other.filename_ && source_ == other.source_);
}

bool PythonScript::operator!=(const voreen::PythonScript& other) const {
    return !(*this == other);
}

bool PythonScript::checkCompileError(bool logErrors) {

    log_ = "";
    errorLine_ = -1;
    errorCol_ = -1;
    if (!PyErr_Occurred())
        return true;

    PyObject *errtype, *errvalue, *traceback;
    PyErr_Fetch(&errtype, &errvalue, &traceback);

    char* msg = 0;
    PyObject* obj = 0;
    if (PyArg_ParseTuple(errvalue, "sO", &msg, &obj)) {
        int line, col;
        char *code = 0;
        char *mod = 0;
        if (PyArg_ParseTuple(obj, "siis", &mod, &line, &col, &code)) {
            errorLine_ = line;
            errorCol_ = col;
            log_ = "[" + itos(line) + ":" + itos(col) + "] " + std::string(msg) + ": " + std::string(code);
            if (!getFilename().empty()) {
                log_ += "  File \"" + tgt::FileSystem::fileName(getFilename()) + "\"";
            }
        }
    }

    // convert error to string, if it could not be parsed
    if (log_.empty()) {
        if (logErrors)
            LWARNING("Failed to parse exception, printing as string:");
        PyObject* s = PyObject_Str(errvalue);
        if (s && PyUnicode_Check(s)) {
            log_ = PyUnicodeAsString(s);
            Py_XDECREF(s);
        }
    }

    if (printErrorsToStdOut_) {
        PyErr_Restore(errtype, errvalue, traceback);
        PyErr_Print();
    }
    else {
        Py_XDECREF(errtype);
        Py_XDECREF(errvalue);
        Py_XDECREF(traceback);
    }

    if (logErrors)
        LERROR(log_);

    return false;
}

bool PythonScript::checkRuntimeError(bool logErrors) {
    using std::string;

    log_ = "";
    errorLine_ = -1;
    errorCol_ = -1;
    if (!PyErr_Occurred())
        return true;

    std::string pyException = "";
    PyObject* pyError_type = 0;
    PyObject* pyError_value = 0;
    PyObject* pyError_traceback = 0;
    PyObject* pyError_string = 0;
    PyErr_Fetch(&pyError_type, &pyError_value, &pyError_traceback);

    // 1. Read and append error type
    /*if (pyError_type != 0 && (pyError_string = PyObject_Str(pyError_type)) != 0 && (PyUnicode_Check(pyError_string))) {
        pyException.append(PyUnicodeAsString(pyError_string));
        Py_XDECREF(pyError_string);
        pyError_string = 0;
    }
    else {
        pyException.append("<Exception unknown>");
    }
    pyException.append("\n"); */

    // 2. Read traceback (reverting order of elements in order to list the most recent one first)
    std::string stacktraceStr;
    if (pyError_traceback) {
        PyTracebackObject* traceback = (PyTracebackObject*)pyError_traceback;
        while (traceback) {
            PyFrameObject* frame = traceback->tb_frame;
            std::string stacktraceLine;
            if (frame && frame->f_code) {
                PyCodeObject* codeObject = frame->f_code;
                if (PyUnicode_Check(codeObject->co_filename))
                    stacktraceLine.append(string("  File \"") + PyUnicodeAsString(codeObject->co_filename) + string("\", "));

                errorLine_ = PyCode_Addr2Line(codeObject, frame->f_lasti);
                stacktraceLine.append(string("line ") + itos(errorLine_));

                if (PyUnicode_Check(codeObject->co_name))
                    stacktraceLine.append(string(", in ") + PyUnicodeAsString(codeObject->co_name));
            }
            stacktraceLine.append("\n");
            stacktraceStr = stacktraceLine + stacktraceStr;

            traceback = traceback->tb_next;
        }
    }

    // 3. Read and append error string
    std::stringstream s;
    s << errorLine_;
    pyException.append(string("[") + s.str() + string("] "));
    if (pyError_value && (pyError_string = PyObject_Str(pyError_value)) != 0 && (PyUnicode_Check(pyError_string))) {
        pyException.append(PyUnicodeAsString(pyError_string));
        Py_XDECREF(pyError_string);
        pyError_string = 0;
    }
    else {
        pyException.append("<No data available>");
    }
    pyException.append("\n");

    // finally append stacktrace string
    if (!stacktraceStr.empty()) {
        pyException.append("Stacktrace (most recent call first):\n");
        pyException.append(stacktraceStr);
    }
    else {
        pyException.append("<No stacktrace available>");
        if (logErrors)
            LWARNING("Failed to parse traceback");
    }

    if (printErrorsToStdOut_) {
        PyErr_Restore(pyError_type, pyError_value, pyError_traceback);
        PyErr_Print();
    }
    else {
        Py_XDECREF(pyError_type);
        Py_XDECREF(pyError_value);
        Py_XDECREF(pyError_traceback);
    }

    log_ = pyException;

    if (logErrors)
        LERROR(log_);

    return false;
}

} // namespace voreen
