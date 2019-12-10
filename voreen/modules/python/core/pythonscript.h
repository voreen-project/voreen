/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_PYTHONSCRIPT_H
#define VRN_PYTHONSCRIPT_H

#include <string>

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/serialization/serializable.h"

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

namespace voreen {

/**
 * Stores a Python script, runs it and checks for errors.
 * Additionally, the script source can be compiled to byte code,
 * resulting in a faster script execution and clearer error messages.
 *
 * @note The caller has to make sure that the Python interpreter is initialized
 *  during all interactions with this class.
 *
 * @note Convenient loading of Python scripts with search path handling
 *  is provided by the PythonModule class.
 */
class VRN_CORE_API PythonScript : public Serializable {

public:
    PythonScript();

    /**
     * Frees the stored byte code. Make sure that the
     * Python interpreter is still initialized
     * when deleting the script.
     */
    ~PythonScript();

    PythonScript(const PythonScript& other);
    PythonScript& operator=(const PythonScript& other);

    /**
     * Returns the uuid of the script. This is mainly used for output forwarding.
     */
    const std::string& getId() const;

    /**
     * Loads a Python script from file and compiles it to byte code by default.
     *
     * @return true, if the script has been loaded successfully.
     *  Note: Does not check the script for validity.
     */
    bool load(const std::string& filename, bool compile = true);

    /**
     * Assigns the Python script source, replacing the current source.
     */
    void setSource(const std::string& source);

    /**
     * Returns the script's source.
     */
    const std::string& getSource() const;

    /**
     * Assigns a filename to the script. This is optionally,
     * since it is only used in debug/error messages.
     */
    void setFilename(const std::string& filename);

    /**
     * Returns the script's filename. Usually, this is
     * the name of file the script has been loaded from,
     * but it may have been overwritten or even empty.
     */
    const std::string& getFilename() const;

    /**
     * Compiles the script source to byte code, which speeds up script execution
     * and is generally recommended, since it also produces more clear error messages.
     *
     * If an error occurs, the error message is stored and can be retrieved through getLog().
     *
     * @param logErrors if true, error messages are also passed to the logger.
     *  The internal log buffer is not affected by this parameter.
     *
     * @return true, if script compilation has been successful
     */
    bool compile(bool logErrors = true);

    /**
     * Runs the script source as string or runs the byte code,
     * if the script has been compiled.
     *
     * If an error occurs, the error message is stored and can be retrieved through getLog().
     *
     * @param logErrors if true, error messages are also passed to the logger.
     *  The internal log buffer is not affected by this parameter.
     *
     * @return true, if script execution has been successful
     */
    bool run(bool logErrors = true);

    /**
     * Returns the error that has occurred during the last operation (compilation or execution).
     * If the last operation has been successful, an empty string is returned.
     */
    const std::string& getLog() const;

    /**
     * Returns the source line number where the last error has occurred.
     * Returns -1, if the last operation has been successful or the line number could not be retrieved.
     */
    int getErrorLine() const;

    /**
     * Returns the source column number where the last error has occurred.
     * Returns -1, if the last operation has been successful or the col number could not be retrieved.
     */
    int getErrorCol() const;

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Two scripts are considered equal, if both their source and filename are equal.
     */
    bool operator==(const PythonScript& other) const;

    /**
     * Two scripts are not considered equal, if any of their source or filename are not equal.
     */
    bool operator!=(const PythonScript& other) const;

private:
    bool checkCompileError(bool logErrors = true);
    bool checkRuntimeError(bool logErrors = true);

    std::string     id_;
    std::string     source_;
    std::string     filename_;
    PyObject*       byteCode_;
    bool            compiled_;
    std::string     log_;
    int             errorLine_;
    int             errorCol_;

    bool printErrorsToStdOut_;  //< only for debugging

    static const std::string loggerCat_;
};

}

#endif // VRN_PYTHONSCRIPT_H
