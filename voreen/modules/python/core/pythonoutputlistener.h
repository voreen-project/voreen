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

#ifndef VRN_PYTHONOUTPUTLISTENER_H
#define VRN_PYTHONOUTPUTLISTENER_H

#include "voreen/core/voreencoreapi.h"

/**
 * Implement this interface in order to receive the output of
 * Python scripts.
 *
 * @see PythonModule::addOutputListener
 */
class VRN_CORE_API PythonOutputListener {
public:
    /// Receives sys.stdout
    virtual void pyStdout(const std::string& out, const std::string& id) = 0;

    /// Receives sys.stderr
    virtual void pyStderr(const std::string& err, const std::string& id) = 0;
};

#endif // VRN_PYTHONOUTPUTLISTENER_H

