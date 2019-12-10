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

#ifndef VRN_DEVILMODULE_H
#define VRN_DEVILMODULE_H

#include "voreen/core/voreenmodule.h"

namespace voreen {

class VRN_CORE_API DevILModule : public VoreenModule {

public:
    DevILModule(const std::string& modulePath);

    virtual std::string getDescription() const {
        return "Provides image handling functionality via the DevIL library.";
    }

    /**
     * Returns the last DevIL error as string
     * and pops it off the error stack.
     */
    static std::string getDevILError();

    /**
     * Returns the image file extensions that are supported for reading.
     */
    static std::vector<std::string> getSupportedReadExtensions();

    /**
     * Returns the image file extensions that are supported for writing.
     */
    static std::vector<std::string> getSupportedWriteExtensions();

protected:
    virtual void initialize();

    static std::vector<std::string> readExtensions_;  //< supported image file extensions (reading)
    static std::vector<std::string> writeExtensions_; //< supported image file extensions (writing)

};

} // namespace

#endif // VRN_DEVILMODULE_H
