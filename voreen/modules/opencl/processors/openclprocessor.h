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

#ifndef VRN_OPENCLPROCESSOR_H
#define VRN_OPENCLPROCESSOR_H

#include "voreen/core/voreencoreapi.h"
#include "modules/opencl/openclmodule.h"

namespace voreen {
namespace cl {

class VRN_CORE_API OpenCLProcessorBase {
public:
    OpenCLProcessorBase();
    virtual ~OpenCLProcessorBase();

    //! Returns if OpenCL data is currently initialized.
    virtual bool isInitializedCL() const;

    //! Requests the deriving processor to initialize CL related data.
    virtual void requestInitializeCL();

    //! Requests the deriving processor to deinitialize CL related data.
    virtual void requestDeinitializeCL();

    //! Determines, whether changing the device on the fly is supported by this processor. Default is false.
    virtual bool isDeviceChangeSupported() const;

protected:

    //! Initializes the actual OpenCL data used by deriving processors.
    virtual void initializeCL() = 0;

    //! Deinitializes the actual OpenCL data used by deriving processors.
    virtual void deinitializeCL() = 0;

    bool initializedCL_; ///< Determines, whether the opencl data is currently initialized.
};

template<typename T>
class VRN_CORE_API OpenCLProcessor : public T, public OpenCLProcessorBase {
    static_assert(std::is_base_of<Processor, T>::value, "T must derive from Processor");
public:

    OpenCLProcessor() : T(), OpenCLProcessorBase() {}
    virtual ~OpenCLProcessor() {}

    virtual bool isReady() const {
        if (!isInitializedCL() && !isDeviceChangeSupported()) {
            T::setNotReadyErrorMessage("OpenCL device changed");
            return false;
        }
        return T::isReady();
    }

protected:

    virtual void initialize() {
        if (!OpenCLModule::getInstance()->getCLContext())
            throw VoreenException("No OpenCL context created");

        OpenCLModule::getInstance()->registerProcessor(this);
        T::initialize();
        requestInitializeCL();
    }
    virtual void deinitialize() {
        requestDeinitializeCL();
        T::deinitialize();
        OpenCLModule::getInstance()->deregisterProcessor(this);
    }
};

}
}

#endif // VRN_OPENCLPROCESSOR_H
