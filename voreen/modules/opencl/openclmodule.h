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

#ifndef VRN_OPENCLMODULE_H
#define VRN_OPENCLMODULE_H

#include "voreen/core/voreenmodule.h"

#include "modules/opencl/utils/clwrapper.h"

#include "voreen/core/properties/stringproperty.h"

#include <set>

namespace voreen {

namespace cl {
    class OpenCLProcessorBase;
    //class OpenCL;
    //class Context;
    //class CommandQueue;
    //class Program;
    //class Device;
}

class VRN_CORE_API OpenCLModule : public VoreenModule {

public:
    OpenCLModule(const std::string& modulePath);

    virtual std::string getDescription() const {
        return "Processors utilizing OpenCL.";
    }

    /**
    * Allocates OpenCL resources.
    */
    virtual void initializeGL();

    /**
     * Frees the allocated OpenCL resources.
     */
    virtual void deinitializeGL();

    /**
     * Returns the OpenCL wrapper.
     */
    cl::OpenCL* getOpenCL() const;

    /**
     * Returns the OpenCL context.
     */
    cl::Context* getCLContext() const;

    /**
     * Returns the OpenCL command queue.
     */
    cl::CommandQueue* getCLCommandQueue() const;

    /**
     * Returns the OpenCL platform.
     */
    cl::Platform& getCLPlatform();

    /**
     * Returns the OpenCL device.
     */
    cl::Device& getCLDevice();

    /**
     * \brief Enables or disables sharing of resources between
     *  OpenGL and OpenCL. Default: enabled
     */
    void setGLSharing(bool enabled);

    /**
     * Returns whether OpenGL sharing is enabled.
     */
    bool getGLSharing() const;

    /**
     * Loads an OpenCL kernel from file.
     *
     * @throw VoreenException, if the kernel could not be loaded.
     */
    cl::Program* loadProgram(const std::string& path);

    /**
     * Returns the global instance of this class.
     *
     * @note Does not create the instance. If the module class has not been
     *       instantiated yet, the null pointer is returned.
     */
    static OpenCLModule* getInstance();

    /// Register and unregister a processor in order to invalidate if device changed.
    void registerProcessor(cl::OpenCLProcessorBase* processor);
    void deregisterProcessor(cl::OpenCLProcessorBase* processor);

private:

    /// Sets up the choosen device in case it has been changed.
    void setupDevice(); ///< Callen by deviceProp_

    // OpenCL resources
    cl::OpenCL* opencl_;
    cl::Context* context_;
    cl::CommandQueue* queue_;
    cl::Platform platform_;
    cl::Device device_;

    std::vector<std::pair<int, int>> devices_;
    OptionProperty<int> deviceProp_;
    int currentDeviceIdx_;
    
    std::set<cl::OpenCLProcessorBase*> processors_;

    bool glSharing_;    ///< Determines whether OpenGL sharing is enabled

    static OpenCLModule* instance_;

    static const std::string loggerCat_;
};

} // namespace

#endif
