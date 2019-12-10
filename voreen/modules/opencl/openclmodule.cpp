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

#include "openclmodule.h"
#include "utils/clwrapper.h"

#include "processors/grayscale_cl.h"
#include "processors/raycaster_cl.h"
#include "processors/dynamicclprocessor.h"
#include "processors/raytracingentryexitpoints.h"
#include "processors/singleoctreeraycastercl.h"
#include "processors/volumegradient_cl.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

OpenCLModule* OpenCLModule::instance_ = 0;
const std::string OpenCLModule::loggerCat_ = "voreen.OpenCLModule";

OpenCLModule::OpenCLModule(const std::string& modulePath)
    : VoreenModule(modulePath)
    , opencl_(nullptr)
    , context_(nullptr)
    , queue_(nullptr)
    , glSharing_("sharingProp", "GL Sharing", true)
    , deviceProp_("deviceProp", "Device:", Processor::INVALID_RESULT, true)
    , currentDeviceIdx_(-1)
{
    setID("OpenCL");
    setGuiName("OpenCL");
    instance_ = this;

    registerSerializableType(new DynamicCLProcessor());
    registerSerializableType(new GrayscaleCL());
    registerSerializableType(new RaycasterCL());
    registerSerializableType(new RaytracingEntryExitPoints());
    registerSerializableType(new SingleOctreeRaycasterCL());
    registerSerializableType(new VolumeGradientCL());

    addShaderPath(getModulePath("glsl"));
    addProperty(deviceProp_);
    addProperty(glSharing_);
}

void OpenCLModule::initializeGL() {

    if (tgt::Singleton<tgt::GpuCapabilities>::isInited()) {
        VoreenModule::initializeGL();
    }
    else if (glSharing_.get()) {
        LWARNING("No OpenGL context available: disabling OpenGL sharing");
        glSharing_.set(false);
    }

    // Initialize OpenCL in order to set up platform property.
#ifdef WIN32
    _putenv_s("CUDA_CACHE_DISABLE", "1");
#else
    setenv("CUDA_CACHE_DISABLE", "1", 1);
#endif

    opencl_ = new cl::OpenCL();

    // Fill plattform.
    std::deque<Option<int>> options;
    const std::vector<cl::Platform>& platforms = opencl_->getPlatforms();
    for (size_t i = 0; i < platforms.size(); i++) {
        const std::vector<cl::Device>& devices = platforms.at(i).getDevices();
        for (size_t j = 0; j < devices.size(); j++) {
            // Key format is: <platform>:<device>(devicenum).
            std::string key = platforms.at(i).getName() + ":" + devices.at(j).getName() + "(" + std::to_string(j) + ")";
            options.push_back(Option<int>( key, key, static_cast<int>(devices_.size()) ));
            devices_.push_back(std::make_pair(static_cast<int>(i), static_cast<int>(j)));
        }
    }

    // There might be no valid device.
    if (devices_.empty()) {
        deviceProp_.setReadOnlyFlag(true);
        throw VoreenException("Found no OpenCL device");
    }

    bool devicesChanged = deviceProp_.getOptions().size() != options.size();
    for (size_t i = 0; i < deviceProp_.getOptions().size() && !devicesChanged; i++) {
        devicesChanged = deviceProp_.getOptions().at(i).key_ != options.at(i).key_;
    }

    deviceProp_.setOptions(options);

    if (devicesChanged) {

        // Notify user.
        //VoreenApplication::app()->showMessageBox("OpenCL Device changed", "OpenCL devices changed since the last start. An automatic selection will be performed.");
        LWARNING("OpenCL devices changed since the last start. An automatic selection will be performed.");
        deviceProp_.selectByIndex(0);

        // Prefer GPU.
        if (tgt::Singleton<tgt::GpuCapabilities>::isInited()) {
            const cl::Platform platform = opencl_->getPlatformByVendor(GpuCaps.getVendor());
            const std::vector<cl::Device> devices = platform.getDevices();
            for (size_t i = 0; i < devices.size(); i++) {
                if (devices.at(i).getType() == CL_DEVICE_TYPE_GPU) {
                    std::string key = platform.getName() + ":" + devices.at(i).getName() + "(" + std::to_string(i) + ")";
                    if (deviceProp_.hasKey(key))
                        deviceProp_.selectByKey(key);
                    else {
                        LERROR("Preferred device on platform not available. Taking the first available.");
                        break;
                    }
                }
            }
        }
    }

    glSharing_.onChange(MemberFunctionCallback<OpenCLModule>(this, &OpenCLModule::setupDevice));
    deviceProp_.onChange(MemberFunctionCallback<OpenCLModule>(this, &OpenCLModule::setupDevice));
    deviceProp_.updateWidgets();
    setupDevice();
}

void OpenCLModule::deinitializeGL() {
    devices_.clear();

    delete queue_;
    queue_ = 0;
    delete context_;
    context_ = 0;
    delete opencl_;
    opencl_ = 0;

    if (tgt::Singleton<tgt::GpuCapabilities>::isInited()) {
        VoreenModule::deinitializeGL();
    }
}

void OpenCLModule::setupDevice() {

    // Clear old resources.
    if (context_) {
        std::vector<cl::OpenCLProcessorBase*> reinitialized;
        for (cl::OpenCLProcessorBase* processor : processors_) {
            processor->requestDeinitializeCL();
            if (!processor->isDeviceChangeSupported())
                reinitialized.push_back(processor);
        }

        delete context_;
        context_ = nullptr;
        delete queue_;
        queue_ = nullptr;

        if (!reinitialized.empty())
            VoreenApplication::app()->showMessageBox("OpenCL Device changed", "The workspace needs to be reloaded in order to apply the selected device.");
    }

    // Apply selection.
    currentDeviceIdx_ = deviceProp_.getValue();
    platform_ = opencl_->getPlatforms().at(devices_.at(currentDeviceIdx_).first);
    device_ = platform_.getDevices().at(devices_.at(currentDeviceIdx_).second);

    // Log infos of selected device.
    device_.logInfos();

    // Create context.
    if (glSharing_.get()) {
        LINFO("OpenGL sharing: enabled");
        context_ = new cl::Context(cl::Context::generateGlSharingProperties(), device_);
    }
    else {
        LINFO("OpenGL sharing: disabled");
        context_ = new cl::Context(device_);
    }

    // Create command queue.
#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    cl_command_queue_properties properties = CL_QUEUE_PROFILING_ENABLE;
#else
    const cl_queue_properties properties[] = { CL_QUEUE_PROPERTIES, CL_QUEUE_PROFILING_ENABLE, 0};
#endif
    queue_ = new cl::CommandQueue(context_, device_, properties);

    // Reinitialize processors.
    for (cl::OpenCLProcessorBase* processor : processors_) {
        if (processor->isDeviceChangeSupported())
            processor->requestInitializeCL();
    }
}

cl::OpenCL* OpenCLModule::getOpenCL() const {
    if (!opencl_)
        LERROR("No OpenCL wrapper. Call initCL first!");
    return opencl_;
}

cl::Context* OpenCLModule::getCLContext() const {
    if (!context_)
        LERROR("No OpenCL context. Call initCL first!");
     return context_;
}

cl::CommandQueue* OpenCLModule::getCLCommandQueue() const {
    if (!queue_)
        LERROR("No OpenCL queue. Call initCL first!");
    return queue_;
}

cl::Platform& OpenCLModule::getCLPlatform() {
    if (!opencl_)
        LERROR("No OpenCL platform. Call initCL first!");
    return platform_;
}

cl::Device& OpenCLModule::getCLDevice() {
    if (device_.getId() == 0)
        LERROR("No OpenCL device. Call initCL first!");
    return device_;
}

void OpenCLModule::setGLSharing(bool enabled) {
    glSharing_.set(enabled);
}

bool OpenCLModule::getGLSharing() const {
    return glSharing_.get();
}

OpenCLModule* OpenCLModule::getInstance() {
    return instance_;
}

cl::Program* OpenCLModule::loadProgram(const std::string& path) {
    if (!opencl_)
        throw VoreenException("OpenCLModule: no OpenCL wrapper. Call initCL first!");

    LINFO("Loading program " << path);
    cl::Program* prog = new cl::Program(context_);
    if (!prog->loadSource(path)) {
        delete prog;
        throw VoreenException("Failed to load OpenCL program: " + path);
    }

    if (!prog->build(device_)) {
        delete prog;
        throw VoreenException("Failed to build OpenCL program: " + path);
    }

    return prog;
}

void OpenCLModule::registerProcessor(cl::OpenCLProcessorBase* processor) {
    processors_.insert(processor);
}
void OpenCLModule::deregisterProcessor(cl::OpenCLProcessorBase* processor) {
    processors_.erase(processor);
}

} // namespace
