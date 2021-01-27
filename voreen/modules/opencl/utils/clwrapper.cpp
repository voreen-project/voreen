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

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "clwrapper.h"
#include "modules/opencl/openclmodule.h"
#include "tgt/filesystem.h"
#include "tgt/memory.h"
#include "voreen/core/utils/stringutils.h"

#include <typeinfo>

#ifdef WIN32
#include <windows.h>
#endif

#if !(defined(WIN32) || defined(__APPLE__))
#include <GL/glx.h>
#define CL_GLX_DISPLAY_KHR 0x200A
#endif

//TODO: remove defines?
#define CL_GL_CONTEXT_KHR 0x2008
#define CL_EGL_DISPLAY_KHR 0x2009
#define CL_WGL_HDC_KHR 0x200B
#define CL_CGL_SHAREGROUP_KHR 0x200C

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenGL/OpenGL.h>
#include <OpenCL/cl_gl_ext.h>
#endif

namespace voreen {

namespace cl {

void CL_API_CALL contextCallback(const char *errinfo, const void* /*private_info*/, size_t /*cb*/, void* /*user_data*/) {
    LINFOC("contextCallback", std::string(errinfo));
}

const std::string OpenCL::loggerCat_ = "voreen.OpenCL";

OpenCL::OpenCL() {
    LINFO("OpenCL init");

    cl_uint numPlatforms;
    LCL_ERROR(clGetPlatformIDs(0, 0, &numPlatforms));

    if (numPlatforms == 0)
        return;

    cl_platform_id* platforms = new cl_platform_id[numPlatforms];
    LCL_ERROR(clGetPlatformIDs(numPlatforms, platforms, 0));

    for(cl_uint i=0; i<numPlatforms; ++i) {
        platforms_.push_back(Platform(platforms[i]));
        LINFO("Platform " << i << ": " << platforms_.back().getName());
    }
    delete[] platforms;
}

cl::Platform OpenCL::getPlatformByVendor(tgt::GpuCapabilities::GpuVendor vendor) const {
    if (platforms_.empty()) {
        LWARNING("No OpenCL platforms found");
        return Platform(0);
    }
    else if (platforms_.size() == 1) {
        return platforms_.front();
    }
    else { // more than one platform available
        switch(vendor) {
        case tgt::GpuCapabilities::GPU_VENDOR_NVIDIA:
            for (const Platform& platform : platforms_) {
                if (toLower(platform.getVendor()).find("nvidia") != std::string::npos)
                    return platform;
            }
            break;
        case tgt::GpuCapabilities::GPU_VENDOR_ATI:
            for (const Platform& platform : platforms_) {
                if ((toLower(platform.getVendor()).find("amd") != std::string::npos) ||
                    (toLower(platform.getVendor()).find("ati") != std::string::npos))
                    return platform;
            }
            break;
        case tgt::GpuCapabilities::GPU_VENDOR_INTEL:
            for (const Platform& platform : platforms_) {
                if (toLower(platform.getVendor()).find("intel") != std::string::npos)
                    return platform;
            }
            break;
        case tgt::GpuCapabilities::GPU_VENDOR_UNKNOWN:
        default:
            break;
        }

        // no matching platform found
        return platforms_.front();
    }
}

Device OpenCL::getCurrentDeviceForGlContext() {
    //cl_device_id dev;

    //std::vector<ContextProperty> properties = Context::generateGlSharingProperties();
    //int numProps = (properties.size() * 2)+1;
    //cl_context_properties* props = new cl_context_properties[numProps];
    //for(size_t i=0; i<properties.size(); ++i) {
        //props[i*2] = properties[i].name_;
        //props[(i*2)+1] = properties[i].value_;
    //}
    //props[numProps-1] = 0;     //terminate list
    ////clGetGLContextInfoKHR_fn fn = (cl_int (*CL_API_CALL) (const cl_context_properties*, cl_gl_context_info , size_t , void*, size_t*)) clGetExtensionFunctionAddress("clGetGLContextInfoKHR");
    ////clGetGLContextInfoKHR = fn*;
    ////clGetGLContextInfoKHR = (cl_int (*CL_API_CALL) (const cl_context_properties*, cl_gl_context_info , size_t , void*, size_t*)) clGetExtensionFunctionAddress("clGetGLContextInfoKHR");
    ////if(fn == 0)
        ////LERROR("Could not find function pointer!");
    ////LCL_ERROR(clGetGLContextInfoKHR(props, CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR, sizeof(dev), &dev, 0));
    ////LCL_ERROR((*fn)(props, CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR, sizeof(dev), &dev, 0));
    //delete props;
    //if(dev != 0)
        //return new Device(dev);
    //else
        return 0;
}

std::vector<Device> OpenCL::getDevicesForGlContext() {
    std::vector<Device> devices;

    /*std::vector<ContextProperty> properties = Context::generateGlSharingProperties();
    int numProps = (properties.size() * 2)+1;
    cl_context_properties* props = new cl_context_properties[numProps];
    for(size_t i=0; i<properties.size(); ++i) {
        props[i*2] = properties[i].name_;
        props[(i*2)+1] = properties[i].value_;
    }
    props[numProps-1] = 0;     //terminate list

    typedef CL_API_ENTRY cl_int (CL_API_CALL *clGetGLContextInfoKHR_fn)(
        const cl_context_properties*,
        cl_gl_context_info,
        size_t,
        void*,
        size_t*);
    clGetGLContextInfoKHR_fn fn = (clGetGLContextInfoKHR_fn)clGetExtensionFunctionAddress("clGetGLContextInfoKHR");

    size_t retSize;
    LCL_ERROR(fn(props, CL_DEVICES_FOR_GL_CONTEXT_KHR, 0, 0, &retSize));
    int numDevs = retSize/sizeof(cl_device_id);
    cl_device_id* devs = new cl_device_id[numDevs];
    LCL_ERROR(fn(props, CL_DEVICES_FOR_GL_CONTEXT_KHR, retSize, devs, 0));
    delete props;
    for (int i=0; i<numDevs; ++i) {
        if(devs[i])
            devices.push_back(Device(devs[i]));
    }
    delete devs; */

    return devices;
}

size_t OpenCL::shrRoundUp(size_t localWorkSize, size_t numItems) {
    size_t result = localWorkSize;
    while (result < numItems)
        result += localWorkSize;
    tgtAssert((result >= numItems) && ((result % localWorkSize) == 0), "invalid post-condition");
    return result;
}

//---------------------------------------------------------

const std::string Platform::loggerCat_ = "voreen.OpenCL.Platform";

Platform::Platform(cl_platform_id id) : id_(id), profile_(UNKNOWN) {
    std::string profileString = getInfo<std::string>(CL_PLATFORM_PROFILE); //FULL_PROFILE or EMBEDDED_PROFILE

    if (profileString == "FULL_PROFILE")
        profile_ = FULL_PROFILE;
    else if(profileString == "EMBEDDED_PROFILE") {
        profile_ = EMBEDDED_PROFILE;
        LWARNING("Embedded profile!");
    }
    else {
        LERROR("Unknown profile!");
    }

    name_ = getInfo<std::string>(CL_PLATFORM_NAME);
    vendor_ = getInfo<std::string>(CL_PLATFORM_VENDOR);
    extensionString_ = getInfo<std::string>(CL_PLATFORM_EXTENSIONS);
    versionString_ = getInfo<std::string>(CL_PLATFORM_VERSION);     //OpenCL<space><major_version.minor_version><space><platform-specific information>
    version_.parseVersionString(versionString_);

    //explode extensions string with space as delimiter:
    std::string str = extensionString_;
    size_t found;
    found = str.find_first_of(" ");
    while(found != std::string::npos){
        if (found > 0){
            extensions_.insert(str.substr(0,found));
        }
        str = str.substr(found+1);
        found = str.find_first_of(" ");
    }
    if (str.length() > 0){
        extensions_.insert(str);
    }

    //LINFO("Extensions:");
    //for( std::set<std::string>::const_iterator iter = extensions_.begin(); iter != extensions_.end(); ++iter ) {
        //LINFO(">" << *iter << "<");
    //}

    //find devices:
    cl_uint numDevices;
    cl_int error_code = clGetDeviceIDs(id_, CL_DEVICE_TYPE_ALL, 0, 0, &numDevices);
    switch(error_code) {
        case CL_SUCCESS:
            break;
        case CL_DEVICE_NOT_FOUND:
            // => 0 devices, see standard
            return;
        default:
            LCL_ERROR(error_code);
    }

    if (numDevices == 0)
        return;

    cl_device_id* devices = new cl_device_id[numDevices];
    LCL_ERROR(clGetDeviceIDs(id_, CL_DEVICE_TYPE_ALL, numDevices, devices, 0));

    for(cl_uint i=0; i<numDevices; ++i) {
        devices_.push_back(Device(devices[i]));
    }

    delete[] devices;
}

Platform::Platform() {
    id_ = 0;
}

Platform::~Platform() {
    //while(!devices_.empty()) {
        //delete devices_.back();
        //devices_.pop_back();
    //}
}

bool Platform::isExtensionSupported(std::string ext) const {
    if(extensions_.find(ext) != extensions_.end())
        return true;
    else
        return false;
}

template<>
std::string Platform::getInfo(cl_platform_info info) const {
    size_t retSize;
    LCL_ERROR(clGetPlatformInfo(id_, info, 0, 0, &retSize));
    char* buffer = new char[retSize];
    LCL_ERROR(clGetPlatformInfo(id_, info, retSize, buffer, 0));
    std::string ret(buffer);
    delete[] buffer;
    return ret;
}

void Platform::logInfos() const {
    LINFO("Name: " << name_);
    LINFO("Profile: " << getInfo<std::string>(CL_PLATFORM_PROFILE););
    LINFO("Version: " << versionString_);
    LINFO("Parsed Version: " << version_);
    LINFO("Vendor: " << vendor_);
    LINFO("Extensions: " << extensionString_);

#if defined (__APPLE__) || defined(MACOSX)
    static const std::string clGlSharingExtString = "cl_APPLE_gl_sharing";
#else
    static const std::string clGlSharingExtString = "cl_khr_gl_sharing";
#endif

    if (isExtensionSupported(clGlSharingExtString))
        LINFO("GL sharing is supported.");
    else
        LWARNING("GL sharing is not supported.");

    LINFO("Number of devices: " << devices_.size());
}

//-----------------------------------------------------------------------------------

const Platform::ClVersion Platform::ClVersion::VRN_CL_VERSION_1_0(1,0);

Platform::ClVersion::ClVersion(int major, int minor)
  : major_(major), minor_(minor)
{}

bool Platform::ClVersion::parseVersionString(const std::string& st) {
    major_ = -1;
    minor_ = -1;

    std::string str;
    //remove "OpenCL ":
    str = st.substr(7);

    //ignore vendor specific part of the string:
    size_t spacePos = str.find_first_of(" ");
    if (spacePos != std::string::npos) {
        str = str.substr(0, spacePos);
    }

    //explode version string with delimiter ".":
    std::vector<std::string> exploded;
    size_t found;
    found = str.find_first_of(".");
    while(found != std::string::npos){
        if (found > 0){
            exploded.push_back(str.substr(0,found));
        }
        str = str.substr(found+1);
        found = str.find_first_of(".");
    }
    if (str.length() > 0){
        exploded.push_back(str);
    }

    // parse numbers
    if (exploded.size() < 2) {
        LWARNING("Version string too short to parse!");
        return false;
    }

    std::stringstream vstr;

    vstr << exploded[0];
    vstr >> major_;
    if (vstr.fail()) {
        LERROR("Failed to parse major version! Version string: " << st);
        major_ = -1;
        return false;
    }

    vstr.clear();
    vstr.str("");

    vstr << exploded[1];
    vstr >> minor_;
    if (vstr.fail()) {
        LERROR("Failed to parse minor version! Version string: " << st);
        major_ = -1;
        minor_ = -1;
        return false;
    }

    vstr.clear();
    vstr.str("");

    return true;
}


bool operator==(const Platform::ClVersion& x, const Platform::ClVersion& y) {
    if ((x.major_ == y.major_) && (x.minor_ == y.minor_))
        return true;
    else
        return false;
}

bool operator!=(const Platform::ClVersion& x, const Platform::ClVersion& y) {
    if ((x.major_ != y.major_) || (x.minor_ != y.minor_))
        return true;
    else
        return false;
}

bool operator<(const Platform::ClVersion& x, const Platform::ClVersion& y) {
    if (x.major_ < y.major_)
        return true;
    else if (x.major_ == y.major_)
        if (x.minor_ < y.minor_)
            return true;

    return false;
}

bool operator<=(const Platform::ClVersion& x, const Platform::ClVersion& y) {
    if (x.major_ < y.major_)
        return true;
    else if (x.major_ == y.major_)
        if (x.minor_ < y.minor_)
            return true;

    return false;
}

bool operator>(const Platform::ClVersion& x, const Platform::ClVersion& y) {
    if (x.major_ > y.major_)
        return true;
    else if (x.major_ == y.major_)
        if (x.minor_ > y.minor_)
            return true;

    return false;
}

bool operator>=(const Platform::ClVersion& x, const Platform::ClVersion& y) {
    if (x.major_ > y.major_)
        return true;
    else if (x.major_ == y.major_)
        if (x.minor_ > y.minor_)
            return true;

    return false;
}

std::ostream& operator<<(std::ostream& s, const Platform::ClVersion& v) {
    if (v.major_ == -1)
        return (s << "unknown");
    else
        return (s << v.major_ << "." << v.minor_);
}

//---------------------------------------------------------

const std::string Device::loggerCat_ = "voreen.OpenCL.Device";

Device::Device() : id_(0) {
}

Device::Device(cl_device_id id) : id_(id) {
    name_ = getInfo<std::string>(CL_DEVICE_NAME);

    extensionString_ = getInfo<std::string>(CL_DEVICE_EXTENSIONS);
    //Returns a space separated list of extension names (the extension names themselves do not contain any spaces). The list of extension names returned currently can include one or more of the following approved extension names:
    //cl_khr_fp64
    //cl_khr_select_fprounding_mode
    //cl_khr_global_int32_base_atomics
    //cl_khr_global_int32_extended_atomics
    //cl_khr_local_int32_base_atomics
    //cl_khr_local_int32_extended_atomics
    //cl_khr_int64_base_atomics
    //cl_khr_int64_extended_atomics
    //cl_khr_3d_image_writes
    //cl_khr_byte_addressable_store
    //cl_khr_fp16
    std::vector<std::string> extensionVec = strSplit(extensionString_, ' ');
    extensions_.insert(extensionVec.begin(), extensionVec.end());

    imageSupport_ = getInfo<cl_bool>(CL_DEVICE_IMAGE_SUPPORT);
    if (imageSupport_) {
        maxImageSize2D_.x = static_cast<int>(getInfo<size_t>(CL_DEVICE_IMAGE2D_MAX_WIDTH));
        maxImageSize2D_.y = static_cast<int>(getInfo<size_t>(CL_DEVICE_IMAGE2D_MAX_HEIGHT));

        maxImageSize3D_.x = static_cast<int>(getInfo<size_t>(CL_DEVICE_IMAGE3D_MAX_WIDTH));
        maxImageSize3D_.y = static_cast<int>(getInfo<size_t>(CL_DEVICE_IMAGE3D_MAX_HEIGHT));
        maxImageSize3D_.z = static_cast<int>(getInfo<size_t>(CL_DEVICE_IMAGE3D_MAX_DEPTH));
    }
}

//template specialization for strings:
template<>
std::string Device::getInfo(cl_device_info info) const {
    size_t retSize;
    LCL_ERROR(clGetDeviceInfo(id_, info, 0, 0, &retSize));
    char* buffer = new char[retSize];
    LCL_ERROR(clGetDeviceInfo(id_, info, retSize, buffer, 0));
    std::string ret(buffer);
    delete[] buffer;
    return ret;
}

cl_uint Device::getMaxWorkGroupSize() const {
    cl_uint result = static_cast<cl_uint>(getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE));
    return result;
}

void Device::logInfos() const {
    //General info about the device:-------------------------------------------------------------------
    LINFO("Name: " << name_);
    //Device name string.

    LINFO("Vendor: " << getInfo<std::string>(CL_DEVICE_VENDOR    ));
    //Vendor name string.

    LINFO("Vendor ID: " << getInfo<cl_uint>(CL_DEVICE_VENDOR_ID));
    //A unique device vendor identifier. An example of a unique device identifier could be the PCIe ID.

    LINFO("CL_DEVICE_VERSION: " << getInfo<std::string>(CL_DEVICE_VERSION));
    //OpenCL version string. Returns the OpenCL version supported by the device. This version string has the following format:
    //OpenCL<space><major_version.minor_version><space><vendor-specific information>
    //The major_version.minor_version value returned will be 1.0.

    LINFO("CL_DRIVER_VERSION: " << getInfo<std::string>(CL_DRIVER_VERSION));
    //OpenCL software driver version string in the form major_number.minor_number.

    //LINFO("CL_DEVICE_PLATFORM: " << getInfo<cl_platform_id>(CL_DEVICE_PLATFORM));
    //The platform associated with this device.

    LINFO("CL_DEVICE_TYPE: " << getInfo<cl_device_type>(CL_DEVICE_TYPE));
    //The OpenCL device type. Currently supported values are one of or a combination of: CL_DEVICE_TYPE_CPU, CL_DEVICE_TYPE_GPU, CL_DEVICE_TYPE_ACCELERATOR, or CL_DEVICE_TYPE_DEFAULT.

    LINFO("CL_DEVICE_AVAILABLE: " << getInfo<cl_bool>(CL_DEVICE_AVAILABLE));
    //Is CL_TRUE if the device is available and CL_FALSE if the device is not available.

    LINFO("CL_DEVICE_COMPILER_AVAILABLE: " << getInfo<cl_bool>(CL_DEVICE_COMPILER_AVAILABLE));
    //Is CL_FALSE if the implementation does not have a compiler available to compile the program source. Is CL_TRUE if the compiler is available. This can be CL_FALSE for the embededed platform profile only.

    LINFO("Extensions: " << extensionString_);
    //Returns a space separated list of extension names (the extension names themselves do not contain any spaces). The list of extension names returned currently can include one or more of the following approved extension names:
    //cl_khr_fp64
    //cl_khr_select_fprounding_mode
    //cl_khr_global_int32_base_atomics
    //cl_khr_global_int32_extended_atomics
    //cl_khr_local_int32_base_atomics
    //cl_khr_local_int32_extended_atomics
    //cl_khr_int64_base_atomics
    //cl_khr_int64_extended_atomics
    //cl_khr_3d_image_writes
    //cl_khr_byte_addressable_store
    //cl_khr_fp16

    LINFO("CL_DEVICE_PROFILE: " << getInfo<std::string>(CL_DEVICE_PROFILE));
    //OpenCL profile string. Returns the profile name supported by the device (see note). The profile name returned can be one of the following strings:
    //FULL_PROFILE - if the device supports the OpenCL specification (functionality defined as part of the core specification and does not require any extensions to be supported).
    //EMBEDDED_PROFILE - if the device supports the OpenCL embedded profile.

    LINFO("CL_DEVICE_MAX_CLOCK_FREQUENCY: " << getInfo<cl_uint>(CL_DEVICE_MAX_CLOCK_FREQUENCY) << " MHz");
    //Maximum configured clock frequency of the device in MHz.

    LINFO("CL_DEVICE_PROFILING_TIMER_RESOLUTION: " << getInfo<size_t>(CL_DEVICE_PROFILING_TIMER_RESOLUTION) << " ns");
    //Describes the resolution of device timer. This is measured in nanoseconds.

    //Image support:----------------------------------------------------------------------------------

    if (imageSupport_) {
        LINFO("Image support: yes; Max sizes: 2D " << maxImageSize2D_ <<  ", 3D " << maxImageSize3D_);

        LINFO("CL_DEVICE_MAX_READ_IMAGE_ARGS: " << getInfo<cl_uint>(CL_DEVICE_MAX_READ_IMAGE_ARGS     ));
        //Max number of simultaneous image objects that can be read by a kernel. The minimum value is 128 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.

        LINFO("CL_DEVICE_MAX_WRITE_IMAGE_ARGS: " << getInfo<cl_uint>(CL_DEVICE_MAX_WRITE_IMAGE_ARGS     ));
        //Max number of simultaneous image objects that can be written to by a kernel. The minimum value is 8 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.

        LINFO("CL_DEVICE_MAX_SAMPLERS: " << getInfo<cl_uint>(CL_DEVICE_MAX_SAMPLERS     ));
        //Maximum number of samplers that can be used in a kernel. The minimum value is 16 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE. (Also see sampler_t.)
    }
    else {
        LINFO("Image support: no");
    }

    //-------------------------------------------------------------------------------------

    LINFO("CL_DEVICE_ADDRESS_BITS: " << getInfo<cl_uint>(CL_DEVICE_ADDRESS_BITS       ));
    //The default compute device address space size specified as an unsigned integer value in bits. Currently supported values are 32 or 64 bits.

    //LINFO("CL_DEVICE_DOUBLE_FP_CONFIG: " << getInfo<cl_device_fp_config>(CL_DEVICE_DOUBLE_FP_CONFIG     ));
    //Describes the OPTIONAL double precision floating-point capability of the OpenCL device. This is a bit-field that describes one or more of the following values:
        //* CL_FP_DENORM - denorms are supported.
        //* CL_FP_INF_NAN - INF and NaNs are supported.
        //* CL_FP_ROUND_TO_NEAREST - round to nearest even rounding mode supported.
        //* CL_FP_ROUND_TO_ZERO - round to zero rounding mode supported.
        //* CL_FP_ROUND_TO_INF - round to +ve and -ve infinity rounding modes supported.
        //* CP_FP_FMA - IEEE754-2008 fused multiply-add is supported.
    //The mandated minimum double precision floating-point capability is CL_FP_FMA | CL_FP_ROUND_TO_NEAREST | CL_FP_ROUND_TO_ZERO | CL_FP_ROUND_TO_INF | CL_FP_INF_NAN | CL_FP_DENORM.

    LINFO("CL_DEVICE_ENDIAN_LITTLE: " << getInfo<cl_bool>(CL_DEVICE_ENDIAN_LITTLE     ));
    //Is CL_TRUE if the OpenCL device is a little endian device and CL_FALSE otherwise.

    LINFO("CL_DEVICE_ERROR_CORRECTION_SUPPORT: " << getInfo<cl_bool>(CL_DEVICE_ERROR_CORRECTION_SUPPORT     ));
    //Is CL_TRUE if the device implements error correction for the memories, caches, registers etc. in the device. Is CL_FALSE if the device does not implement error correction. This can be a requirement for certain clients of OpenCL.

    LINFO("CL_DEVICE_EXECUTION_CAPABILITIES: " << getInfo<cl_device_exec_capabilities>(CL_DEVICE_EXECUTION_CAPABILITIES     ));
    //Describes the execution capabilities of the device. This is a bit-field that describes one or more of the following values:
    //CL_EXEC_KERNEL - The OpenCL device can execute OpenCL kernels.
    //CL_EXEC_NATIVE_KERNEL - The OpenCL device can execute native kernels.
    //The mandated minimum capability is CL_EXEC_KERNEL.

    LINFO("CL_DEVICE_GLOBAL_MEM_CACHE_SIZE: " << getInfo<cl_ulong>(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE     ));
    //Size of global memory cache in bytes.

    LINFO("CL_DEVICE_GLOBAL_MEM_CACHE_TYPE: " << getInfo<cl_device_mem_cache_type>(CL_DEVICE_GLOBAL_MEM_CACHE_TYPE     ));
    //Type of global memory cache supported. Valid values are: CL_NONE, CL_READ_ONLY_CACHE, and CL_READ_WRITE_CACHE.

    LINFO("CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE: " << getInfo<cl_uint>(CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE     ));
    //Size of global memory cache line in bytes.

    cl_ulong globalMemSize;
    globalMemSize = getInfo<cl_ulong>(CL_DEVICE_GLOBAL_MEM_SIZE);
    LINFO("CL_DEVICE_GLOBAL_MEM_SIZE: " <<  globalMemSize << "b (" << (globalMemSize / (1024*1024))<< " mb)");
    //Size of global device memory in bytes.

    //LINFO("CL_DEVICE_HALF_FP_CONFIG: " << getInfo<cl_device_fp_config>(CL_DEVICE_HALF_FP_CONFIG     ));
    //Describes the OPTIONAL half precision floating-point capability of the OpenCL device. This is a bit-field that describes one or more of the following values:
        //* CL_FP_DENORM - denorms are supported.
        //* CL_FP_INF_NAN - INF and NaNs are supported.
        //* CL_FP_ROUND_TO_NEAREST - round to nearest even rounding mode supported.
        //* CL_FP_ROUND_TO_ZERO - round to zero rounding mode supported.
        //* CL_FP_ROUND_TO_INF - round to +ve and -ve infinity rounding modes supported.
        //* CP_FP_FMA - IEEE754-2008 fused multiply-add is supported.
    //The required minimum half precision floating-point capability as implemented by this extension is CL_FP_ROUND_TO_ZERO | CL_FP_ROUND_TO_INF | CL_FP_INF_NAN.

    cl_ulong localMemSize;
    localMemSize = getInfo<cl_ulong>(CL_DEVICE_LOCAL_MEM_SIZE);
    LINFO("CL_DEVICE_LOCAL_MEM_SIZE: " << localMemSize << "b (" << (localMemSize / 1024) << "kb)");
    //Size of local memory arena in bytes. The minimum value is 16 KB.

    LINFO("CL_DEVICE_LOCAL_MEM_TYPE: " << getInfo<cl_device_local_mem_type>(CL_DEVICE_LOCAL_MEM_TYPE));
    //Type of local memory supported. This can be set to CL_LOCAL implying dedicated local memory storage such as SRAM, or CL_GLOBAL.

    LINFO("CL_DEVICE_MAX_COMPUTE_UNITS: " << getInfo<cl_uint>(CL_DEVICE_MAX_COMPUTE_UNITS     ));
    //The number of parallel compute cores on the OpenCL device. The minimum value is 1.

    LINFO("CL_DEVICE_MAX_CONSTANT_ARGS: " << getInfo<cl_uint>(CL_DEVICE_MAX_CONSTANT_ARGS     ));
    //Max number of arguments declared with the __constant qualifier in a kernel. The minimum value is 8.

    LINFO("CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE: " << getInfo<cl_ulong>(CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE     ));
    //Max size in bytes of a constant buffer allocation. The minimum value is 64 KB.

    LINFO("CL_DEVICE_MAX_MEM_ALLOC_SIZE: " << getInfo<cl_ulong>(CL_DEVICE_MAX_MEM_ALLOC_SIZE     ));
    //Max size of memory object allocation in bytes. The minimum value is max (1/4th of CL_DEVICE_GLOBAL_MEM_SIZE, 128*1024*1024)

    LINFO("CL_DEVICE_MAX_PARAMETER_SIZE: " << getInfo<size_t>(CL_DEVICE_MAX_PARAMETER_SIZE     ));
    //Max size in bytes of the arguments that can be passed to a kernel. The minimum value is 256.

    LINFO("CL_DEVICE_MAX_WORK_GROUP_SIZE: " << getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE     ));
    //Maximum number of work-items in a work-group executing a kernel using the data parallel execution model. (Refer to clEnqueueNDRangeKernel). The minimum value is 1.

    cl_uint maxWorkItemDimensions = getInfo<cl_uint>(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);
    //LINFO("CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: " << maxWorkItemDimensions);
    //Maximum dimensions that specify the global and local work-item IDs used by the data parallel execution model. (Refer to clEnqueueNDRangeKernel). The minimum value is 3.

    std::vector<size_t> maxWorkItemSizes;
    maxWorkItemSizes.resize(maxWorkItemDimensions);
    LCL_ERROR(clGetDeviceInfo(id_, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*maxWorkItemDimensions, &maxWorkItemSizes[0], 0));
    std::stringstream temp;
    for(size_t i=0; i<maxWorkItemDimensions; ++i)
        temp << maxWorkItemSizes[i] << " ";
    LINFO("CL_DEVICE_MAX_WORK_ITEM_SIZES: " << temp.str());
    //Maximum number of work-items that can be specified in each dimension of the work-group to clEnqueueNDRangeKernel.
    //Returns n size_t entries, where n is the value returned by the query for CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS. The minimum value is (1, 1, 1).

    LINFO("CL_DEVICE_MEM_BASE_ADDR_ALIGN: " << getInfo<cl_uint>(CL_DEVICE_MEM_BASE_ADDR_ALIGN     ));
    //Describes the alignment in bits of the base address of any allocated memory object.

    LINFO("CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE: " << getInfo<cl_uint>(CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE     ));
    //The smallest alignment in bytes which can be used for any data type.

    LINFO("CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR: " << getInfo<cl_uint>(CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR));
    LINFO("CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT: " << getInfo<cl_uint>(CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT));
    LINFO("CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT: " << getInfo<cl_uint>(CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT));
    LINFO("CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG: " << getInfo<cl_uint>(CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG));
    LINFO("CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT: " << getInfo<cl_uint>(CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT));
    LINFO("CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE: " << getInfo<cl_uint>(CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE));
    //Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.
    //If the cl_khr_fp64 extension is not supported, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE must return 0.

    LINFO("CL_DEVICE_QUEUE_PROPERTIES: " << getInfo<cl_command_queue_properties>(CL_DEVICE_QUEUE_PROPERTIES     ));
    //Describes the command-queue properties supported by the device. This is a bit-field that describes one or more of the following values:
    //CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
    //CL_QUEUE_PROFILING_ENABLE
    //These properties are described in the table for clCreateCommandQueue. The mandated minimum capability is CL_QUEUE_PROFILING_ENABLE.

    LINFO("CL_DEVICE_SINGLE_FP_CONFIG: " << getInfo<cl_device_fp_config>(CL_DEVICE_SINGLE_FP_CONFIG     ));
    //Describes single precision floating-point capability of the device. This is a bit-field that describes one or more of the following values:
    //CL_FP_DENORM - denorms are supported
    //CL_FP_INF_NAN - INF and quiet NaNs are supported
    //CL_FP_ROUND_TO_NEAREST - round to nearest even rounding mode supported
    //CL_FP_ROUND_TO_ZERO - round to zero rounding mode supported
    //CL_FP_ROUND_TO_INF - round to +ve and -ve infinity rounding modes supported
    //CL_FP_FMA - IEEE754-2008 fused multiply-add is supported
    //The mandated minimum floating-point capability is CL_FP_ROUND_TO_NEAREST | CL_FP_INF_NAN.
}

bool Device::isExtensionSupported(std::string ext) const {
    if (extensions_.find(ext) != extensions_.end())
        return true;
    else
        return false;
}

bool Device::hasImageSupport() const {
    return imageSupport_;
}

tgt::ivec2 Device::getMaxImageSize2D() const  {
    return maxImageSize2D_;
}

tgt::ivec3 Device::getMaxImageSize3D() const {
    return maxImageSize3D_;
}

//-------------------------------------------------------------------------------------

const std::string Context::loggerCat_ = "voreen.OpenCL.Context";

Context::Context(const Device& device) {
    cl_int err;
    cl_device_id id = device.getId();
    //id_ = clCreateContext(0, 1, &id, contextCallback, 0, &err);
    id_ = clCreateContext(0, 1, &id, 0, 0, &err);
    LCL_ERROR(err);
}

Context::Context(const std::vector<ContextProperty>& properties, const Device& device) {
    cl_int err;
    cl_device_id id = device.getId();
    int numProps = static_cast<int>((properties.size() * 2)+1);
    cl_context_properties* props = new cl_context_properties[numProps];
    for(size_t i=0; i<properties.size(); ++i) {
        props[i*2] = properties[i].name_;
        props[(i*2)+1] = properties[i].value_;
    }
    props[numProps-1] = 0;     //terminate list

    //id_ = clCreateContext(props, 1, &id, contextCallback, 0, &err);
    id_ = clCreateContext(props, 1, &id, 0, 0, &err);
    LCL_ERROR(err);
    delete[] props;
}

Context::~Context() {
    LCL_ERROR(clReleaseContext(id_));
}

//template specialization for strings:
template<>
std::string Context::getInfo(cl_context_info info) const {
    size_t retSize;
    LCL_ERROR(clGetContextInfo(id_, info, 0, 0, &retSize));
    char* buffer = new char[retSize];
    LCL_ERROR(clGetContextInfo(id_, info, retSize, buffer, 0));
    std::string ret(buffer);
    delete[] buffer;
    return ret;
}

std::vector<ContextProperty> Context::generateGlSharingProperties() {
    std::vector<ContextProperty> ret;
#ifdef WIN32
    //Windows:
    HGLRC glCtx = wglGetCurrentContext();
    ret.push_back(ContextProperty(CL_GL_CONTEXT_KHR, (intptr_t) glCtx));
    ret.push_back(ContextProperty(CL_WGL_HDC_KHR, (intptr_t) wglGetCurrentDC()));
#elif defined(__APPLE__)
    //Apple:
    CGLContextObj curCGLContext = CGLGetCurrentContext();
    CGLShareGroupObj curCGLShareGroup = CGLGetShareGroup(curCGLContext);
    ret.push_back(ContextProperty(CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE, (cl_context_properties)curCGLShareGroup));
#else
    //Linux:
    GLXContext glCtx = glXGetCurrentContext();
    ret.push_back(ContextProperty(CL_GL_CONTEXT_KHR, (cl_context_properties ) glCtx));
    ret.push_back(ContextProperty(CL_GLX_DISPLAY_KHR, (intptr_t) glXGetCurrentDisplay()));
#endif
    return ret;
}

//---------------------------------------------------------

const std::string CommandQueue::loggerCat_ = "voreen.OpenCL.CommandQueue";

#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
CommandQueue::CommandQueue(const Context* context, const Device& device, cl_command_queue_properties properties) : properties_(properties) {
    cl_int err;
    id_ = clCreateCommandQueue(context->getId(), device.getId(), properties, &err);
    LCL_ERROR(err);
}
#else
CommandQueue::CommandQueue(const Context* context, const Device& device, const cl_queue_properties* properties) : properties_(properties) {
    cl_int err;
    id_ = clCreateCommandQueueWithProperties(context->getId(), device.getId(), properties, &err);
    LCL_ERROR(err);
}
#endif

CommandQueue::~CommandQueue() {
    clReleaseCommandQueue(id_);
}

void CommandQueue::enqueue(const Kernel* kernel, size_t globalWorkSize, size_t localWorkSize, Event* event /*= 0*/) {
    if (globalWorkSize < localWorkSize) {
        LWARNING("globalWorkSize (" << globalWorkSize << ") must be greater than localWorkSize (" << localWorkSize << ")");
    }
    if (globalWorkSize % localWorkSize != 0) {
        LWARNING("globalWorkSize (" << globalWorkSize << ") must be a multiple of localWorkSize (" << localWorkSize << ")");
    }
    LCL_ERROR(clEnqueueNDRangeKernel(id_, kernel->getId(), 1, 0, &globalWorkSize, &localWorkSize, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueue(const Kernel* kernel, size_t globalWorkSize, Event* event /*= 0*/) {
    LCL_ERROR(clEnqueueNDRangeKernel(id_, kernel->getId(), 1, 0, &globalWorkSize, 0, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueue(const Kernel* kernel, tgt::svec2 globalWorkSizes, tgt::svec2 localWorkSizes, Event* event /*= 0*/) {
    LCL_ERROR(clEnqueueNDRangeKernel(id_, kernel->getId(), 2, 0, &globalWorkSizes[0], &localWorkSizes[0], 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueue(const Kernel* kernel, tgt::svec2 globalWorkSizes, Event* event /*= 0*/) {
    LCL_ERROR(clEnqueueNDRangeKernel(id_, kernel->getId(), 2, 0, &globalWorkSizes[0], 0, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueue(const Kernel* kernel, tgt::svec3 globalWorkSizes, tgt::svec3 localWorkSizes, Event* event /*= 0*/) {
    LCL_ERROR(clEnqueueNDRangeKernel(id_, kernel->getId(), 3, 0, &globalWorkSizes[0], &localWorkSizes[0], 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueue(const Kernel* kernel, tgt::svec3 globalWorkSizes, Event* event /*= 0*/) {
    LCL_ERROR(clEnqueueNDRangeKernel(id_, kernel->getId(), 3, 0, &globalWorkSizes[0], 0, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueue(const Kernel* kernel, const std::vector<size_t>& globalWorkSizes, const std::vector<size_t>& localWorkSizes, Event* event /*= 0*/) {
    tgtAssert(globalWorkSizes.size() == localWorkSizes.size(), "Dimension mismatch!");
    LCL_ERROR(clEnqueueNDRangeKernel(id_, kernel->getId(), static_cast<cl_uint>(globalWorkSizes.size()),
        0, &globalWorkSizes[0], &localWorkSizes[0], 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueueReadBuffer(const Buffer* buffer, void* data, bool blocking, Event* event /*= 0*/) {
    LCL_ERROR(clEnqueueReadBuffer(id_, buffer->getId(), blocking, 0, buffer->getSize(), data, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueueReadBuffer(const Buffer* buffer, size_t byteOffset, size_t numBytes, void* data, bool blocking /*= true*/, Event* event /*= 0*/) {
    tgtAssert(byteOffset < buffer->getSize(), "offset outside buffer");
    tgtAssert(byteOffset+numBytes <= buffer->getSize(), "offset+numBytes outside buffer");

    LCL_ERROR(clEnqueueReadBuffer(id_, buffer->getId(), blocking, byteOffset, numBytes, data, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueueWriteBuffer(const Buffer* buffer, void* data, bool blocking, Event* event /*= 0*/) {
    LCL_ERROR(clEnqueueWriteBuffer(id_, buffer->getId(), blocking, 0, buffer->getSize(), data, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueueWriteBuffer(const Buffer* buffer, size_t byteOffset, size_t numBytes, void* data, bool blocking /*= true*/, Event* event /*= 0*/ ) {
    tgtAssert(byteOffset < buffer->getSize(), "offset outside buffer");
    tgtAssert(byteOffset+numBytes <= buffer->getSize(), "offset+numBytes outside buffer");

    LCL_ERROR(clEnqueueWriteBuffer(id_, buffer->getId(), blocking, byteOffset, numBytes, data, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueueCopyBuffer(const Buffer &src,
                        const Buffer &dst,
                        size_t src_offset, size_t dst_offset, size_t size, Event* event /*= 0*/)
{
    LCL_ERROR(
        clEnqueueCopyBuffer(id_, src.getId(), dst.getId(), src_offset, dst_offset, size, 0, 0, event ? &event->id_ : 0)
    );
}

void CommandQueue::enqueueCopyImage(const ImageObject2D &src,
                           const ImageObject2D &dst,
                           const size_t src_origin[3],
                           const size_t dst_origin[3],
                           const size_t region[3],
                           Event* event /*= 0*/)
{
    LCL_ERROR(
        clEnqueueCopyImage(id_, src.getId(), dst.getId(), (const size_t *)src_origin,
                            (const size_t *)dst_origin, (const size_t *)region,
                            0, 0, event ? &event->id_ : 0)
    );
}

void CommandQueue::enqueueWriteImage(const ImageObject2D &img,
                        tgt::Texture *tex,
                        const size_t origin[3],
                        const size_t region[3],
                        bool blocking,
                        Event* event /*= 0*/)
{
    LCL_ERROR(
        clEnqueueWriteImage(id_, img.getId(), blocking, (const size_t *)origin,
                            (const size_t *)region, 0, 0, tex->getCpuTextureData(), 0, 0,
                            event ? &event->id_ : 0)
    );
}

void CommandQueue::enqueueReadImage(const ImageObject2D &img,
                       tgt::Texture *tex,
                       const size_t origin[3],
                       const size_t region[3],
                       bool blocking,
                       Event* event /*= 0*/)
{
    LCL_ERROR(
        clEnqueueReadImage(id_, img.getId(), blocking, (const size_t *)origin,
                            (const size_t *)region, 0, 0, tex->getCpuTextureData(), 0, 0,
                            event ? &event->id_ : 0)
    );
}

void CommandQueue::enqueueReadImage(const ImageObject2D &img,
                       tgt::Texture *tex,
                       bool blocking,
                       Event* event /*= 0*/)
{
    size_t region[] = {static_cast<size_t>(tex->getDimensions().x), static_cast<size_t>(tex->getDimensions().y), 1};
    size_t origin[] = {0, 0, 0};

    LCL_ERROR(
        clEnqueueReadImage(id_, img.getId(), blocking, (const size_t *)origin,
                            (const size_t *)region, 0, 0, tex->getCpuTextureData(), 0, 0,
                            event ? &event->id_ : 0)
    );
}

void CommandQueue::enqueueCopyImageToBuffer(const ImageObject2D &src,
                                             const Buffer &dst,
                                             const size_t origin[3],
                                             const size_t region[3],
                                             size_t   dst_offset /* = 0 */,
                                             Event* event /*= 0*/)
{
    LCL_ERROR(
        clEnqueueCopyImageToBuffer(id_, src.getId(), dst.getId(), origin, region, dst_offset, 0, 0, event ? &event->id_ : 0)
    );
}

void CommandQueue::enqueueCopyBufferToImage(const Buffer &src,
                                             const ImageObject2D &dst,
                                             size_t src_offset,
                                             const size_t dst_origin[3],
                                             const size_t region[3],
                                             Event* event /*= 0*/)
{
    LCL_ERROR(
        clEnqueueCopyBufferToImage(id_, src.getId(), dst.getId(), src_offset, dst_origin, region, 0, 0, event ? &event->id_ : 0)
    );
}

void CommandQueue::enqueueAcquireGLObject(const MemoryObject* obj, Event* event /*= 0*/) {
    cl_mem mem = obj->getId();
    LCL_ERROR(clEnqueueAcquireGLObjects(id_, 1, &mem, 0, 0, event ? &event->id_ : 0));
}

void CommandQueue::enqueueReleaseGLObject(const MemoryObject* obj, Event* event /*= 0*/) {
    cl_mem mem = obj->getId();
    LCL_ERROR(clEnqueueReleaseGLObjects(id_, 1, &mem, 0, 0, event ? &event->id_ : 0));
}

#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
void CommandQueue::enqueueMarker(Event* event) {
    LCL_ERROR(clEnqueueMarker(id_, event ? &event->id_ : 0));
}

void CommandQueue::enqueueBarrier() {
    LCL_ERROR(clEnqueueBarrier(id_));
}

void CommandQueue::enqueueWaitForEvent(Event* event) {
    cl_event e = event->getId();
    LCL_ERROR(clEnqueueWaitForEvents(id_, 1, &e));
}

void CommandQueue::enqueueWaitForEvents(const std::vector<Event*>& events) {
    cl_event* e = new cl_event[events.size()];
    for (size_t i = 0; i<events.size(); ++i)
        e[i] = events[i]->getId();

    LCL_ERROR(clEnqueueWaitForEvents(id_, static_cast<cl_uint>(events.size()), e));
    delete[] e;
}
#else
void CommandQueue::enqueueMarker(const std::vector<Event*>& waitlist, Event* event) {
    cl_event* e = new cl_event[waitlist.size()];
    for (size_t i = 0; i<waitlist.size(); ++i)
        e[i] = waitlist[i]->getId();

    LCL_ERROR(clEnqueueMarkerWithWaitList(id_, static_cast<cl_uint>(waitlist.size()), e, event ? &event->id_ : 0));
    delete[] e;
}

void CommandQueue::enqueueBarrier(const std::vector<Event*>& waitlist, Event* event) {
    cl_event* e = new cl_event[waitlist.size()];
    for (size_t i = 0; i<waitlist.size(); ++i)
        e[i] = waitlist[i]->getId();

    LCL_ERROR(clEnqueueBarrierWithWaitList(id_, static_cast<cl_uint>(waitlist.size()), e, event ? &event->id_ : 0));
    delete[] e;
}
#endif

//template specialization for strings:
template<>
std::string CommandQueue::getInfo(cl_command_queue_info info) const {
    size_t retSize;
    LCL_ERROR(clGetCommandQueueInfo(id_, info, 0, 0, &retSize));
    char* buffer = new char[retSize];
    LCL_ERROR(clGetCommandQueueInfo(id_, info, retSize, buffer, 0));
    std::string ret(buffer);
    delete[] buffer;
    return ret;
}

//---------------------------------------------------------

const std::string Program::loggerCat_ = "voreen.OpenCL.Program";

Program::Program(const Context& context) : id_(0), context_(&context), header_("") {
}

Program::Program(const Context* context) : id_(0), context_(context), header_("") {
}

Program::~Program() {
    releaseProgram();
}

void Program::releaseProgram() {
    clearKernels();
    if (id_)
        LCL_ERROR(clReleaseProgram(id_));
    id_ = 0;
}

void Program::createWithSource() {
    if(source_.empty()) {
        LERROR("Trying to create program without source!");
        return;
    }

    cl_int err;
    size_t numSources = source_.size() + 1;
    const char** s = new const char*[numSources];
    size_t* l = new size_t[numSources];

    std::string header;
    if(sizeof(size_t) == 8)
        header = "typedef unsigned long vrn_size_t;\n";
    else
        header = "typedef unsigned int vrn_size_t;\n";

    header += header_;

    s[0] = header.c_str();
    l[0] = header.length();

    for(size_t i=1; i<numSources; ++i)  {
        s[i] = source_[i-1].c_str();
        l[i] = source_[i-1].length();
    }

    id_ = clCreateProgramWithSource(context_->getId(), static_cast<cl_uint>(numSources), s, l, &err);
    LCL_ERROR(err);

    delete[] l;
    delete[] s;
}

void Program::setSource(const std::string& source) {
    if(id_)
        releaseProgram();

    filenames_.clear();
    source_.clear();

    source_.push_back(source);
    createWithSource();
}

void Program::setSource(const std::vector<std::string>& source) {
    if(id_)
        releaseProgram();

    filenames_.clear();
    source_.clear();
    source_ = source;
    createWithSource();
}

bool Program::loadSource(const std::string& filename) {
    tgt::File* file = FileSys.open(filename);

    // check if file is open
    if (!file || !file->isOpen()) {
        LERROR("File not found: " << filename);
        return false;
    }

    source_.clear();

    source_.push_back(file->getAsString());
    file->close();
    delete file;
    filenames_.clear();
    filenames_.push_back(filename);
    createWithSource();
    return true;
}

bool Program::loadSource(const std::vector<std::string>& filenames) {
    source_.clear();

    for(size_t i=0; i<filenames.size(); ++i) {
        tgt::File* file = FileSys.open(filenames[i]);

        // check if file is open
        if (!file || !file->isOpen()) {
            LERROR("File not found: " << filenames[i]);
            delete file;
            return false;
        }
        source_.push_back(file->getAsString());
        file->close();
        delete file;
    }
    filenames_.clear();
    filenames_ = filenames;
    createWithSource();
    return true;
}

std::string Program::getHeader() const {
    return header_;
}

void Program::setHeader(const std::string& header) {
    header_ = header;
}

bool Program::reload() {
    if(!filenames_.empty()) {
        if(filenames_.size() == 1)
            return loadSource(filenames_[0]);
        else
            return loadSource(filenames_);
    }
    LWARNING("Trying to reload program that has not been loaded from disk!");
    return true;
}

bool Program::build() {
    if(id_ == 0)
        return false;

    //TODO: check if program needs to be re-created
    createWithSource();
    clearKernels();
    cl_int err = LCL_ERROR(clBuildProgram(id_, 0, 0, buildOptions_.c_str(), 0, 0));
    if(err == CL_SUCCESS) {
        LINFO("Successfully built program.");
        return true;
    }
    else {
        LWARNING("Failed to build program");
        return false;
    }
}

bool Program::build(const Device& device) {
    if(id_ == 0)
        return false;

    if(getBuildStatus(device) != CL_BUILD_NONE) {
        //to rebuild we need a new program
        createWithSource();
    }

    clearKernels();
    cl_device_id devId = device.getId();
    cl_int err = LCL_ERROR(clBuildProgram(id_, 1, &devId, buildOptions_.c_str(), 0, 0));
    if(err == CL_SUCCESS) {
        LINFO("Successfully built program for device: " << device.getName());
        createKernels();
        return true;
    }
    else {
        LWARNING("Failed to build program: " << statusToString(getBuildStatus(device)) << " (device: " << device.getName() << ")");
        LWARNING("Build Log:\n" << getBuildLog(device));
        return false;
    }
}

bool Program::build(const std::vector<Device>& devices) {
    if(!id_)
        return false;

    //TODO: check if program needs to be re-created
    createWithSource();
    clearKernels();
    cl_device_id* devIds = new cl_device_id[devices.size()];
    for(size_t i=0; i<devices.size(); ++i)
        devIds[i] = devices[i].getId();

    cl_int err = LCL_ERROR(clBuildProgram(id_, static_cast<cl_uint>(devices.size()), devIds, buildOptions_.c_str(), 0, 0));
    if(err == CL_SUCCESS) {
        std::string devs = "";
        for(size_t i=0; i<devices.size(); ++i)
            devs += devices[i].getName() + ";";
        LINFO("Succesfully build program for devices: " << devs);
        return true;
    }
    else {
        for(size_t i=0; i<devices.size(); ++i) {
            LWARNING("Failed to build program: " << statusToString(getBuildStatus(devices[i])) << " (device: " << devices[i].getName() << ")");
            LWARNING("Build Log:\n" << getBuildLog(devices[i]));
        }
        return false;
    }
    delete[] devIds;
}

std::string Program::statusToString(cl_build_status status) {
    switch(status) {
        case CL_BUILD_NONE:    return "CL_BUILD_NONE";
        case CL_BUILD_ERROR: return "CL_BUILD_ERROR";
        case CL_BUILD_SUCCESS: return "CL_BUILD_SUCCESS";
        case CL_BUILD_IN_PROGRESS: return "CL_BUILD_IN_PROGRESS";
        default:
            return "unknown";
    }
}

std::string Program::getBuildLog(const Device& device) const {
    return getBuildInfo<std::string>(device, CL_PROGRAM_BUILD_LOG);
}

//template specialization for strings:
template<>
std::string Program::getBuildInfo(const Device& device, cl_program_build_info info) const {
    size_t retSize;
    LCL_ERROR(clGetProgramBuildInfo(id_, device.getId(), info, 0, 0, &retSize));
    char* buffer = new char[retSize+1];
    LCL_ERROR(clGetProgramBuildInfo(id_, device.getId(), info, retSize, buffer, 0));
    buffer[retSize] = '\0';
    std::string s(buffer);
    delete[] buffer;
    return s;
}

template<>
std::string Program::getInfo(cl_program_info info) const {
    size_t retSize;
    LCL_ERROR(clGetProgramInfo(id_, info, 0, 0, &retSize));
    char* buffer = new char[retSize];
    LCL_ERROR(clGetProgramInfo(id_, info, retSize, buffer, 0));
    std::string s(buffer);
    delete[] buffer;
    return s;
}

cl_build_status Program::getBuildStatus(const Device& device) const {
    return getBuildInfo<cl_build_status>(device, CL_PROGRAM_BUILD_STATUS);
}

std::string Program::getBuildOptions(const Device& device) const {
    return getBuildInfo<std::string>(device, CL_PROGRAM_BUILD_OPTIONS);
}
std::map<std::string, Kernel*> Program::getCurrentKernels() const {
    std::lock_guard<std::mutex> guard(kernelsMutex_);

    std::map<std::string, Kernel*> res;
    for(auto& kv : kernels_) {
       res[kv.first.first] = kv.second.get();
    }
    return res;
}

Kernel* Program::getKernel(const std::string& name) {
    std::lock_guard<std::mutex> guard(kernelsMutex_);

    std::pair<std::string, std::thread::id> kernelIdentifier(name, std::this_thread::get_id());

    auto kernelIt = kernels_.find(kernelIdentifier);
    if(kernelIt != kernels_.end()) {
        return kernelIt->second.get();
    }
    else {
        if(id_ == 0)
            return nullptr;

        cl_int err;
        cl_kernel kernel = clCreateKernel(id_, name.c_str(), &err);
        LCL_ERROR(err);
        if((kernel != 0) && (err == CL_SUCCESS)) {
            auto kernelPtr = tgt::make_unique<Kernel>(kernel);
            auto res = kernelPtr.get();
            kernels_[kernelIdentifier] = std::move(kernelPtr);
            return res;
        }
        else
            return nullptr;
    }
}

void Program::createKernels() {
    clearKernels();

    cl_uint numKernels;
    LCL_ERROR(clCreateKernelsInProgram(id_, 0, 0, &numKernels));
    LINFO("Found " << numKernels << " kernels in program.");
    cl_kernel* kernels = new cl_kernel[numKernels];
    LCL_ERROR(clCreateKernelsInProgram(id_, numKernels, kernels, 0));

    for(size_t i=0; i<numKernels; ++i) {
        auto k = tgt::make_unique<Kernel>(kernels[i]);
        std::string name = k->getName();
        std::pair<std::string, std::thread::id> kernelIdentifier(name, std::this_thread::get_id());
        LINFO("Found kernel '" << name << "' with " << k->getNumArgs() << " arguments." );
        kernels_[kernelIdentifier] = std::move(k);
    }
    delete[] kernels;
}

void Program::clearKernels() {
    std::lock_guard<std::mutex> guard(kernelsMutex_);

    kernels_.clear();
}

//---------------------------------------------------------

const std::string Kernel::loggerCat_ = "voreen.OpenCL.Kernel";

Kernel::Kernel(cl_kernel id) : id_(id) {
}

Kernel::~Kernel() {
    LCL_ERROR(clReleaseKernel(id_));
}

std::string Kernel::getName() const {
   return getInfo<std::string>(CL_KERNEL_FUNCTION_NAME);
}

//template specialization for strings:
template<>
std::string Kernel::getInfo(cl_kernel_info info) const {
    size_t retSize;
    LCL_ERROR(clGetKernelInfo(id_, info, 0, 0, &retSize));
    char* buffer = new char[retSize];
    LCL_ERROR(clGetKernelInfo(id_, info, retSize, buffer, 0));
    std::string ret(buffer);
    delete[] buffer;
    return ret;
}

template<>
std::string Kernel::getWorkGroupInfo(cl_kernel_work_group_info info, const Device& device) const {
    size_t retSize;
    LCL_ERROR(clGetKernelWorkGroupInfo(id_, device.getId(), info, 0, 0, &retSize));
    char* buffer = new char[retSize];
    LCL_ERROR(clGetKernelWorkGroupInfo(id_, device.getId(), info, retSize, buffer, 0));
    std::string ret(buffer);
    delete[] buffer;
    return ret;
}

cl_uint Kernel::getNumArgs() const {
    return getInfo<cl_uint>(CL_KERNEL_NUM_ARGS);
}

cl_uint Kernel::getWorkGroupSize(const Device& device) const {
    return getWorkGroupInfo<cl_uint>(CL_KERNEL_WORK_GROUP_SIZE, device);
}

//---------------------------------------------------------

MemoryObject::MemoryObject() : id_(0) {
}

MemoryObject::~MemoryObject() {
    if(id_)
        LCL_ERROR(clReleaseMemObject(id_));
}

//---------------------------------------------------------

Buffer::Buffer(const Context* context, cl_mem_flags flags, size_t size, const void* hostPtr,  cl_int* err) : MemoryObject(), size_(size) {
    if(err) {
        //FIXME: const_cast (stefan)
        id_ = clCreateBuffer(context->getId(), flags, size, const_cast<void*>(hostPtr), err);
    }
    else {
        cl_int tmpErr;
        //FIXME: const_cast (stefan)
        id_ = clCreateBuffer(context->getId(), flags, size, const_cast<void*>(hostPtr), &tmpErr);
        LCL_ERROR(tmpErr);
    }
}

//Buffer::Buffer(const Context *context, cl_mem_flags flags, size_t size, void* hostPtr) : MemoryObject(), size_(size) {
//    cl_int err;
//    id_ = clCreateBuffer(context->getId(), flags, size, hostPtr, &err);
//    LCL_ERROR(err);
//}

Buffer::~Buffer() {
}

//---------------------------------------------------------

SharedTexture::SharedTexture(const Context* context, cl_mem_flags flags, tgt::Texture* tex) : MemoryObject(), tex_(tex) {
    //check if tex really is a texture:
    //#ifdef DEBUG this?
    if(!glIsTexture(tex->getId())) {
        LERRORC("voreen.OpenCL.SharedTexture", "tex is not a texture!");
        return;
    }

#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    cl_int err;
    switch (tex->getType()) {
    case GL_TEXTURE_3D:
        id_ = clCreateFromGLTexture3D(context->getId(), flags, tex->getType(), 0, tex->getId(), &err);
        LCL_ERROR(err);
        break;
    case GL_TEXTURE_2D:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_X:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_Y:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_Z:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_X:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_Y:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_Z:
    #ifdef GL_TEXTURE_RECTANGLE
    case GL_TEXTURE_RECTANGLE:
    #else
    case GL_TEXTURE_RECTANGLE_ARB:
    #endif
        id_ = clCreateFromGLTexture2D(context->getId(), flags, tex->getType(), 0, tex->getId(), &err);
        LCL_ERROR(err);
        break;
    default:
        LERRORC("voreen.OpenCL.SharedTexture", "Unknown texture type!");
    }
#else
    cl_int err;
    id_ = clCreateFromGLTexture(context->getId(), flags, tex->getType(), 0, tex->getId(), &err);
    LCL_ERROR(err);
#endif
}

SharedTexture::~SharedTexture() {
}

//---------------------------------------------------------

ImageFormat::ImageFormat(cl_channel_order order, cl_channel_type type) {
    image_channel_order = order;
    image_channel_data_type = type;
}


//---------------------------------------------------------

ImageObject2D::ImageObject2D(const Context *context,
                             cl_mem_flags flags,
                             tgt::Texture *tex)
    : MemoryObject()
{
    if (!tex) {
        LERRORC("voreen.OpenCL.ImageObject2D", "tex is not a valid texture");
        return ;
    }
    if (!tex->getCpuTextureData()) {
      tex->downloadTexture();
    }
    if (!tex->getCpuTextureData()) {
        LERRORC("voreen.OpenCL.ImageObject2D", "tex is not a valid texture");
        return ;
    }

    cl_image_format image_format;
    tgt::ivec3 dims = tex->getDimensions();
    GLenum texDataType = tex->getGLDataType();
    GLint  texFormat = tex->getGLFormat();
    //int bpp = 1;

    switch(texDataType) {
    case GL_BYTE:
        image_format.image_channel_data_type = CL_SNORM_INT8;
        break;
    case GL_UNSIGNED_BYTE:
        image_format.image_channel_data_type = CL_UNORM_INT8;
        break;
    case GL_SHORT:
        image_format.image_channel_data_type = CL_SNORM_INT16;
        break;
    case GL_UNSIGNED_SHORT:
        image_format.image_channel_data_type = CL_UNORM_INT16;
        break;
    case GL_INT:
        image_format.image_channel_data_type = CL_SIGNED_INT32;
        break;
    case GL_UNSIGNED_INT:
        image_format.image_channel_data_type = CL_UNSIGNED_INT32;
        break;
    case GL_FLOAT:
        image_format.image_channel_data_type = CL_FLOAT;
        break;
    }

    switch(texFormat) {
    case GL_RGBA:
        image_format.image_channel_order = CL_RGBA;
        break;
    case GL_RGB:
        image_format.image_channel_order = CL_RGB;
        break;
    case GL_DEPTH_COMPONENT:
    case GL_ALPHA:
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    case GL_LUMINANCE:
#else
    case GL_RED:
#endif
        image_format.image_channel_order = CL_INTENSITY;
        break;
    }

#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    cl_int err;
    id_ = clCreateImage2D(context->getId(),
                          flags,
                          &image_format,
                          dims.x, dims.y,
                          0,
                          tex->getCpuTextureData(),
                          &err);
    LCL_ERROR(err);
#else
    cl_image_desc desc;
    desc.image_type = CL_MEM_OBJECT_IMAGE2D;
    desc.image_width = dims.x;
    desc.image_height = dims.y;
    desc.image_row_pitch = 0;
    desc.image_slice_pitch = 0;
    desc.num_mip_levels = 0;
    desc.num_samples = 0;
    desc.buffer = nullptr;

    cl_int err;
    id_ = clCreateImage(context->getId(), flags, &image_format, &desc, tex->getCpuTextureData(), &err);
    LCL_ERROR(err);
#endif
}

ImageObject2D::ImageObject2D(const Context *context, cl_mem_flags flags,
                             ImageFormat format,
                             size_t width,
                             size_t height,
                             size_t row_pitch /* = 0 */,
                             void *host_ptr /* = NULL */)
{
#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    cl_int err;
    id_ = clCreateImage2D(context->getId(), flags, &format, width, height, row_pitch, host_ptr, &err);
    LCL_ERROR(err);
#else
    cl_image_desc desc;
    desc.image_type = CL_MEM_OBJECT_IMAGE2D;
    desc.image_width = width;
    desc.image_height = height;
    desc.image_row_pitch = row_pitch;
    desc.image_slice_pitch = 0;
    desc.num_mip_levels = 0;
    desc.num_samples = 0;
    desc.buffer = nullptr;

    cl_int err;
    id_ = clCreateImage(context->getId(), flags, &format, &desc, host_ptr, &err);
    LCL_ERROR(err);
#endif
}

ImageObject2D::~ImageObject2D() {
}

//---------------------------------------------------------

ImageObject3D::ImageObject3D(const Context* context, cl_mem_flags flags, const VolumeRAM* vol)
    : MemoryObject()
{
    if(!vol || !vol->getData()) {
        LERRORC("voreen.OpenCL.ImageObject3D", "vol is not a valid volume!");
        return;
    }

    cl_image_format volume_format;

    const std::type_info& volumeType(typeid(*vol));

    // MOST OF THESE HAVE NOT YET BEEN TESTED, USE CAREFULLY!
    // VolumeRAM_UIntX
    if (volumeType == typeid(VolumeRAM_UInt8)) {
        volume_format.image_channel_order = CL_INTENSITY;
        volume_format.image_channel_data_type = CL_UNORM_INT8;
    }
    else if (volumeType == typeid(VolumeRAM_UInt16)) {
        volume_format.image_channel_order = CL_INTENSITY;
        volume_format.image_channel_data_type = CL_UNORM_INT16;
    }
    else if (volumeType == typeid(VolumeRAM_UInt32)) {
        volume_format.image_channel_order = CL_INTENSITY;
        volume_format.image_channel_data_type = CL_UNSIGNED_INT32;
    }
    // VolumeRAM_IntX
    else if (volumeType == typeid(VolumeRAM_Int8)) {
        volume_format.image_channel_order = CL_INTENSITY;
        volume_format.image_channel_data_type = CL_SNORM_INT8;
    }
    else if (volumeType == typeid(VolumeRAM_Int16)) {
        volume_format.image_channel_order = CL_INTENSITY;
        volume_format.image_channel_data_type = CL_SNORM_INT16;
    }
    else if (volumeType == typeid(VolumeRAM_Int32)) {
        volume_format.image_channel_order = CL_INTENSITY;
        volume_format.image_channel_data_type = CL_SIGNED_INT32;
    }
    // VolumeRAM_Float
    else if (volumeType == typeid(VolumeRAM_Float)) {
        volume_format.image_channel_order = CL_INTENSITY;
        volume_format.image_channel_data_type = CL_FLOAT;
    }
    // VolumeRAM_3x with int8 types
    else if (volumeType == typeid(VolumeRAM_3xUInt8)) {
        volume_format.image_channel_order = CL_RGB;
        volume_format.image_channel_data_type = CL_UNORM_INT8;
    }
    else if (volumeType == typeid(VolumeRAM_3xInt8)) {
        volume_format.image_channel_order = CL_RGB;
        volume_format.image_channel_data_type = CL_UNORM_INT16;
    }
    // VolumeRAM_4x with int8 types
    else if (volumeType == typeid(VolumeRAM_4xUInt8)) {
        volume_format.image_channel_order = CL_RGBA;
        volume_format.image_channel_data_type = CL_UNORM_INT8;
    }
    else if (volumeType == typeid(VolumeRAM_4xInt8)) {
        volume_format.image_channel_order = CL_RGBA;
        volume_format.image_channel_data_type = CL_SNORM_INT8;
    }
    // VolumeRAM_3x with int16 types
    else if (volumeType == typeid(VolumeRAM_3xUInt16)) {
        volume_format.image_channel_order = CL_RGB;
        volume_format.image_channel_data_type = CL_UNORM_INT16;
    }
    else if (volumeType == typeid(VolumeRAM_3xInt16)) {
        volume_format.image_channel_order = CL_RGB;
        volume_format.image_channel_data_type = CL_SNORM_INT16;
    }
    // VolumeRAM_4x with int16 types
    else if (volumeType == typeid(VolumeRAM_4xUInt16)) {
        volume_format.image_channel_order = CL_RGBA;
        volume_format.image_channel_data_type = CL_UNORM_INT16;
    }
    else if (volumeType == typeid(VolumeRAM_4xInt16)) {
        volume_format.image_channel_order = CL_RGBA;
        volume_format.image_channel_data_type = CL_SNORM_INT16;
    }
    // VolumeRAM_3x with real types
    else if (volumeType == typeid(VolumeRAM_3xFloat)) {
        volume_format.image_channel_order = CL_RGB;
        volume_format.image_channel_data_type = CL_FLOAT;
    }
    // VolumeRAM_4x with real types
    else if (volumeType == typeid(VolumeRAM_4xFloat)) {
        volume_format.image_channel_order = CL_RGBA;
        volume_format.image_channel_data_type = CL_FLOAT;
    }
    // -> not found
    else {
        LERRORC("voreen.OpenCL.ImageObject", "Unknown texture type!");
    }

#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    tgt::ivec3 dims = vol->getDimensions();

    cl_int err;
    /*id_ = clCreateImage3D(context->getId(), flags, &volume_format,
                            dims.x, dims.y, dims.z,
                            dims.x * vol->getBytesPerVoxel(),
                            dims.x * dims.y * vol->getBytesPerVoxel(),
                            vol->getData(), &err);*/
    id_ = clCreateImage3D(context->getId(), flags, &volume_format,
                          dims.x, dims.y, dims.z,
                          0,
                          0,
                          const_cast<void*>(vol->getData()), &err);
    LCL_ERROR(err);
#else
    cl_image_desc desc;
    desc.image_type = CL_MEM_OBJECT_IMAGE3D;
    desc.image_width  = vol->getDimensions().x;
    desc.image_height = vol->getDimensions().y;
    desc.image_depth  = vol->getDimensions().z;
    desc.image_row_pitch = 0;
    desc.image_slice_pitch = 0;
    desc.num_mip_levels = 0;
    desc.num_samples = 0;
    desc.buffer = nullptr;

    cl_int err;
    id_ = clCreateImage(context->getId(), flags, &volume_format, &desc, const_cast<void*>(vol->getData()), &err);
    LCL_ERROR(err);
#endif
}

ImageObject3D::~ImageObject3D() {
}

//---------------------------------------------------------
#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
Sampler::Sampler(const Context* context, cl_bool normalized_coords, cl_addressing_mode addressing_mode, cl_filter_mode filter_mode) {
    cl_int err;
    id_ = clCreateSampler(context->getId(), normalized_coords, addressing_mode, filter_mode, &err);
    LCL_ERROR(err);
}
#else
Sampler::Sampler(const Context* context, const cl_sampler_properties* properties) {
    cl_int err;
    id_ = clCreateSamplerWithProperties(context->getId(), properties, &err);
    LCL_ERROR(err);
}
#endif

Sampler::~Sampler() {
    clReleaseSampler(id_);
}

//---------------------------------------------------------

std::string clErrorToString(cl_int err) {
    std::string msg;
    switch(err) {
        case(CL_SUCCESS): msg = "CL_SUCCESS"; break;
        case(CL_BUILD_PROGRAM_FAILURE): msg = "CL_BUILD_PROGRAM_FAILURE"; break;
        case(CL_COMPILER_NOT_AVAILABLE): msg = "CL_COMPILER_NOT_AVAILABLE"; break;
        case(CL_DEVICE_NOT_AVAILABLE): msg = "CL_DEVICE_NOT_AVAILABLE"; break;
        case(CL_DEVICE_NOT_FOUND): msg = "CL_DEVICE_NOT_FOUND"; break;
        case(CL_IMAGE_FORMAT_MISMATCH): msg = "CL_IMAGE_FORMAT_MISMATCH"; break;
        case(CL_IMAGE_FORMAT_NOT_SUPPORTED): msg = "CL_IMAGE_FORMAT_NOT_SUPPORTED"; break;
        case(CL_INVALID_ARG_INDEX): msg = "CL_INVALID_ARG_INDEX"; break;
        case(CL_INVALID_ARG_SIZE): msg = "CL_INVALID_ARG_SIZE"; break;
        case(CL_INVALID_ARG_VALUE): msg = "CL_INVALID_ARG_VALUE"; break;
        case(CL_INVALID_BINARY): msg = "CL_INVALID_BINARY"; break;
        case(CL_INVALID_BUFFER_SIZE): msg = "CL_INVALID_BUFFER_SIZE"; break;
        case(CL_INVALID_BUILD_OPTIONS): msg = "CL_INVALID_BUILD_OPTIONS"; break;
        case(CL_INVALID_COMMAND_QUEUE): msg = "CL_INVALID_COMMAND_QUEUE"; break;
        case(CL_INVALID_CONTEXT): msg = "CL_INVALID_CONTEXT"; break;
        case(CL_INVALID_DEVICE): msg = "CL_INVALID_DEVICE"; break;
        case(CL_INVALID_DEVICE_TYPE): msg = "CL_INVALID_DEVICE_TYPE"; break;
        case(CL_INVALID_EVENT): msg = "CL_INVALID_EVENT"; break;
        case(CL_INVALID_EVENT_WAIT_LIST): msg = "CL_INVALID_EVENT_WAIT_LIST"; break;
        case(CL_INVALID_GLOBAL_OFFSET): msg = "CL_INVALID_GLOBAL_OFFSET"; break;
//        case(CL_INVALID_GLOBAL_WORK_SIZE): msg = "CL_INVALID_GLOBAL_WORK_SIZE"; break;
        case(CL_INVALID_GL_OBJECT): msg = "CL_INVALID_GL_OBJECT"; break;
        case(CL_INVALID_HOST_PTR): msg = "CL_INVALID_HOST_PTR"; break;
        case(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR): msg = "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"; break;
        case(CL_INVALID_IMAGE_SIZE): msg = "CL_INVALID_IMAGE_SIZE"; break;
        case(CL_INVALID_KERNEL): msg = "CL_INVALID_KERNEL"; break;
        case(CL_INVALID_KERNEL_ARGS): msg = "CL_INVALID_KERNEL_ARGS"; break;
        case(CL_INVALID_KERNEL_DEFINITION): msg = "CL_INVALID_KERNEL_DEFINITION"; break;
        case(CL_INVALID_KERNEL_NAME): msg = "CL_INVALID_KERNEL_NAME"; break;
        case(CL_INVALID_MEM_OBJECT): msg = "CL_INVALID_MEM_OBJECT"; break;
        case(CL_INVALID_MIP_LEVEL): msg = "CL_INVALID_MIP_LEVEL"; break;
        case(CL_INVALID_OPERATION): msg = "CL_INVALID_OPERATION"; break;
        case(CL_INVALID_PLATFORM): msg = "CL_INVALID_PLATFORM"; break;
        case(CL_INVALID_PROGRAM): msg = "CL_INVALID_PROGRAM"; break;
        case(CL_INVALID_PROGRAM_EXECUTABLE): msg = "CL_INVALID_PROGRAM_EXECUTABLE"; break;
        case(CL_INVALID_QUEUE_PROPERTIES): msg = "CL_INVALID_QUEUE_PROPERTIES"; break;
        case(CL_INVALID_SAMPLER): msg = "CL_INVALID_SAMPLER"; break;
        case(CL_INVALID_VALUE): msg = "CL_INVALID_VALUE"; break;
        case(CL_INVALID_WORK_DIMENSION): msg = "CL_INVALID_WORK_DIMENSION"; break;
        case(CL_INVALID_WORK_GROUP_SIZE): msg = "CL_INVALID_WORK_GROUP_SIZE"; break;
        case(CL_INVALID_WORK_ITEM_SIZE): msg = "CL_INVALID_WORK_ITEM_SIZE"; break;
        case(CL_MAP_FAILURE): msg = "CL_MAP_FAILURE"; break;
        case(CL_MEM_COPY_OVERLAP): msg = "CL_MEM_COPY_OVERLAP"; break;
        case(CL_MEM_OBJECT_ALLOCATION_FAILURE): msg = "CL_MEM_OBJECT_ALLOCATION_FAILURE"; break;
        case(CL_OUT_OF_HOST_MEMORY): msg = "CL_OUT_OF_HOST_MEMORY"; break;
        case(CL_OUT_OF_RESOURCES): msg = "CL_OUT_OF_RESOURCES"; break;
        case(CL_PROFILING_INFO_NOT_AVAILABLE): msg = "CL_PROFILING_INFO_NOT_AVAILABLE"; break;
        case(CL_INVALID_PROPERTY): msg = "CL_INVALID_PROPERTY"; break;
        default:
           msg = "unknown";
    }
    std::ostringstream tmp;
    tmp << msg << " (" << (int)err << ")";
    return tmp.str();
}

cl_int _lCLError(cl_int err, int line, const char* file) {
    if (err != CL_SUCCESS) {
        std::ostringstream tmp2, loggerCat;
        if (file) {
            tmp2 << " File: " << file << "@" << line;
            loggerCat << "cl-error:" << file << ':' << line;
        }
        LogMgr.log(loggerCat.str(), tgt::Error, clErrorToString(err), tmp2.str());
        //tgtAssert(false, "");
    }
    return err;
}

VolumeWriteBuffer::VolumeWriteBuffer(const Context* context, VolumeRAM* vol) : infoBuffer_(0), dataBuffer_(0), refCount_(0) {

    VolumeWriteBufferCL b;

    b.dimensions_ = tgt::Vector4<cl_uint>(vol->getDimensions(), 1);
    b.numChannels_ = (cl_uchar)(vol->getNumChannels());
    b.numBitsPerChannel_ = (cl_uchar)(vol->getBitsAllocated() / vol->getNumChannels());

    if(vol->elementRange().x < 0.f)
        b.isSigned_ = 1;
    else
        b.isSigned_ = 0;

    const std::type_info& volumeType(typeid(*vol));

    if (volumeType == typeid(VolumeRAM_Float) || volumeType == typeid(VolumeRAM_3xFloat) || volumeType == typeid(VolumeRAM_4xFloat))
        b.isFloat_ = 1;
    else
        b.isFloat_ = 0;

    b.data_ = 0;

    infoBuffer_ = new Buffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(VolumeWriteBuffer), &b);
    dataBuffer_ = new Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, vol->getNumBytes());

    volume_ = vol;
}


void fillGradientModesPropertyCL(StringOptionProperty& gradientMode) {
    gradientMode.addOption("none",                 "none"                  );
    gradientMode.addOption("forward-differences",  "Forward Differences"   );
    gradientMode.addOption("central-differences",  "Central Differences"   );
    gradientMode.addOption("sobel",                "Sobel"                 );
    gradientMode.addOption("filtered",             "Filtered"              );
    gradientMode.select("central-differences");
}

void fillShadingModesPropertyCL(StringOptionProperty& shadeMode) {
    shadeMode.addOption("none",                   "none"                   );
    shadeMode.addOption("phong-diffuse",          "Phong (Diffuse)"        );
    shadeMode.addOption("phong-specular",         "Phong (Specular)"       );
    shadeMode.addOption("phong-diffuse-ambient",  "Phong (Diffuse+Amb.)"   );
    shadeMode.addOption("phong-diffuse-specular", "Phong (Diffuse+Spec.)"  );
    shadeMode.addOption("phong",                  "Phong (Full)"           );
    shadeMode.addOption("toon",                   "Toon"                   );
    shadeMode.addOption("cook-torrance",          "Cook-Torrance"          );
    shadeMode.addOption("oren-nayar",             "Oren-Nayar"             );
    shadeMode.addOption("ward",                   "Ward (Isotropic)"       );
    shadeMode.select("phong");
}

void fillCompositingModesPropertyCL(StringOptionProperty& compositingMode) {
    compositingMode.addOption("dvr", "DVR");
    compositingMode.addOption("mip", "MIP");
    compositingMode.addOption("mida", "MIDA");
    compositingMode.addOption("iso", "ISO");
    compositingMode.addOption("fhp", "FHP");
    compositingMode.addOption("fhn", "FHN");
}

void fillClassificationModesPropertyCL(StringOptionProperty* prop) {
    prop->addOption("none", "none");
    prop->addOption("transfer-function", "Transfer Function");
    //prop->addOption("pre-integrated-gpu", "Pre-Integrated TF (GPU)");
    prop->addOption("pre-integrated-approximation", "Pre-integrated TF (Approximation)");
    prop->addOption("pre-integrated-cpu", "Pre-integrated TF (CPU)");
    prop->select("transfer-function");
}

std::string getGradientDefineCL(std::string gradientMode, std::string functionName) {
    std::string headerSource = " -D " + functionName + "(volume,volumeStruct,samplePos)=";
    if (gradientMode == "none")
        headerSource += "(color-(float4)(0.5))*2.0";
    else if (gradientMode == "forward-differences")
        headerSource += "calcGradientAFD(volume,volumeStruct,samplePos)";
    else if (gradientMode == "central-differences")
        headerSource += "calcGradientA(volume,volumeStruct,samplePos)";
    else if (gradientMode == "filtered")
        headerSource += "calcGradientFiltered(volume,volumeStruct,samplePos)";
    else if (gradientMode == "sobel")
        headerSource += "calcGradientSobel(volume,volumeStruct,samplePos)";

    return headerSource;
}

std::string getShaderDefineCL(std::string shadeMode, std::string functionName, std::string n, std::string pos, std::string lPos, std::string cPos, std::string ka, std::string kd, std::string ks) {
    std::string headerSource = " -D " + functionName + "(" + n + "," + pos + "," + lPos + "," + cPos + "," + ka + "," + kd + "," + ks + ",lightSource_" +")=";
    if (shadeMode == "none")
        headerSource += "" + ka + ";";
    else if (shadeMode == "phong-diffuse")
        headerSource += "phongShadingD(" + n + "," + pos + "," + lPos + "," + cPos + "," + kd + ",lightSource_" + ")";
    else if (shadeMode == "phong-specular")
        headerSource += "phongShadingS(" + n + "," + pos + "," + lPos + "," + cPos + "," + ks + ",lightSource_" + ")";
    else if (shadeMode == "phong-diffuse-ambient")
        headerSource += "phongShadingDA(" + n + "," + pos + "," + lPos + "," + cPos + "," + kd + "," + ka + ",lightSource_" + ")";
    else if (shadeMode == "phong-diffuse-specular")
        headerSource += "phongShadingDS(" + n + "," + pos + "," + lPos + "," + cPos + "," + kd + "," + ks + ",lightSource_" + ")";
    else if (shadeMode == "phong")
        headerSource += "phongShading(" + n + "," + pos + "," + lPos + "," + cPos + "," + ka + "," + kd + "," + ks + ",lightSource_"+ ")";
    else if (shadeMode == "toon")
        headerSource += "toonShading(" + n + "," + pos + "," + lPos + "," + cPos + "," + kd + ",3" + ",lightSource_" + ")";
    else if (shadeMode == "cook-torrance")
        headerSource += "cookTorranceShading(" + n + "," + pos + "," + lPos + "," + cPos + "," + ka + "," + kd + "," + ks + ",lightSource_" + ")";
    else if (shadeMode == "oren-nayar")
        headerSource += "orenNayarShading(" + n + "," + pos + "," + lPos + "," + cPos + "," + ka + "," + kd + ",lightSource_" + ")";
    else if (shadeMode == "lafortune")
        headerSource += "lafortuneShading(" + n + "," + pos + "," + lPos + "," + cPos + "," + ka + "," + kd + "," + ks + ",lightSource_" + ")";
    else if (shadeMode == "ward")
        headerSource += "wardShading(" + n + "," + pos + "," + lPos + "," + cPos + "," + ka + "," + kd + "," + ks + ",lightSource_" + ")";

    return headerSource;
}

std::string getCompositingDefineCl(std::string compositingMode, std::string functionName, bool applyDVROpacityCorrection, std::string result, std::string color, std::string samplePos, std::string gradient, std::string t, std::string samplingStepSize, std::string tDepth) {
    std::string headerSource = " -D " + functionName +"(" + result + "," + color + "," +samplePos + "," + gradient + "," + t + "," + samplingStepSize + "," + tDepth + ")=";
    if (compositingMode == "dvr" && applyDVROpacityCorrection)
        headerSource += "compositeDVR(result,color,t,samplingStepSize,tDepth)";
    else if (compositingMode == "dvr")
        headerSource += "compositeDVRNoOpacityCorrection(result,color,t,tDepth)";
    else if (compositingMode == "mip")
        headerSource += "compositeMIP(result,color,t,tDepth)";
    else if (compositingMode == "mida")
        headerSource += "compositeMIDA(result,voxel,color,f_max_i,t,samplingStepSize,tDepth,gammaValue_)";
    else if (compositingMode == "iso")
        headerSource += "compositeISO(result,color,t,tDepth,isoValue_)";
    else if (compositingMode == "fhp")
        headerSource += "compositeFHP(samplePos,result,t,tDepth)";
    else if (compositingMode == "fhn")
        headerSource += "compositeFHN(gradient,result,t,tDepth)";

    return headerSource;
}

std::string getShaderDefineFunctionCL(const std::string mode, const std::string& defineName) {
    std::string headerSource = " -D " + defineName + "(transFunc,transFuncTex,voxel,lastIntensity)=";

    if (mode == "none")
        headerSource += "voxel";
    else if (mode == "transfer-function")
        headerSource += "applyTF(transFunc,transFuncTex,voxel)";
    else if (startsWith(mode, "pre-integrated"))
        headerSource += "applyTFpi(transFunc,transFuncTex,voxel,lastIntensity)";

    return headerSource;
}

} //namespace cl

} //namespace voreen
