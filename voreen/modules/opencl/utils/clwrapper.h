/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_CLWRAPPER_H
#define VRN_CLWRAPPER_H

#include "tgt/tgt_gl.h"

#ifndef __APPLE__
#include <CL/cl.h>
#include <CL/cl_gl.h>
#else
#include <OpenCL/cl.h>
#include <OpenCL/cl_gl.h>
#endif

#include "voreen/core/properties/optionproperty.h"
#include <string>
#include <vector>
#include <set>
#include <map>
#include "tgt/vector.h"
#include "tgt/texture.h"
#include "voreen/core/datastructures/volume/volumeram.h"

namespace voreen {

namespace cl {

std::string clErrorToString(cl_int err);
cl_int _lCLError(cl_int err, int line, const char* file);
#define LCL_ERROR(err) _lCLError(err, __LINE__, __FILE__)

class Device;
class Platform;
class Kernel;
class Buffer;
class SharedTexture;
class ImageObject2D;
class ImageObject3D;
class MemoryObject;
class CommandQueue;

// Nvidia OpenCL SDK versions prior to 3.0 declared functions with __cdecl.
// Newer versions use __stdcall. This will adapt our code to the one used.
void CL_API_CALL contextCallback(const char *errinfo, const void* private_info, size_t cb, void* user_data);

//---------------------------------------------------------

class VRN_CORE_API OpenCL {
public:
    OpenCL();
    const std::vector<Platform>& getPlatforms() const { return platforms_; }

    /**
     * Selects the platform that matches the passed vendor. If no matching
     * platform is found, the first of the available platforms is returned.
     */
    Platform getPlatformByVendor(tgt::GpuCapabilities::GpuVendor vendor) const;

    ///Needs CL extension, doesn´t link
    static Device getCurrentDeviceForGlContext();
    ///Needs CL extension, doesn´t link
    static std::vector<Device> getDevicesForGlContext();
    /// Computes an appropriate global work size as the multiple
    /// of local workSize that is nearest (greater) to numItems
    static size_t shrRoundUp(size_t localWorkSize, size_t numItems);

protected:
    std::vector<Platform> platforms_;

    static const std::string loggerCat_;
};

//---------------------------------------------------------

class VRN_CORE_API Platform {
public:
    enum Profile {
        FULL_PROFILE = 0,
        EMBEDDED_PROFILE,
        UNKNOWN
    };

    /**
     * Specifies the version of the OpenCL implementation.
     * CL_VERSION_x_y denotes OpenCL version x.y.
     */
    class ClVersion {
        public:
        ClVersion(int major = 0, int minor = 0);

        /**
         * Parse OpenCL version string as specified:
         *
         * The CL_PLATFORM_VERSION strings begins with "OpenCL " and is followed by a version number.
         * The version number uses the following format:
         *      major_version.minor_version
         *
         * Vendor-specific information may follow the version number.
         * Its format depends on the implementation, but a space always separates the version number and the vendor-specific information.
         */
        bool parseVersionString(const std::string& st);

        int major_;
        int minor_;

        // major and minor are in brackets to avoid confusion with (although now deprecated) glibc macro definitions
        int (major)() const { return major_; }
        int (minor)() const { return minor_; }

        friend bool operator==(const ClVersion& x, const ClVersion& y);
        friend bool operator!=(const ClVersion& x, const ClVersion& y);
        friend bool operator<(const ClVersion& x, const ClVersion& y);
        friend bool operator<=(const ClVersion& x, const ClVersion& y);
        friend bool operator>(const ClVersion& x, const ClVersion& y);
        friend bool operator>=(const ClVersion& x, const ClVersion& y);
        friend std::ostream& operator<<(std::ostream& s, const ClVersion& v);

        static const ClVersion VRN_CL_VERSION_1_0;
    };

    Platform();
    Platform(cl_platform_id id);
    ~Platform();

    cl_platform_id getID() const { return id_; }
    const std::vector<Device>& getDevices() const { return devices_; }
    Profile getProfile() const { return profile_; }
    std::string getName() const { return name_; }
    std::string getVendor() const { return vendor_; }
    std::string getVersionString() const { return versionString_; }
    std::string getExtensionString() const { return extensionString_; }
    const std::set<std::string>& getExtensions() const { return extensions_; }
    bool isExtensionSupported(std::string ext) const;
    bool supportsGlSharing() const { return isExtensionSupported("cl_khr_gl_sharing"); }

    template<class T>
    T getInfo(cl_platform_info info) const {
        T ret;
        LCL_ERROR(clGetDeviceInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    void logInfos() const;

protected:

    cl_platform_id id_;

    //info:
    Profile profile_;
    std::string name_;
    std::string vendor_;
    std::string extensionString_;
    std::set<std::string> extensions_;
    std::string versionString_;
    ClVersion version_;

    std::vector<Device> devices_;

    static const std::string loggerCat_;
};

//template specialization for strings:
template<>
std::string Platform::getInfo(cl_platform_info info) const;

//---------------------------------------------------------

class VRN_CORE_API Device {
public:
    Device();
    Device(cl_device_id id);

    cl_device_id getId() const { return id_; }

    template<class T>
    T getInfo(cl_device_info info) const {
        T ret;
        LCL_ERROR(clGetDeviceInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    void logInfos() const;

    std::string getName() const { return name_; }
    cl_device_type getType() const { return getInfo<cl_device_type>(CL_DEVICE_TYPE); }
    std::string getExtensionString() const { return extensionString_; }
    const std::set<std::string>& getExtensions() const { return extensions_; }
    bool isExtensionSupported(std::string ext) const;
    bool supportsGlSharing() const { return isExtensionSupported("cl_khr_gl_sharing"); }
    cl_uint getMaxWorkGroupSize() const;

    bool hasImageSupport() const;
    tgt::ivec2 getMaxImageSize2D() const;
    tgt::ivec3 getMaxImageSize3D() const;

protected:
    cl_device_id id_;

    std::string name_;
    std::string extensionString_;
    std::set<std::string> extensions_;

    bool imageSupport_;
    tgt::ivec2 maxImageSize2D_;
    tgt::ivec3 maxImageSize3D_;

    static const std::string loggerCat_;
};

//template specialization for strings:
template<>
std::string Device::getInfo(cl_device_info info) const;

//---------------------------------------------------------

struct VRN_CORE_API ContextProperty {
    cl_context_properties name_;
    cl_context_properties value_;
    ContextProperty(cl_context_properties name, cl_context_properties value) : name_(name), value_(value) {}
};

class VRN_CORE_API Context {
public:
    Context(const Device& device);
    Context(const std::vector<ContextProperty>& properties, const Device& device);
    //TODO:
    //Context(const std::vector<Device>& devices);
    //Context(const Platform& platform, const std::vector<Device>& devices);
    //Context(const Platform& platform, const Device& device);
    //TODO: create context from type
    //TODO: additional context properties
    //TODO: notification functions, user data
    ~Context();

    ///Tested for Linux, code for windows is untested, code for apple is commented out
    static std::vector<ContextProperty> generateGlSharingProperties();

    cl_context getId() const { return id_; }

    template<class T>
    T getInfo(cl_context_info info) const {
        T ret;
        LCL_ERROR(clGetContextInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    //TODO: getDevices()
    //TODO: getContextProperties()
protected:
    cl_context id_;

    static const std::string loggerCat_;
};

//template specialization for strings:
template<>
std::string Context::getInfo(cl_context_info info) const;

//---------------------------------------------------------

class VRN_CORE_API Event {
    friend class CommandQueue;

public:
    Event() : id_(0) {
    }

    Event(cl_event id) : id_(id) {
        //if(id_)
            //LCL_ERROR(clRetainEvent(id_));
    }

    Event(const Event& e) {
        id_ = e.getId();
        //if(id_)
            //LCL_ERROR(clRetainEvent(id_));
    }

    ~Event() {
        //if(id_)
            //LCL_ERROR(clReleaseEvent(id_));
    }

    void wait() const { LCL_ERROR(clWaitForEvents(1, &id_)); }

    void releaseEvent() {
        if (id_)
            clReleaseEvent(id_);
    }

    template<class T>
    T getInfo(cl_event_info info) const {
        T ret;
        LCL_ERROR(clGetEventInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    ///Possible values for info: CL_PROFILING_COMMAND_QUEUED CL_PROFILING_COMMAND_SUBMIT CL_PROFILING_COMMAND_START CL_PROFILING_COMMAND_END
    cl_ulong getProfilingInfo(cl_profiling_info info) const {
        cl_ulong ret;
        LCL_ERROR(clGetEventProfilingInfo(id_, info, sizeof(cl_ulong), &ret, 0));
        return ret;
    }

    cl_ulong getProfilingStartToEnd() const { return (getProfilingInfo(CL_PROFILING_COMMAND_END) - getProfilingInfo(CL_PROFILING_COMMAND_START)); }

    /**
     * Return the execution status of the command identified by event. The valid values are:
     * CL_QUEUED (command has been enqueued in the command-queue),
     * CL_SUBMITTED (enqueued command has been submitted by the host to the device associated with the command-queue),
     * CL_RUNNING (device is currently executing this command),
     * CL_COMPLETE (the command has completed), or Error code given by a negative integer value.
     */
    cl_int getExecutionStatus() {
        return getInfo<cl_int>(CL_EVENT_COMMAND_EXECUTION_STATUS);
    }

    /**
     * Return the command associated with event. Can be one of the following values:
     * CL_COMMAND_NDRANGE_KERNEL
     * CL_COMMAND_TASK
     * CL_COMMAND_NATIVE_KERNEL
     * CL_COMMAND_READ_BUFFER
     * CL_COMMAND_WRITE_BUFFER
     * CL_COMMAND_COPY_BUFFER
     * CL_COMMAND_READ_IMAGE
     * CL_COMMAND_WRITE_IMAGE
     * CL_COMMAND_COPY_IMAGE
     * CL_COMMAND_COPY_BUFFER_TO_IMAGE
     * CL_COMMAND_COPY_IMAGE_TO_BUFFER
     * CL_COMMAND_MAP_BUFFER
     * CL_COMMAND_MAP_IMAGE
     * CL_COMMAND_UNMAP_MEM_OBJECT
     * CL_COMMAND_MARKER
     * CL_COMMAND_ACQUIRE_GL_OBJECTS
     * CL_COMMAND_RELEASE_GL_OBJECTS
     */
    cl_command_type getCommandType() {
        return getInfo<cl_command_type>(CL_EVENT_COMMAND_TYPE);
    }

    cl_event getId() const { return id_; }
protected:
    cl_event id_;
};

//---------------------------------------------------------

class VRN_CORE_API CommandQueue {
public:
#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    CommandQueue(const Context* context, const Device& device, cl_command_queue_properties properties = 0);
#else
    CommandQueue(const Context* context, const Device& device, const cl_queue_properties* properties = 0);
#endif
    ~CommandQueue();

    //TODO: set command queue properties
#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    cl_command_queue_properties getProperties() const { return properties_; }
#else
    const cl_queue_properties* getProperties() const { return properties_; }
#endif
    cl_int flush() { return LCL_ERROR(clFlush(id_)); }
    cl_int finish() { return LCL_ERROR(clFinish(id_)); }

    void enqueue(const Kernel* kernel, size_t globalWorkSize, size_t localWorkSizes, Event* event = 0);
    void enqueue(const Kernel* kernel, size_t globalWorkSize, Event* event = 0);
    void enqueue(const Kernel* kernel, tgt::svec2 globalWorkSizes, tgt::svec2 localWorkSizes, Event* event = 0);
    void enqueue(const Kernel* kernel, tgt::svec2 globalWorkSizes, Event* event = 0);
    void enqueue(const Kernel* kernel, tgt::svec3 globalWorkSizes, tgt::svec3 localWorkSizes, Event* event = 0);
    void enqueue(const Kernel* kernel, tgt::svec3 globalWorkSizes, Event* event = 0);
    void enqueue(const Kernel* kernel, const std::vector<size_t>& globalWorkSizes, const std::vector<size_t>& localWorkSizes, Event* event = 0);
    //TODO: add method with offsets
    //TODO: wrap clEnqueueNativeKernel

    void enqueueReadBuffer(const Buffer* buffer, void* data, bool blocking = true, Event* event = 0);
    void enqueueReadBuffer(const Buffer* buffer, size_t byteOffset, size_t numBytes, void* data, bool blocking = true, Event* event = 0);

    void enqueueWriteBuffer(const Buffer* buffer, void* data, bool blocking = true, Event* event = 0);
    void enqueueWriteBuffer(const Buffer* buffer, size_t byteOffset, size_t numBytes, void* data, bool blocking = true, Event* event = 0);

    //FIND: k_nguyen add enqueueCopyImage
    void enqueueCopyBuffer(const Buffer &src, const Buffer &dst, size_t src_offset, size_t dst_offset, size_t size, Event* event = 0);

    void enqueueCopyImage(const ImageObject2D &src,
                           const ImageObject2D &dst,
                           const size_t src_origin[3],
                           const size_t dst_origin[3],
                           const size_t region[3],
                           Event* event = 0);
    void enqueueWriteImage(const ImageObject2D &img,
                            tgt::Texture *tex,
                            const size_t origin[3],
                            const size_t region[3],
                            bool blocking=true,
                            Event* event = 0);
    void enqueueReadImage(const ImageObject2D &img,
                           tgt::Texture *tex,
                           const size_t origin[3],
                           const size_t region[3],
                           bool blocking=true,
                           Event* event = 0);
    void enqueueReadImage(const ImageObject2D &img,
                           tgt::Texture *tex,
                           bool blocking=true,
                           Event* event = 0);

    void enqueueCopyImageToBuffer(const ImageObject2D &src,
                                   const Buffer &dst,
                                   const size_t origin[3],
                                   const size_t region[3],
                                   size_t   dst_offset = 0,
                                   Event* event = 0);

    void enqueueCopyBufferToImage(const Buffer &src,
                                   const ImageObject2D &dst,
                                   size_t src_offset,
                                   const size_t dst_origin[3],
                                   const size_t region[3],
                                   Event* event = 0);

    //

    //TODO: add methods with offset + size
    //TODO: wrap clEnqueueCopyBuffer

    void enqueueAcquireGLObject(const MemoryObject* obj, Event* event = 0);
    void enqueueReleaseGLObject(const MemoryObject* obj, Event* event = 0);

#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    void enqueueMarker(Event* event = 0);
    void enqueueBarrier();
    void enqueueWaitForEvent(Event* event);
    void enqueueWaitForEvents(const std::vector<Event*>& event);
#else
    void enqueueMarker(const std::vector<Event*>& waitlist = std::vector<Event*>(), Event* event = 0);
    void enqueueBarrier(const std::vector<Event*>& waitlist = std::vector<Event*>(), Event* event = 0);
#endif

    template<class T>
    T getInfo(cl_command_queue_info info) const {
        T ret;
        LCL_ERROR(clGetCommandQueueInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    cl_command_queue getId() const { return id_; }

protected:
    cl_command_queue id_;
#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    cl_command_queue_properties properties_;
#else
    const cl_queue_properties* properties_;
#endif

    static const std::string loggerCat_;
};

//template specialization for strings:
template<>
std::string CommandQueue::getInfo(cl_command_queue_info info) const;

//---------------------------------------------------------

class VRN_CORE_API Program {
public:
    Program(const Context& context);
    Program(const Context* context);
    ~Program();

    void setSource(const std::string& source);
    void setSource(const std::vector<std::string>& source);
    bool loadSource(const std::string& filename);
    bool loadSource(const std::vector<std::string>& filenames);

    std::string getHeader() const;
    ///Does NOT automatically recompile the program.
    void setHeader(const std::string& header);

    ///Reload source from disk (if program has been loaded from disk)
    bool reload();

    std::string getBuildOptions() const { return buildOptions_; }
    ///This does NOT automatically recompile the program.
    void setBuildOptions(const std::string& buildOptions) { buildOptions_ = buildOptions; }

    bool build();
    bool build(const Device& device);
    bool build(const std::vector<Device>& devices);

    /**
     * Get a Kernel.
     * \warning Do not save this pointer to use it later. Kernels get invalid and are deleted if the Program is rebuild.
     *
     * @param name Name of the Kernel function
     *
     * @return The Kernel named 'name' or 0 if no such Kernel exists or the Program is invalid or not build.
     */
    Kernel* getKernel(const std::string& name);

    std::string getBuildLog(const Device& device) const;
    ///Returns one of the following: CL_BUILD_NONE CL_BUILD_ERROR CL_BUILD_SUCCESS CL_BUILD_IN_PROGRESS
    cl_build_status getBuildStatus(const Device& device) const;
    std::string getBuildOptions(const Device& device) const;

    const std::map<std::string, Kernel*>& getCurrentKernels() const {
        return kernels_;
    }

    template<class T>
    T getInfo(cl_program_info info) const {
        T ret;
        LCL_ERROR(clGetProgramdInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    //TODO: get binaries
    template<class T>
    T getBuildInfo(const Device& device, cl_program_build_info info) const {
        T ret;
        LCL_ERROR(clGetProgramBuildInfo(id_, device.getId(), info, sizeof(ret), &ret, 0));
        return ret;
    }

    static std::string statusToString(cl_build_status status);

    cl_program getId() const { return id_; }
protected:

    void createWithSource();
    void createKernels();
    void clearKernels();
    void releaseProgram();

    cl_program id_;
    const Context* context_;

    std::string buildOptions_;
    std::vector<std::string> source_;
    std::string header_;
    std::vector<std::string> filenames_;
    std::map<std::string, Kernel*> kernels_;
    static const std::string loggerCat_;
};

//template specialization for strings:
template<>
std::string Program::getBuildInfo(const Device& device, cl_program_build_info info) const;

template<>
std::string Program::getInfo(cl_program_info info) const;

//---------------------------------------------------------

class VRN_CORE_API MemoryObject {
public:
    MemoryObject();
    virtual ~MemoryObject();

    cl_mem getId() const { return id_; }
protected:
    cl_mem id_;
};

//---------------------------------------------------------

class VRN_CORE_API Buffer : public MemoryObject {
public:
    ///flags: CL_MEM_READ_WRITE CL_MEM_WRITE_ONLY CL_MEM_READ_ONLY CL_MEM_USE_HOST_PTR CL_MEM_ALLOC_HOST_PTR CL_MEM_COPY_HOST_PTR
    //Buffer(const Context& context, cl_mem_flags flags, size_t size, void* hostPtr = 0);
    //TODO: const hostPtr
    Buffer(const Context* context, cl_mem_flags flags, size_t size, const void* hostPtr = 0, cl_int* err = 0);
    virtual ~Buffer();

    size_t getSize() const { return size_; }
protected:
    size_t size_;
};

//---------------------------------------------------------

class VRN_CORE_API SharedTexture : public MemoryObject {
public:
    ///flags: CL_MEM_READ_WRITE CL_MEM_WRITE_ONLY CL_MEM_READ_ONLY CL_MEM_USE_HOST_PTR CL_MEM_ALLOC_HOST_PTR CL_MEM_COPY_HOST_PTR
    SharedTexture(const Context* context, cl_mem_flags flags, tgt::Texture* tex);
    virtual ~SharedTexture();

    size_t getSize() const { return size_; }

protected:
    size_t size_;
    tgt::Texture* tex_;
};


//---------------------------------------------------------
struct VRN_CORE_API ImageFormat : public cl_image_format {
    ImageFormat(cl_channel_order order, cl_channel_type type);
};


//---------------------------------------------------------
class VRN_CORE_API ImageObject2D : public MemoryObject {
public:
    ///flags: CL_MEM_READ_WRITE CL_MEM_WRITE_ONLY CL_MEM_READ_ONLY CL_MEM_USE_HOST_PTR CL_MEM_ALLOC_HOST_PTR CL_MEM_COPY_HOST_PTR
    ImageObject2D(const Context *context, cl_mem_flags flags, tgt::Texture *tex);
    ImageObject2D(const Context *context, cl_mem_flags flags,
                  ImageFormat format,
                  size_t width,
                  size_t height,
                  size_t row_pitch=0,
                  void *host_ptr=NULL);

    virtual ~ImageObject2D();
};

//---------------------------------------------------------

class VRN_CORE_API ImageObject3D : public MemoryObject {
public:
    ///flags: CL_MEM_READ_WRITE CL_MEM_WRITE_ONLY CL_MEM_READ_ONLY CL_MEM_USE_HOST_PTR CL_MEM_ALLOC_HOST_PTR CL_MEM_COPY_HOST_PTR
    ImageObject3D(const Context* context, cl_mem_flags flags, const VolumeRAM* vol);
    virtual ~ImageObject3D();
};

//---------------------------------------------------------

class VRN_CORE_API Sampler {
public:
    
#ifdef CL_USE_DEPRECATED_OPENCL_1_2_APIS
    /**
     * @param context Must be a valid OpenCL context.
     * @param normalized_coords Determines if the image coordinates specified are normalized (if normalized_coords is CL_TRUE) or not (if normalized_coords is CL_FALSE).
     * @param addressing_mode Specifies how out-of-range image coordinates are handled when reading from an image. This can be set to CL_ADDRESS_REPEAT, CL_ADDRESS_CLAMP_TO_EDGE, CL_ADDRESS_CLAMP, and CL_ADDRESS_NONE.
     * @param filtering_mode Specifies the type of filter that must be applied when reading an image. This can be CL_FILTER_NEAREST or CL_FILTER_LINEAR.
     */
    Sampler(const Context* context, cl_bool normalized_coords, cl_addressing_mode addressing_mode, cl_filter_mode filter_mode);
#else
    /**
    * @param context Must be a valid OpenCL context.
    * @param properties Properties of the sampler. See documentation for appropriate values.
    */
    Sampler(const Context* context, const cl_sampler_properties* properties);
#endif
    ~Sampler();

    cl_sampler getId() const { return id_; }

    template<class T>
    T getInfo(cl_sampler_info info) const {
        T ret;
        LCL_ERROR(clGetSamplerInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    ///Return the value specified by addressing_mode argument to clCreateSampler.
    cl_addressing_mode getAdressingMode() { return getInfo<cl_addressing_mode>(CL_SAMPLER_ADDRESSING_MODE); }

    ///Return the value specified by addressing_mode argument to clCreateSampler.
    cl_filter_mode getFilterMode() { return getInfo<cl_filter_mode>(CL_SAMPLER_FILTER_MODE); }

    ///Return the value specified by filter_mode argument to clCreateSampler.
    cl_bool    getNormalizedCoords() { return getInfo<cl_bool>(CL_SAMPLER_NORMALIZED_COORDS); }

protected:
    cl_sampler id_;
};

//---------------------------------------------------------

class VRN_CORE_API Kernel {
public:
    Kernel(cl_kernel id);
    ~Kernel();

    ///Be sure this does what you want!
    template<class T>
    cl_int setArgTemplate(cl_uint index, const T& data) {
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(data), &data));
    }

#ifdef __APPLE__
    cl_int setArg(cl_uint index, size_t data) { return setArgTemplate(index, data); }
#endif
    cl_int setArg(cl_uint index, float data) { return setArgTemplate(index, data); }

    cl_int setArg(cl_uint index, int data) { return setArgTemplate(index, data); }

    cl_int setArg(cl_uint index, uint32_t data) { return setArgTemplate(index, data); }

    cl_int setArg(cl_uint index, uint64_t data) { return setArgTemplate(index, data); }

    cl_int setArg(cl_uint index, tgt::vec2 data) {
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(float)*2, data.elem));
       }

    cl_int setArg(cl_uint index, tgt::vec3 data) {
        return LCL_ERROR(setArg(index, tgt::vec4(data, 0.0f) ));
       }

    cl_int setArg(cl_uint index, tgt::vec4 data) {
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(float)*4, data.elem));
       }

    cl_int setArg(cl_uint index, tgt::ivec2 data) {
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(int)*2, data.elem));
       }

    cl_int setArg(cl_uint index, tgt::ivec3 data) {
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(int)*3, data.elem));
       }

    cl_int setArg(cl_uint index, tgt::ivec4 data) {
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(int)*4, data.elem));
       }

    cl_int setArg(cl_uint index, tgt::mat4 data) {
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(float)*16, data.elem));
       }

    cl_int setArg(cl_uint index, const MemoryObject* memObj) {
        cl_mem id = memObj->getId();
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(cl_mem), &id));
    }

    cl_int setArg(cl_uint index, const MemoryObject& memObj) {
        cl_mem id = memObj.getId();
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(cl_mem), &id));
    }

    cl_int setArg(cl_uint index, const Sampler* sampler) {
        cl_sampler id = sampler->getId();
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(cl_sampler), &id));
    }

    cl_int setArg(cl_uint index, const Sampler& sampler) {
        cl_sampler id = sampler.getId();
        return LCL_ERROR(clSetKernelArg(id_, index, sizeof(cl_sampler), &id));
    }

    cl_int setArg(cl_uint index, size_t size, const void* data) {
        return LCL_ERROR(clSetKernelArg(id_, index, size, data));
    }

    template<class T>
    T getInfo(cl_kernel_info info) const {
        T ret;
        LCL_ERROR(clGetKernelInfo(id_, info, sizeof(ret), &ret, 0));
        return ret;
    }

    template<class T>
    T getWorkGroupInfo(cl_kernel_work_group_info info, const Device& device) const{
        T ret;
        LCL_ERROR(clGetKernelWorkGroupInfo(id_, device.getId(), info, sizeof(ret), &ret, 0));
        return ret;
    }

    cl_kernel getId() const { return id_; }
    cl_uint getNumArgs() const;
    cl_uint getWorkGroupSize(const Device& device) const;
    std::string getName() const;
protected:
    cl_kernel id_;

    static const std::string loggerCat_;
};

//template specialization for strings:
template<>
std::string Kernel::getInfo(cl_kernel_info info) const;

template<>
std::string Kernel::getWorkGroupInfo(cl_kernel_work_group_info info, const Device& device) const;

//-----------------------------------------------------------------------------------------------

struct VRN_CORE_API VolumeWriteBufferCL {
    tgt::Vector4<cl_uint> dimensions_;
    cl_uchar numChannels_;
    cl_uchar numBitsPerChannel_;
    cl_uchar isSigned_;
    cl_uchar isFloat_;
    intptr_t data_;
};

struct VRN_CORE_API VolumeWriteBuffer {

    VolumeWriteBuffer(const Context* context, VolumeRAM* vol);

    VolumeWriteBuffer(const VolumeWriteBuffer& buf) : infoBuffer_(buf.infoBuffer_), dataBuffer_(buf.dataBuffer_), refCount_(buf.refCount_), volume_(buf.volume_) {
        refCount_++;
    }

    ~VolumeWriteBuffer() {

        if(refCount_ == 0) {
            delete infoBuffer_;
            delete dataBuffer_;
            infoBuffer_ = 0;
            dataBuffer_ = 0;
        } else
            refCount_--;
    }

    void readBackData(CommandQueue* q) {
        q->enqueueReadBuffer(dataBuffer_, volume_->getData(), true);
    }

    Buffer* infoBuffer_;
    Buffer* dataBuffer_;
    unsigned short refCount_;

    VolumeRAM* volume_;
};

VolumeWriteBuffer createVolumeWriteBuffer(const Context* context, VolumeRAM* vol);

//Some useful typedefs:

typedef tgt::vec4 float4;
typedef tgt::vec2 float2;

void fillGradientModesPropertyCL(StringOptionProperty& gradientMode);
void fillShadingModesPropertyCL(StringOptionProperty& shadeMode);
void fillCompositingModesPropertyCL(StringOptionProperty& compositingMode);
void fillClassificationModesPropertyCL(StringOptionProperty* prop);
std::string getGradientDefineCL(std::string gradientMode, std::string functionName);
std::string getShaderDefineCL(std::string shadeMode, std::string functionName, std::string n = "n", std::string pos = "pos", std::string lPos = "lPos", std::string cPos = "cPos", std::string ka = "ka", std::string kd = "kd", std::string ks = "ks");
std::string getCompositingDefineCl(std::string compositingMode, std::string functionName, bool applyDVROpacityCorrection = true, std::string result = "result", std::string color = "color", std::string samplePos = "samplePos", std::string gradient = "gradient", std::string t = "t", std::string samplingStepSize = "samplingStepSize", std::string tDepth = "tDepth");
std::string getShaderDefineFunctionCL(const std::string mode, const std::string& defineName);
} //namespace cl

} //namespace voreen

#endif // VRN_CLWRAPPER_H
