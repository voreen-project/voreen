/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/gpucapabilities.h"

#include "tgt/logmanager.h"

#ifdef WIN32
    #include <windows.h>
    #include <WinBase.h>
    #include <tchar.h>
    #ifdef _MSC_VER
        #include <strsafe.h>
    #endif
#else // WIN32 -> must be POSIX
    #include <sys/utsname.h> // for uname
#endif // WIN32

#include <string>
#include <sstream>

using std::string;
using std::stringstream;

//-----------------------------------------------------------------------------

namespace { // anonymous helper functions

#ifdef WIN32

/// Converts a STL string to a STL wide-string.
std::wstring str2wstr(const std::string& s);

/// Converts a STL wide-string to a STL string.
std::string wstr2str(const std::wstring& s);

bool parseFileVersionString(tgt::GpuCapabilities::FileVersion &fileVersion);
bool createFileVersionFromDigits(tgt::GpuCapabilities::FileVersion &fileVersion);
bool createVersionStringFromDigits(tgt::GpuCapabilities::FileVersion &fileVersion);

#endif

}

//-----------------------------------------------------------------------------

namespace tgt {

const std::string GpuCapabilities::loggerCat_("tgt.GpuCapabilities");

GpuCapabilities::GpuCapabilities(bool detectCaps) {
    if (detectCaps) {
        detectCapabilities();
        detectOS();
    }

#ifdef TGT_WITH_WMI
    pIWbemLocator_ = 0;
    pWbemServices_ = 0;
    if (!WMIinit()) {
        LWARNING("WMI initialization failed");
    }
#endif

}

GpuCapabilities::~GpuCapabilities() {
#ifdef TGT_WITH_WMI
    if (isWMIinited())
        WMIdeinit();
#endif
}

GpuCapabilities::GlVersion GpuCapabilities::getGlVersion() {
    return glVersion_;
}

bool GpuCapabilities::isOpenGlVersionSupported(GpuCapabilities::GlVersion version) {
    return (glVersion_ >= version);
}

GpuCapabilities::GpuVendor GpuCapabilities::getVendor() {
    return vendor_;
}

std::string GpuCapabilities::getVendorAsString() {
    switch (vendor_) {
    case GPU_VENDOR_NVIDIA:
        return "NVIDIA";
    case GPU_VENDOR_ATI:
        return "ATI";
    case GPU_VENDOR_INTEL:
        return "Intel";
    default:
        return "unknown";
    }
}

bool GpuCapabilities::isExtensionSupported(string extension) const {
    return (glExtensionsString_.find(extension) != string::npos);
}

string GpuCapabilities::getGlVersionString() {
    return glVersionString_;
}

string GpuCapabilities::getGlVendorString() {
    return glVendorString_;
}

string GpuCapabilities::getGlRendererString() {
    return glRendererString_;
}

string GpuCapabilities::getShadingLanguageVersionString() {
    return glslVersionString_;
}

string GpuCapabilities::getGlExtensionsString() {
    return glExtensionsString_;
}

bool GpuCapabilities::areShadersSupported() {
    return shaderSupport_;
}

bool GpuCapabilities::areShadersSupportedARB() {
    return shaderSupportARB_;
}

GpuCapabilities::GlVersion GpuCapabilities::getShaderVersion() {
    return shaderVersion_;
}

GpuCapabilities::ShaderModel GpuCapabilities::getShaderModel() {
    return shaderModel_;
}

bool GpuCapabilities::isShaderModelSupported(GpuCapabilities::ShaderModel shaderModel) {
    return (shaderModel_ >= shaderModel);
}

int GpuCapabilities::getMaxTextureSize() {
    return maxTexSize_;
}

bool GpuCapabilities::is3DTexturingSupported() {
    return texturing3D_;
}

int GpuCapabilities::getMax3DTextureSize() {
    return max3DTexSize_;
}

int GpuCapabilities::retrieveAvailableTextureMemory() const {
    int availableTexMem = -1;

    if (vendor_ == GPU_VENDOR_NVIDIA) {
        if(isExtensionSupported("GL_NVX_gpu_memory_info")){
#ifndef GL_GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX
#define GL_GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX 0x9049
#endif
            GLint retVal;
            glGetIntegerv(GL_GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX, &retVal);
            availableTexMem = static_cast<int>(retVal);
        } else {
            LDEBUG("No GL_NVX_gpu_memory_info support!!!");
        }
    }
    else if (vendor_ == GPU_VENDOR_ATI) {
         /*HACK: the extension GL_ATI_meminfo is  no longer supported by AMD
                However, the enum seems to be still available. Only the GL_ERROR
                has to be removed.*/
        //if(isExtensionSupported("GL_ATI_meminfo")) {
#ifndef GL_TEXTURE_FREE_MEMORY_ATI
#define GL_TEXTURE_FREE_MEMORY_ATI 0x87FC
#endif
            GLint retVal[4];
            glGetIntegerv(GL_TEXTURE_FREE_MEMORY_ATI, retVal);
            availableTexMem = static_cast<int>(retVal[1]); //0=total 1=availabel
        //} else {
        //    LDEBUG("No GL_ATI_meminfo support");
        //}
        //reset error flag to GL_NO_ERROR
        glGetError();
    }

    return availableTexMem;
}

int GpuCapabilities::retrieveTotalTextureMemory() const {
    int availableTexMem = -1;

    if (vendor_ == GPU_VENDOR_NVIDIA) {
        if(isExtensionSupported("GL_NVX_gpu_memory_info")){
#ifndef GL_GPU_MEMORY_INFO_DEDICATED_VIDMEM_NVX
#define GL_GPU_MEMORY_INFO_DEDICATED_VIDMEM_NVX 0x9047
#endif
            GLint retVal;
            glGetIntegerv(GL_GPU_MEMORY_INFO_DEDICATED_VIDMEM_NVX, &retVal);
            availableTexMem = static_cast<int>(retVal);
        } else {
            LDEBUG("No GL_NVX_gpu_memory_info support!!!");
        }
    }
    else if (vendor_ == GPU_VENDOR_ATI) {
        /*HACK: the extension GL_ATI_meminfo is  no longer supported by AMD
                However, the enum seems to be still available. Only the GL_ERROR
                has to be removed.*/
        //if(isExtensionSupported("GL_ATI_meminfo")) {
#ifndef GL_TEXTURE_FREE_MEMORY_ATI
#define GL_TEXTURE_FREE_MEMORY_ATI 0x87FC
#endif
        GLint retVal[4] = {-1,-1,-1,-1};
        glGetIntegerv(GL_TEXTURE_FREE_MEMORY_ATI, retVal);
        availableTexMem = static_cast<int>(retVal[0]); //0=total 1=availabel
        //} else {
        //    LDEBUG("No GL_ATI_meminfo support");
        //}
        //reset error flag to GL_NO_ERROR
        glGetError();
    }

    return availableTexMem;
}

int GpuCapabilities::getNumTextureUnits() {
    return numTextureUnits_;
}

bool GpuCapabilities::isNpotSupported() {
    return npotTextures_;
}

bool GpuCapabilities::areTextureRectanglesSupported() {
    return textureRectangles_;
}

bool GpuCapabilities::isAnisotropicFilteringSupported() {
    return anisotropicFiltering_;
}

float GpuCapabilities::getMaxTextureAnisotropy() {
    return maxTextureAnisotropy_;
}

int GpuCapabilities::getMaxGeometryShaderVertices() {
    return maxGeometryShaderVertices_;
}

bool GpuCapabilities::isTextureCompressionSupported() {
    return textureCompression_;
}

bool GpuCapabilities::areFramebufferObjectsSupported() {
    return framebufferObjects_;
}

int GpuCapabilities::getMaxColorAttachments() {
    return maxColorAttachments_;
}

void GpuCapabilities::logCapabilities(bool extensionsString, bool osString) {
    if (osString)
        LINFO("OS version:          " << osVersionString_);

    LINFO("OpenGL Version:      " << glVersionString_);
    LINFO("OpenGL Renderer:     " << glRendererString_);
    if (glRendererString_.find("Cheetah") != string::npos) {
        LWARNING("It seems that NVIDIA GPU emulation is running on this system.");
        LWARNING("This sometimes leads to strange errors.");
        LWARNING("It is therefore strongly recommended to turn the emulation off!");
    }
    LINFO("GPU Vendor:          " << glVendorString_ << " (" << getVendorAsString() << ")");

    if (extensionsString)
        LINFO("OpenGL Extensions:   " << glExtensionsString_);

    stringstream features;
    features << "Texturing:           " << (isOpenGlVersionSupported(GlVersion::TGT_GL_VERSION_1_1) ? "yes" : "no");
    if (isOpenGlVersionSupported(GlVersion::TGT_GL_VERSION_1_1))
        features << ", max size: " << getMaxTextureSize();
    features << ", 3D: " << (is3DTexturingSupported() ? "yes" : "no");
    if (is3DTexturingSupported())
        features << ", max 3D size: " << getMax3DTextureSize();
    LINFO(features.str());

    features.clear();
    features.str("");
    features << "Texture features:    " << getNumTextureUnits() << " units, "
             << (isNpotSupported() ? "" : "no ") << "NPOT, "
             << (areTextureRectanglesSupported() ? "" : "no") << "rectangles, "
             << (isTextureCompressionSupported() ? "" : "no") << "compression, ";
    if (isAnisotropicFilteringSupported())
        features << getMaxTextureAnisotropy() << "x anisotropic";
    else
        features << "no anisotropic";

    LINFO(features.str());

    features.clear();
    features.str("");
    features << "Framebuffer Objects: " << (areFramebufferObjectsSupported() ? "yes" : "no");
    if(areFramebufferObjectsSupported())
        features << ", max " << getMaxColorAttachments() << " color attachments";
    LINFO(features.str());

    features.clear();
    features.str("");
    features << "Shaders:             " << (areShadersSupported() ? "yes (OpenGL 2.0)" : "no");
    if (!areShadersSupported()) {
        LINFO(features.str());
    }
    else {
        features << ", GLSL Version " << shaderVersion_;
        features << ", Shader Model ";
        if (shaderModel_ == SHADER_MODEL_5)
            features << "5.0";
        else if (shaderModel_ == SHADER_MODEL_4)
            features << "4.0";
        else if (shaderModel_ == SHADER_MODEL_3)
            features << "3.0";
        else if (shaderModel_ == SHADER_MODEL_2)
            features << "2.0";
        else
            features << "unknown";

        if(maxGeometryShaderVertices_ >= 0)
            features << ", max GeometryShader Output Vertices: " << maxGeometryShaderVertices_;

        LINFO(features.str());

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
        if (GLEW_NV_fragment_program2) {
            features.clear();
            features.str("");
            features << "Shader Capabilities: ";

            GLint i = -1;
            glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_LOOP_COUNT_NV, &i);
            features << "MAX_PROGRAM_LOOP_COUNT=" << i;
            glGetProgramivARB(GL_FRAGMENT_PROGRAM_ARB, GL_MAX_PROGRAM_EXEC_INSTRUCTIONS_NV, &i);
            features << ", MAX_PROGRAM_EXEC_INSTRUCTIONS=" << i;
            LDEBUG(features.str());
        }
#endif
    }

    if(retrieveTotalTextureMemory() != -1)
        LINFO("Total Graphics Memory Size:        " << retrieveTotalTextureMemory()/1024 << " MB");
    else
        LINFO("Total Graphics Memory Size:        No Information");

    if(retrieveAvailableTextureMemory() != -1)
        LINFO("Available Graphics Memory Size:    " << retrieveAvailableTextureMemory()/1024 << " MB");
    else
        LINFO("Available Graphics Memory Size:    No Information");

#ifdef TGT_WITH_WMI
    GraphicsDriverInformation driverInfo = getGraphicsDriverInformation();
    LINFO("Graphics Driver Version: " << driverInfo.driverVersion.versionString);
    LINFO("Graphics Driver Date:    " << driverInfo.driverDate);
#endif

}

GpuCapabilities::OSVersion GpuCapabilities::getOSVersion() const{
    return osVersion_;
}

std::string GpuCapabilities::getOSVersionString() {
    return osVersionString_;
}

void GpuCapabilities::detectCapabilities() {
    glVersionString_ = string(reinterpret_cast<const char*>(glGetString(GL_VERSION)));
    glVendorString_  = string(reinterpret_cast<const char*>(glGetString(GL_VENDOR)));
    glRendererString_  = string(reinterpret_cast<const char*>(glGetString(GL_RENDERER)));

    // melt all extentions in one string
    GLint numExtensions;
    glGetIntegerv(GL_NUM_EXTENSIONS, &numExtensions);
    stringstream extensions;
    for (int i = 0; i < numExtensions; ++i) {
        extensions << glGetStringi(GL_EXTENSIONS, i) << std::endl;
    }
    glExtensionsString_ = extensions.str();

    // Prevent segfault
    const char* glslVS = reinterpret_cast<const char*>(glGetString(GL_SHADING_LANGUAGE_VERSION));
    if (glslVS)
        glslVersionString_ = string(glslVS);
    else
        glslVersionString_ = "";

    if (!glVersion_.parseVersionString(glVersionString_)) {
        LERROR("Malformed OpenGL version string: " << glVersionString_);
    }

    // GPU Vendor
    if (glVendorString_.find("NVIDIA") != std::string::npos)
        vendor_ = GPU_VENDOR_NVIDIA;
    else if (glVendorString_.find("ATI") != std::string::npos)
        vendor_ = GPU_VENDOR_ATI;
    else if (glVendorString_.find("INTEL") != std::string::npos ||
             glVendorString_.find("Intel") != std::string::npos)
        vendor_ = GPU_VENDOR_INTEL;
    else {
        LWARNING("Unknown graphics board vendor: " << glVendorString_);
        vendor_ = GPU_VENDOR_UNKNOWN;
    }

    // Shaders
    shaderSupport_ = (glVersion_ >= GlVersion::TGT_GL_VERSION_2_0);
    shaderSupportARB_ = (isExtensionSupported("GL_ARB_vertex_program") &&
                         isExtensionSupported("GL_ARB_fragment_program"));

    if (!shaderVersion_.parseVersionString(glslVersionString_)) {
        LERROR("Malformed GLSL version string: " << glslVersionString_);
    }

    // Shader model
    if (isExtensionSupported("GL_ARB_tessellation_shader"))
        shaderModel_ = SHADER_MODEL_5;
    else if (glVersion_ >= GlVersion::TGT_GL_VERSION_3_2 ||
             isExtensionSupported("GL_ARB_geometry_shader4") ||
             isExtensionSupported("GL_EXT_geometry_shader4"))
        shaderModel_ = SHADER_MODEL_4;
    else if (isExtensionSupported("GL_NV_vertex_program3")  ||
             isExtensionSupported("GL_ATI_shader_texture_lod"))
        shaderModel_ = SHADER_MODEL_3;
    else if (glVersion_ >= GlVersion::TGT_GL_VERSION_2_0)
        shaderModel_ = SHADER_MODEL_2;
    else
        shaderModel_ = SHADER_MODEL_UNKNOWN;

    //Geometry Shader is supported since ShaderModel 4.0
    if(shaderModel_ >= SHADER_MODEL_4)
        glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES,&maxGeometryShaderVertices_);
    else
        maxGeometryShaderVertices_ = -1;

    // Texturing
    glGetIntegerv(GL_MAX_TEXTURE_SIZE, (GLint *) &maxTexSize_);
    texturing3D_ = (glVersion_ >= GlVersion::TGT_GL_VERSION_1_2 ||
                    isExtensionSupported("GL_EXT_texture3D"));

    if (texturing3D_) {
        if (glVersion_ >= GlVersion::TGT_GL_VERSION_2_0)
            glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE, (GLint *) &max3DTexSize_);
        else
            glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE, (GLint *) &max3DTexSize_);
    }
    else {
        max3DTexSize_ = 0;
    }

    numTextureUnits_ = -1;
    if (isExtensionSupported("GL_ARB_fragment_program"))
        glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, (GLint *) &numTextureUnits_);

    if (numTextureUnits_ < 0)
        glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, (GLint *) &numTextureUnits_);

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    npotTextures_ = (isExtensionSupported("GL_ARB_texture_non_power_of_two"));
    textureRectangles_ = (isExtensionSupported("GL_ARB_texture_rectangle"));
    textureCompression_ = (isExtensionSupported("GL_ARB_texture_compression"));
    framebufferObjects_ = (isExtensionSupported("GL_EXT_framebuffer_object"));
#else
    npotTextures_ = true;
    textureRectangles_ = true;
    textureCompression_ = (isExtensionSupported("GL_EXT_texture_compression_s3tc"));
    framebufferObjects_ = true;
#endif

    anisotropicFiltering_ = (isExtensionSupported("GL_EXT_texture_filter_anisotropic"));
    if (anisotropicFiltering_)
        glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &maxTextureAnisotropy_);
    else
        maxTextureAnisotropy_ = 1.0;

    if(framebufferObjects_) {
        glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS, &maxColorAttachments_);
    }
    else
        maxColorAttachments_ = -1;

}

void GpuCapabilities::detectOS() {
    osVersion_ = OS_UNKNOWN;
    osVersionString_ = "unknown";
    std::ostringstream oss;

#ifdef WIN32
    OSVERSIONINFOEX osvi;
    ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);

    if (!GetVersionEx((OSVERSIONINFO*)&osvi)) {
        oss << "unknown: GetVersionEx() failed";
    }
    else if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0) {
        osVersion_ = OS_WIN_2000;
        oss << "Windows 2000";
    }
    else if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1) {
        osVersion_ = OS_WIN_XP;
        oss << "Windows XP";
    }
    else if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 2) {
        osVersion_ = OS_WIN_SERVER_2003;
        oss << "Windows Server 2003";
    }
    else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 0) {
        if (osvi.wProductType == VER_NT_WORKSTATION) {
            osVersion_ = OS_WIN_VISTA;
            oss << "Windows Vista";
        }
        else {
            osVersion_ = OS_WIN_SERVER_2008;
            oss << "Windows Server 2008";
        }
    }
    else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 1) {
        if (osvi.wProductType == VER_NT_WORKSTATION) {
            osVersion_ = OS_WIN_7;
            oss << "Windows 7";
        }
        else {
            osVersion_ = OS_WIN_SERVER_2008_R2;
            oss << "Windows Server 2008 R2";
        }
    }
    else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 2) {
        if (osvi.wProductType == VER_NT_WORKSTATION) {
            osVersion_ = OS_WIN_8;
            oss << "Windows 8";
        }
        else {
            osVersion_ = OS_WIN_SERVER_2012;
            oss << "Windows Server 2012";
        }
    }
    else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 3) {
        if (osvi.wProductType == VER_NT_WORKSTATION) {
            osVersion_ = OS_WIN_8_1;
            oss << "Windows 8.1";
        }
        else {
            osVersion_ = OS_WIN_SERVER_2012_R2;
            oss << "Windows Server 2012 R2";
        }
    }
    else if (osvi.dwMajorVersion == 10 && osvi.dwMinorVersion == 0) {
        if (osvi.wProductType == VER_NT_WORKSTATION) {
            osVersion_ = OS_WIN_10;
            oss << "Windows 10";
        }
        else {
            osVersion_ = OS_WIN_SERVER_2016;
            oss << "Windows Server 2016";
        }
    }
    else {
        oss << "unknown Windows version " << osvi.dwMajorVersion << "." << osvi.dwMinorVersion;
    }
    // Include service pack (if any) and build number.
    if (strlen(osvi.szCSDVersion) > 0) {
        oss << " " << osvi.szCSDVersion;
    }
    oss << " (build " <<  osvi.dwBuildNumber << ")";
    osVersionString_ = oss.str();

#else // WIN32 -> must be UNIX/POSIX
    osVersion_ = OS_POSIX;

    utsname name;
    if (uname(&name) != 0)
        return; // command not successful

    oss << name.sysname << ' ' << name.release << ' ' << name.version << ' ' << name.machine;

#endif // WIN32

    osVersionString_ = oss.str();
}

bool GpuCapabilities::overrideGLSLVersion(const std::string& versionString){
    GlVersion overrideVersion;
    if (overrideVersion.parseVersionString(versionString)) {
        shaderVersion_ = overrideVersion;
        return true;
    }
    else {
        LERROR("GLSL version string '" << versionString << "' could not be parsed. "
               << "Keeping detected version: " << shaderVersion_);
        return false;
    }
}

//-----------------------------------------------------------------------------------


const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_1_0(1,0,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_1_1(1,1,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_1_2(1,2,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_1_3(1,3,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_1_4(1,4,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_1_5(1,5,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_2_0(2,0,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_2_1(2,1,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_3_0(3,0,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_3_1(3,1,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_3_2(3,2,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_3_3(3,3,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_4_0(4,0,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_4_1(4,1,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_4_2(4,2,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_4_3(4,3,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_4_4(4,4,0);
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::TGT_GL_VERSION_4_5(4,5,0);

const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_110(1,10); ///< GLSL version 1.10
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_120(1,20); ///< GLSL version 1.20
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_130(1,30); ///< GLSL version 1.30
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_140(1,40); ///< GLSL version 1.40
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_150(1,50); ///< GLSL version 1.50
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_330(3,30); ///< GLSL version 3.30
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_400(4, 0); ///< GLSL version 4.00
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_410(4,10); ///< GLSL version 4.10
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_420(4,20); ///< GLSL version 4.20
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_430(4,30); ///< GLSL version 4.30
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_440(4,40); ///< GLSL version 4.40
const GpuCapabilities::GlVersion GpuCapabilities::GlVersion::SHADER_VERSION_450(4,50); ///< GLSL version 4.50

GpuCapabilities::GlVersion::GlVersion(int major, int minor, int release)
  : major_(major), minor_(minor), release_(release)
{}

bool GpuCapabilities::GlVersion::parseVersionString(const string& st) {
    major_ = -1;
    minor_ = -1;
    release_ = -1;

    string str;

    //ignore vendor specific part of the string:
    size_t spacePos = st.find_first_of(" ");
    if (spacePos != string::npos)
        str = st.substr(0, spacePos);
    else
        str = st;

    //explode version string with delimiter ".":
    std::vector<string> exploded;
    size_t found;
    found = str.find_first_of(".");
    while(found != string::npos){
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

    stringstream vstr;

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

    //try to parse release number (not always there):
    if (exploded.size() > 2) {
        vstr << exploded[2];
        vstr >> release_;
        if (vstr.fail())
            release_ = 0;
    }
    else
        release_ = 0;

    return true;
}


bool operator==(const GpuCapabilities::GlVersion& x, const GpuCapabilities::GlVersion& y) {
    if ((x.major_ == y.major_) && (x.minor_ == y.minor_) && (x.release_ == y.release_))
        return true;
    else
        return false;
}

bool operator!=(const GpuCapabilities::GlVersion& x, const GpuCapabilities::GlVersion& y) {
    if ((x.major_ != y.major_) || (x.minor_ != y.minor_) || (x.release_ != y.release_))
        return true;
    else
        return false;
}

bool operator<(const GpuCapabilities::GlVersion& x, const GpuCapabilities::GlVersion& y) {
    if (x.major_ < y.major_)
        return true;
    else if (x.major_ == y.major_) {
        if (x.minor_ < y.minor_)
            return true;
        else if (x.minor_ == y.minor_)
            if (x.release_ < y.release_)
                return true;
    }
    return false;
}

bool operator<=(const GpuCapabilities::GlVersion& x, const GpuCapabilities::GlVersion& y) {
    if (x.major_ < y.major_)
        return true;
    else if (x.major_ == y.major_) {
        if (x.minor_ < y.minor_)
            return true;
        else if (x.minor_ == y.minor_)
            if (x.release_ <= y.release_)
                return true;
    }
    return false;
}

bool operator>(const GpuCapabilities::GlVersion& x, const GpuCapabilities::GlVersion& y) {
    if (x.major_ > y.major_)
        return true;
    else if (x.major_ == y.major_) {
        if (x.minor_ > y.minor_)
            return true;
        else if (x.minor_ == y.minor_)
            if (x.release_ > y.release_)
                return true;
    }
    return false;
}

bool operator>=(const GpuCapabilities::GlVersion& x, const GpuCapabilities::GlVersion& y) {
    if (x.major_ > y.major_)
        return true;
    else if (x.major_ == y.major_) {
        if (x.minor_ > y.minor_)
            return true;
        else if (x.minor_ == y.minor_)
            if (x.release_ >= y.release_)
                return true;
    }
    return false;
}

std::ostream& operator<<(std::ostream& s, const GpuCapabilities::GlVersion& v) {
    if (v.major_ == -1)
        return (s << "unknown");
    else if (v.release_ == 0)
        return (s << v.major_ << "." << v.minor_);
    else
        return (s << v.major_ << "." << v.minor_ << "." << v.release_);
}

//-----------------------------------------------------------------------------

#ifdef TGT_WITH_WMI

GpuCapabilities::GraphicsDriverInformation GpuCapabilities::getGraphicsDriverInformation() {
    GraphicsDriverInformation driverInfo;
    driverInfo.driverVersion = getDriverVersion();
    if (driverInfo.driverVersion.versionString.length() == 0) {
        driverInfo.driverVersion.d1 = 0;
        driverInfo.driverVersion.d2 = 0;
        driverInfo.driverVersion.d3 = 0;
        driverInfo.driverVersion.d4 = 0;
        driverInfo.driverVersion.version = 0;
    }
    driverInfo.driverDate = getDriverDate();

    return driverInfo;
}

GpuCapabilities::FileVersion GpuCapabilities::getDriverVersion() {

    FileVersion fileVersion;
    fileVersion.versionString = "";

    // unknown graphics board vendor or reading NVIDIA/Ati driver version failed:
    // get driver version from WMI
    LDEBUG("Reading driver version from WMI ...");

    // Win32_VideoController class: http://msdn2.microsoft.com/en-us/library/aa394512.aspx
    std::string version = WMIqueryStr("Win32_VideoController", "DriverVersion");

    if (version.length() == 0) {
        LDEBUG("Failed to retrieve driver version from WMI.");
        LWARNING("Failed to detect driver version.");
    }
    else {
        fileVersion.versionString = version;
        // do not parse driver version string since it is not necessarily a file version string
        fileVersion.d1 = 0;
        fileVersion.d2 = 0;
        fileVersion.d3 = 0;
        fileVersion.d4 = 0;
        fileVersion.version = 0;
        LDEBUG("Successfully read driver version from WMI");
    }
    return fileVersion;
}

std::string GpuCapabilities::getDriverDate() {

    LDEBUG("Reading driver date from WMI ...");

    if (isWMIinited()) {
        // Win32_VideoController class: http://msdn2.microsoft.com/en-us/library/aa394512.aspx
        std::string date = WMIqueryStr("Win32_VideoController", "DriverDate");

        if (date.length() > 0) {

            LDEBUG("Reading driver date successful.");

            // convert to yyyy-mm-dd format.
            // see http://msdn2.microsoft.com/en-us/library/aa387237.aspx for specification of CIM_DATETIME
            std::string dateformat = "";
            dateformat.append(date.substr(0,4));
            dateformat.append("-");
            dateformat.append(date.substr(4,2));
            dateformat.append("-");
            dateformat.append(date.substr(6,2));

            return dateformat;
        }
        else {
            LDEBUG("Failed reading driver date.");
            LWARNING("Failed to detect driver date");
            return "";
        }
    }
    else {
        LDEBUG("Unable to read DriverDate from WMI: not inited");
    }

    return "";
}

bool GpuCapabilities::WMIinit() {

    // Code based upon:  "Example: Getting WMI Data from the Local Computer"
    // http://msdn2.microsoft.com/en-us/library/aa390423.aspx

    if (isWMIinited()) {
        LWARNING("The WMI connection has already been inited.");
        return false;
    }

    HRESULT hRes;

    // Step 1: --------------------------------------------------
    // Initialize COM. ------------------------------------------

    hRes = CoInitializeEx(0, COINIT_APARTMENTTHREADED);
    if (FAILED(hRes)) {
        LWARNING("CoInitializeEx() failed");
        return false;
    }

    // Step 2: --------------------------------------------------
    // Set general COM security levels --------------------------
    // Note: If you are using Windows 2000, you need to specify -
    // the default authentication credentials for a user by using
    // a SOLE_AUTHENTICATION_LIST structure in the pAuthList ----
    // parameter of CoInitializeSecurity ------------------------

    hRes =  CoInitializeSecurity(
        NULL,
        -1,                          // COM authentication
        NULL,                        // Authentication services
        NULL,                        // Reserved
        RPC_C_AUTHN_LEVEL_DEFAULT,   // Default authentication
        RPC_C_IMP_LEVEL_IMPERSONATE, // Default Impersonation
        NULL,                        // Authentication info
        EOAC_NONE,                   // Additional capabilities
        NULL                         // Reserved
        );

    if (FAILED(hRes)) {

        // Failure of CoInitializeSecurity is not necessarily critical
        // => do not abort initialization

        if (hRes == RPC_E_TOO_LATE) {
            LDEBUG("CoInitializeSecurity failed. HRESULT = RPC_E_TOO_LATE.");
        }
        else if (hRes == RPC_E_NO_GOOD_SECURITY_PACKAGES) {
            LDEBUG("CoInitializeSecurity failed. HRESULT = RPC_E_NO_GOOD_SECURITY_PACKAGES.");
        }
        /*else if (hres == E_OUT_OF_MEMORY) {
            LDEBUG("CoInitializeSecurity failed. HRESULT = E_OUT_OF_MEMORY.");
        } */
        else {
            LDEBUG("CoInitializeSecurity failed. HRESULT = " << hRes);
        }
    }
    else {
        LDEBUG("CoInitializeSecurity successful");

    }

    // Step 3: ---------------------------------------------------
    // Obtain the initial locator to WMI -------------------------

    hRes = CoCreateInstance(
        CLSID_WbemLocator,
        0,
        CLSCTX_INPROC_SERVER,
        IID_IWbemLocator, (LPVOID *) &pIWbemLocator_);

    if (FAILED(hRes)) {
        LWARNING("Failed to create IWbemLocator object. Err code = 0x" << hRes);
        pIWbemLocator_ = 0;
        CoUninitialize();
        return false;     // error is critical
    }

    LDEBUG("IWbemLocator object successfully created");

    // Step 4: -----------------------------------------------------
    // Connect to WMI through the IWbemLocator::ConnectServer method

    // Connect to the root\cimv2 namespace with
    // the current user and obtain pointer pSvc
    // to make IWbemServices calls.
    hRes = pIWbemLocator_->ConnectServer(
         _bstr_t(L"ROOT\\CIMV2"), // Object path of WMI namespace
         NULL,                    // User name. NULL = current user
         NULL,                    // User password. NULL = current
         0,                       // Locale. NULL indicates current
         NULL,                    // Security flags.
         0,                       // Authority (e.g. Kerberos)
         0,                       // Context object
         &pWbemServices_          // pointer to IWbemServices proxy
         );

    if (FAILED(hRes)) {
        LWARNING("Could not connect to WMI server. Error code = 0x" << hRes);
        pIWbemLocator_->Release();
        CoUninitialize();
        pIWbemLocator_ = 0;
        pWbemServices_ = 0;
        return false;            // error is critical
    }

    LDEBUG("Connected to ROOT\\CIMV2 WMI namespace");


    // Step 5: --------------------------------------------------
    // Set security levels on the proxy -------------------------

    hRes = CoSetProxyBlanket(
       pWbemServices_,              // Indicates the proxy to set
       RPC_C_AUTHN_WINNT,           // RPC_C_AUTHN_xxx
       RPC_C_AUTHZ_NONE,            // RPC_C_AUTHZ_xxx
       NULL,                        // Server principal name
       RPC_C_AUTHN_LEVEL_CALL,      // RPC_C_AUTHN_LEVEL_xxx
       RPC_C_IMP_LEVEL_IMPERSONATE, // RPC_C_IMP_LEVEL_xxx
       NULL,                        // client identity
       EOAC_NONE                    // proxy capabilities
    );

    if (FAILED(hRes)) {
        LWARNING("Could not set proxy blanket. Error code = 0x" << hRes);
        pWbemServices_->Release();
        pIWbemLocator_->Release();
        CoUninitialize();
        pIWbemLocator_ = 0;
        pWbemServices_ = 0;
        return false;               // error is critical
    }
    LDEBUG("WMI successfully inited");
    return true;

}

bool GpuCapabilities::WMIdeinit() {

    if (!isWMIinited()) {
        LWARNING("WMI is not inited");
        return false;
    }

    LDEBUG("Deinitializing WMI.");

    if (pIWbemLocator_)
        pIWbemLocator_->Release();
    /*if (pWbemServices_)
        pWbemServices_->Release(); */

    pIWbemLocator_ = 0;
    pWbemServices_ = 0;

    CoUninitialize();

    LDEBUG("Finished Deinitializing WMI.");

    return true;
}

bool GpuCapabilities::isWMIinited() const {
    return (pIWbemLocator_ && pWbemServices_);
}

VARIANT* GpuCapabilities::WMIquery(std::string wmiClass, std::string attribute) const {

    // Code based upon:  "Example: Getting WMI Data from the Local Computer"
    // http://msdn2.microsoft.com/en-us/library/aa390423.aspx

    if (!isWMIinited()) {
        LWARNING("WMI not initiated");
        return 0;
    }

    HRESULT hres;
    VARIANT* result = 0;

    // Step 6: --------------------------------------------------
    // Use the IWbemServices pointer to make requests of WMI ----
    IEnumWbemClassObject* pEnumerator = NULL;
    std::string query = "SELECT " + attribute + " FROM " + wmiClass;
    hres = pWbemServices_->ExecQuery(
        bstr_t("WQL"),
        bstr_t(query.c_str()),
        WBEM_FLAG_FORWARD_ONLY | WBEM_FLAG_RETURN_IMMEDIATELY,
        NULL,
        &pEnumerator);

    if (FAILED(hres)) {
        LWARNING("ERROR: WMI query failed: " << query);
        return 0;
    }

    // Step 7: -------------------------------------------------
    // Get the data from the query in step 6 -------------------

    IWbemClassObject* pclsObj = 0;
    ULONG uReturn = 0;

    if (pEnumerator) {
        HRESULT hr = pEnumerator->Next(WBEM_INFINITE, 1,
            &pclsObj, &uReturn);

        if (uReturn) {
            // Get the value of the attribute and store it in result
            result = new VARIANT;
            hr = pclsObj->Get(LPCWSTR(str2wstr(attribute).c_str()), 0, result, 0, 0);
        }
    }

    if (!result) {
        LWARNING("No WMI query result");
    }

    // Clean enumerator and pclsObject
    if (pEnumerator)
        pEnumerator->Release();
    if (pclsObj)
        pclsObj->Release();

    return result;

}

std::string GpuCapabilities::WMIqueryStr(std::string wmiClass, std::string attribute) {

    VARIANT* variant = WMIquery(wmiClass, attribute);
    if (!variant) {
        return "";
    }
    else {
        std::string result = wstr2str(std::wstring(variant->bstrVal));
        VariantClear(variant);
        return result;
    }
}

int GpuCapabilities::WMIqueryInt(std::string wmiClass, std::string attribute) const {

    VARIANT* variant = WMIquery(wmiClass, attribute);
    if (!variant) {
        return -1;
    }
    else {
        int result = variant->intVal;
        VariantClear(variant);
        return result;
    }
}

#endif // TGT_WITH_WMI

//-----------------------------------------------------------------------------

#ifdef WIN32

GpuCapabilities::FileVersion GpuCapabilities::getFileVersion(const std::string& filename) {
    LDEBUG("Reading file version of file '" << filename << "' ...");

    FileVersion fileVersion;
    fileVersion.versionString = "";

    LPVOID vData = NULL;
    VS_FIXEDFILEINFO *fInfo = NULL;
    DWORD dwHandle;

    DWORD dwSize = GetFileVersionInfoSize(filename.c_str(), &dwHandle);
    if (dwSize <= 0) {
        LDEBUG("Error: Invalid file version info size. File " << filename.c_str() << " not existing?");
        return fileVersion;
    }

    // try to allocate memory
    vData = malloc(dwSize);
    if (vData == NULL) {
        LDEBUG("Error: Allocating memory for file version info failed.");
        SetLastError(ERROR_OUTOFMEMORY);
        return fileVersion;
    }

    // try to the the version-data
    if (GetFileVersionInfo(filename.c_str(), dwHandle, dwSize, vData) == 0) {
        DWORD gle = GetLastError();
        free(vData);
        SetLastError(gle);
        LDEBUG("Error reading file version info.");
        return fileVersion;
    }

    // try to query the root-version-info
    UINT len;
    if (VerQueryValue(vData, TEXT("\\"), (LPVOID*) &fInfo, &len) == 0) {
        DWORD gle = GetLastError();
        free(vData);
        SetLastError(gle);
        LDEBUG("Error reading file version info.");
        return fileVersion;
    }

    // extract hex digits
    int ls = static_cast<int>(fInfo->dwFileVersionLS);
    int ms = static_cast<int>(fInfo->dwFileVersionMS);
    fileVersion.d1 = (ls % 0x10000);
    fileVersion.d2 = (ls / 0xffff);
    fileVersion.d3 = (ms % 0x10000);
    fileVersion.d4 = (ms / 0xffff);

    if ( !createFileVersionFromDigits(fileVersion) ) {
        fileVersion.versionString = "";
        LDEBUG("Error in file version conversion.");
        return fileVersion;
    }

    if ( !createVersionStringFromDigits(fileVersion) ) {
        fileVersion.versionString = "";
        LDEBUG("Error in file version conversion.");
        return fileVersion;
    }

    LDEBUG("Reading file version successful.");
    return fileVersion;
}

std::string GpuCapabilities::getFileDate(const std::string& pFilename) {

#ifdef _MSC_VER
    std::string filename = pFilename;
    LDEBUG("Retrieving file date (last write date) of system file '" << filename << "' ...");

    // add system path to filename
    LPTSTR lpBuffer = static_cast<LPTSTR>(new char[1024]);
    int len = GetSystemDirectory(
        lpBuffer,
        1024
    );
    if (len == 0) {
        LWARNING("Failed to detect system directory.");
        delete lpBuffer;
        return "";
    }
    std::string systemPath = std::string( static_cast<char*>(lpBuffer), len);
    delete lpBuffer;
    systemPath.append("\\");
    systemPath.append(filename);
    filename = systemPath;

    // create file handle for reading attributes
    LDEBUG("Acquiring handle for file '" << filename << "' ...");
    HANDLE filehandle = CreateFile( filename.c_str(),
        FILE_READ_ATTRIBUTES,
        FILE_SHARE_READ,
        0,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        0 );
    if ( filehandle == INVALID_HANDLE_VALUE ) {
        LDEBUG("Failed to acquire handle for file '" << filename << "'.");
        return "";
    }

    // read file information
    LDEBUG("Reading file information ...");
    LPBY_HANDLE_FILE_INFORMATION fileInformation = static_cast<LPBY_HANDLE_FILE_INFORMATION>(new BY_HANDLE_FILE_INFORMATION);
    BOOL success = GetFileInformationByHandle(filehandle, fileInformation);
    if (success == 0) {
        LDEBUG("Failed to read file information of file '" << filename << "'.");
        delete fileInformation;
        return "";
    }

    // convert file time to local time
    SYSTEMTIME stUTC, stLocal;
    DWORD dwRet;
    FileTimeToSystemTime(&(fileInformation->ftLastWriteTime), &stUTC);
    SystemTimeToTzSpecificLocalTime(NULL, &stUTC, &stLocal);

    // build a string showing the date
    LPTSTR lpszString = static_cast<LPTSTR>(new char[1024]);
    dwRet = StringCchPrintf(lpszString, 1024,
        TEXT("%d-%02d-%02d"),
        stLocal.wYear, stLocal.wMonth, stLocal.wDay);
    std::string result = std::string( static_cast<char*>(lpszString));
    delete lpszString;

    LDEBUG("Reading file date successful.");

    return result;

#else
    LWARNING("File date retrieval only supported for MSVC");
    return "";
#endif

}

#endif // WIN32

} // namespace tgt

//-----------------------------------------------------------------------------

namespace { // anonymous helper functions

#ifdef WIN32

std::wstring str2wstr(const std::string& str) {
    int len;
    int slength = (int)str.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, str.c_str(), slength, 0, 0);
    wchar_t* buf = new wchar_t[len];
    MultiByteToWideChar(CP_ACP, 0, str.c_str(), slength, buf, len);
    std::wstring r(buf);
    delete[] buf;
    return r;
}

std::string wstr2str(const std::wstring& wstr) {
    int len;
    int slength = (int)wstr.length() + 1;
    len = WideCharToMultiByte(CP_ACP, 0, wstr.c_str(), slength, 0, 0, 0, 0);
    char* buf = new char[len];
    WideCharToMultiByte(CP_ACP, 0, wstr.c_str(), slength, buf, len, 0, 0);
    std::string r(buf);
    delete[] buf;
    return r;
}

bool parseFileVersionString(tgt::GpuCapabilities::FileVersion &fileVersion) {
    std::string verStr = std::string( fileVersion.versionString );
    std::string substr;
    std::istringstream converter;

    // fourth (= most significant) hex digit
    std::string::size_type pos = verStr.find( ".", 0 );
    if (pos > 0 && pos < 5) {
        substr = verStr.substr( 0, pos );
        converter.clear();
        converter.str(substr);
        converter >> fileVersion.d4;
        if (verStr.length() == pos) {
            LDEBUGC("tgt.GpuCapabilities", "ERROR: Malformed file version string.");
            return false;
        }
        verStr = verStr.substr(pos+1, verStr.length());
    } else {
        LDEBUGC("tgt.GpuCapabilities", "ERROR: Malformed file version string.");
        return false;
    }

    // third hex digit
    pos = verStr.find( ".", 0 );
    if (pos > 0 && pos < 5) {
        substr = verStr.substr( 0, pos );
        converter.clear();
        converter.str(substr);
        converter >> fileVersion.d3;
        if (verStr.length() == pos) {
            LWARNINGC("tgt.GpuCapabilities", "Malformed file version string.");
            return false;
        }
        verStr = verStr.substr(pos+1, verStr.length());
    } else {
        LDEBUGC("tgt.GpuCapabilities", "ERROR: Malformed file version string.");
        return false;
    }

    // second hex digit
    pos = verStr.find( ".", 0 );
    if (pos > 0 && pos < 5) {
        substr = verStr.substr( 0, pos );
        converter.clear();
        converter.str(substr);
        converter >> fileVersion.d2;
        if (verStr.length() == pos) {
            LDEBUGC("tgt.GpuCapabilities", "ERROR: Malformed file version string.");
            return false;
        }
        verStr = verStr.substr(pos+1, verStr.length());
    } else {
        LDEBUGC("tgt.GpuCapabilities", "ERROR: Malformed file version string.");
        return false;
    }

    // first hex digit
    if ( verStr.length() > 0 && verStr.length() < 5) {
        converter.clear();
        converter.str(verStr);
        converter >> fileVersion.d1;
    } else {
        LDEBUGC("tgt.GpuCapabilities", "ERROR: Malformed file version string.");
        return false;
    }

    return createFileVersionFromDigits(fileVersion);
}

bool createFileVersionFromDigits(tgt::GpuCapabilities::FileVersion &fileVersion) {
    fileVersion.version = ( static_cast<uint64_t>(fileVersion.d4) << 48 )
        + ( static_cast<uint64_t>(fileVersion.d3) << 32 )
        + ( static_cast<uint64_t>(fileVersion.d2) << 16 )
        + ( static_cast<uint64_t>(fileVersion.d1) );

    return true;
}

bool createVersionStringFromDigits(tgt::GpuCapabilities::FileVersion &fileVersion) {
    std::ostringstream converter;

    fileVersion.versionString = "";

    converter.clear();
    converter << fileVersion.d4;
    fileVersion.versionString.append(converter.str());

    fileVersion.versionString.append(".");
    converter.str("");
    converter.clear();
    converter << fileVersion.d3;
    fileVersion.versionString.append(converter.str());

    fileVersion.versionString.append(".");
    converter.str("");
    converter.clear();
    converter << fileVersion.d2;
    fileVersion.versionString.append(converter.str());

    fileVersion.versionString.append(".");
    converter.str("");
    converter.clear();
    converter << fileVersion.d1;
    fileVersion.versionString.append(converter.str());

    return true;
}

#endif // WIN32

} // namespace anonymous
