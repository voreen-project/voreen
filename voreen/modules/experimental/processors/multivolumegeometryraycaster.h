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

#ifndef VRN_MULTIVOLUMEGEOMETRYRAYCASTER_H
#define VRN_MULTIVOLUMEGEOMETRYRAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/interaction/camerainteractionhandler.h"

#include "voreen/core/ports/volumeport.h"

#ifdef GL_ATOMIC_COUNTER_BUFFER //disable compilation for old gl headers

namespace voreen {

/**
 * OpenGL 4.3 based raycaster supporting multiple volumes using polyhedral proxy-geometries and semi-transparent geometries.
 */
class MultiVolumeGeometryRaycaster : public VolumeRaycaster {
public:
    MultiVolumeGeometryRaycaster();
    ~MultiVolumeGeometryRaycaster();
    Processor* create() const;

    std::string getClassName() const    { return "MultiVolumeGeometryRaycaster"; }
    std::string getCategory() const     { return "Raycasting"; }
    CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }

    /**
     * All inports and at least one outport need to be connected.
     */
    bool isReady() const;

protected:
    struct ProxyGeometry {
        tgt::mat4 textureToWorldMatrix_;
        const Geometry* g_;
        int volumeId_;
        tgt::vec3 cameraPositionPhysical_;
        tgt::vec3 lightPositionPhysical_;
    };

    struct ProxyGeometryGPU {
        tgt::mat4 textureToWorldMatrix_;
        tgt::mat4 worldToTextureMatrix_;
        tgt::vec3 cameraPositionPhysical_;
        int32_t volumeId_;
        tgt::vec3 lightPositionPhysical_;
        int32_t unused1_;
        //int32_t unused2_;
        //int32_t unused3_;
    };
// "The structure’s size will be padded out to a multiple of the size of a vec4." - std140 Uniform Buffer Layout (Appendix L of

struct LinkedListStructGPU {
    uint32_t color_; // 4x8 bit color
    float depth_;
    uint32_t next_;
    uint32_t proxyGeometryId_;
};
// "The structure’s size will be padded out to a multiple of the size of a vec4." - std140 Uniform Buffer Layout (Appendix L of

    virtual void setDescriptions() {
        setDescription("OpenGL 4.3 based raycaster supporting multiple volumes using polyhedral proxy-geometries and semi-transparent geometries.");
    }

    void process();
    /// Compile shader etc.
    void beforeProcess();

    void clearDatastructures();
    /// Returns false if linked list buffer is too small
    void renderGeometries(const std::vector<ProxyGeometry>& pgs);
    void sortFragments();
    void performRaycasting();

    void initialize();
    void deinitialize();

    std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    void compile();
private:
    float getVoxelSamplingStepSize(const VolumeBase* vol, float worldSamplingStepSize);
    void adjustPropertyVisibilities();
    void initializeAtomicBuffer(int numAtomics);
    void initializeStorageBuffer(size_t s);
    void initializeProxyGeometryBuffer(size_t s);
    bool adjustStorageBufferSize();
    void gatherProxyGeometriesRecursive(const Geometry* g, std::vector<ProxyGeometry>& pgs, int volumeId, tgt::mat4 transformation, const tgt::mat4& physicalToTextureMatrix);
    void updateProxyGeometryBuffer(const std::vector<ProxyGeometry>& pgs);

    GeometryPort pgPort1_;
    VolumePort volumeInport2_;
    GeometryPort pgPort2_;
    VolumePort volumeInport3_;
    GeometryPort pgPort3_;
    VolumePort volumeInport4_;
    GeometryPort pgPort4_;

    GeometryPort geometryPort_;


    RenderPort tmpPort_;
    RenderPort headPort_;

    ShaderProperty rcShaderProp_;           ///< The shader performing the actual raycasting/compositing
    ShaderProperty oitShaderProp_;          ///< The shader used for geometry rendering
    ShaderProperty oitClearShaderProp_;     ///< The shader used to clear the datastructures on the GPU
    ShaderProperty oitSortShaderProp_;     ///< The shader used to clear the datastructures on the GPU
    StringOptionProperty sortingAlgorithm_;
    IntProperty numPages_;
    IntProperty maxDepthComplexity_;
    BoolProperty benchmark_;
    BoolProperty adjustSize_;
    IntProperty storageBufferSizeProp_;

    StringOptionProperty classificationMode2_;  ///< What type of transfer function should be used for classification
    StringOptionProperty classificationMode3_;  ///< What type of transfer function should be used for classification
    StringOptionProperty classificationMode4_;  ///< What type of transfer function should be used for classification

    GLuint atomicBuffer_;      ///< GLSL Buffer holding the atomic counter
    int atomicBufferSize_;
    bool atomicBufferInitialized_;
    GLuint storageBuffer_;
    size_t storageBufferSize_;
    size_t targetStorageBufferSize_;
    bool storageBufferInitialized_;
    GLuint proxyGeometryBuffer_;

    StringOptionProperty shadeMode1_;       ///< shading mode to use for volume 1
    StringOptionProperty shadeMode2_;       ///< shading mode to use for volume 2
    StringOptionProperty shadeMode3_;       ///< shading mode to use for volume 3
    StringOptionProperty shadeMode4_;       ///< shading mode to use for volume 4
    FloatProperty clippingGradientDepth_;

    TransFunc1DKeysProperty transferFunc1_;       ///< transfer function to apply to volume 1
    TransFunc1DKeysProperty transferFunc2_;       ///< transfer function to apply to volume 2
    TransFunc1DKeysProperty transferFunc3_;       ///< transfer function to apply to volume 3
    TransFunc1DKeysProperty transferFunc4_;       ///< transfer function to apply to volume 4

    CameraProperty camera_;                 ///< the camera used for lighting calculations
    CameraInteractionHandler* cameraHandler_;

};


} // namespace voreen

#endif //GL_ATOMIC_COUNTER_BUFFER

#endif // VRN_OITRAYCASTER_H
