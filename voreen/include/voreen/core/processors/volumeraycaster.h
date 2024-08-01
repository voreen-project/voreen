/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMERAYCASTER_H
#define VRN_VOLUMERAYCASTER_H

#include "tgt/types.h"
#include "tgt/shadermanager.h"
#include "tgt/textureunit.h"
#include "tgt/exception.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/datastructures/volume/volumegl.h"
#include "voreen/core/properties/collectivesettingsproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/utils/voreenqualitymode.h"

namespace voreen {

/**
 * Abstract base class for GPU-based ray-casters.
 *
 * It extends the generateHeader, setGlobalShaderParameters and
 * bindVolumes methods inherited from VolumeRenderer.
 * Additionally, it provides several ray-casting-related properties
 * and handles bricked volumes.
 */
class VRN_CORE_API VolumeRaycaster : public VolumeRenderer, public QualityModeObserver {
public:
    VolumeRaycaster(bool addDefaultProps = true, bool usePortgroup = false, bool addPorts = true);

    /// Switches interaction coarseness on/off by resizing all renderports.
    virtual void qualityModeChanged();

    /**
     * Resizes the RenderPorts, if interaction coarseness is currently active.
     *
     * @see RenderProcessor::invalidate
     */
    virtual void invalidate(int inv = INVALID_RESULT);

protected:

    void addCompositionModes(StringOptionProperty &compositionProp);

    virtual void setDescriptions() override {
        VolumeRenderer::setDescriptions();
        //default properties
        renderingQualityOption_.setDescription("Groups all rendering quality settings, e.g., the sampling rate. The value <i>full</i> \
                                                equals to the best quality but may be slow on certain hardware.");
        interactionQualityOption_.setDescription("Groups all rendering quality settings during interaction, e.g., camera movements. The value <i>full</i> \
                                                  equals to the best quality but may be slow on certain hardware.");
        samplingRate_.setDescription("Defines the step size during raycasting. The default value is 2.0 which equals a step size of half a voxel. Higher values \
                                      improve the quality, but may be slow on certain hardware.");
        interactionCoarseness_.setDescription("The screen (output) resolution will be devided by this value during interaction. An increase of frames per seconds \
                                               will be gernerated through the decreased number of pixels.");
        interactionQuality_.setDescription("The sampling rate will be multiplied by this value during interation.");
        //non default properties
        compositingMode1_.setDescription("Which compositing should be used during raycasting? The default value is (D)irect (V)olume (R)endering.");
        shadeMode_.setDescription("Which shading/light model should be applied?");
        gradientMode_.setDescription("Which type of calculation should be used for on-the-fly gradients required for the shading.");
    }

    virtual void initialize();
    virtual void deinitialize();

    virtual void afterProcess();

    /**
     * Defines ray-casting macros to be used in the raycasting shader.
     *
     * @see VolumeRenderer::generateHeader
     */
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    std::string generateCompositionHeader(std::string macroname, std::string compositionMode, bool usePreIntegrationClassification, bool isActive);

    /**
     * Sets frustum parameters necessary for depth value calculation in shaders.
     * The camera parameter is passed to the super class function.
     *
     * @see VolumeRenderer::setGlobalShaderParameters
     */
    virtual void setGlobalShaderParameters(tgt::Shader* shader, const tgt::Camera* camera = 0, tgt::ivec2 screenDim = tgt::ivec2(-1));

    /**
     * Binds volume textures (inherited from VolumeRenderer) and sets the sampling step size
     * relative to the resolution of the first volume. The camera and light position
     * parameters are passed to the super class function.
     *
     * @see VolumeRenderer::bindVolumes
     */
    virtual bool bindVolumes(tgt::Shader* shader, const std::vector<VolumeStruct> &volumes,
        const tgt::Camera* camera = 0, const tgt::vec4& lightPosition = tgt::vec4(0.f));

    /**
     * Copies over the content of the srcPort to the destPort,
     * thereby implicitly rescaling the image to the dest dimensions.
     * To be used by subclasses for implementing coarseness (i.e., rendering with reduced dimensions in interaction mode).
     */
    void rescaleRendering(RenderPort& srcPort, RenderPort& destPort);

    /// Calculate sampling step size for a given volume using the current sampling rate
    float getSamplingStepSize(const VolumeBase* vol);

    VolumePort volumeInport_;
    RenderPort entryPort_;
    RenderPort exitPort_;

    RenderPort outport1_;
    RenderPort outport2_;
    RenderPort outport3_;

    // we render into internal buffers, which allows to reduce rendering size in interaction mode (coarseness)
    RenderPort internalRenderPort1_;
    RenderPort internalRenderPort2_;
    RenderPort internalRenderPort3_;
    PortGroup internalPortGroup_;

    CollectiveSettingsProperty renderingQualityOption_;   ///< generell settings used to adjust all rendering qualities
    CollectiveSettingsProperty interactionQualityOption_; ///< generell settings used to adjust all interaction qualities

    FloatProperty samplingRate_;  ///< Sampling rate of the raycasting, specified relative to the size of one voxel
    FloatProperty isoValue_;      ///< The used isovalue, when isosurface raycasting is enabled
    FloatProperty gammaValue1_;
    FloatProperty gammaValue2_;
    FloatProperty gammaValue3_;

    // properties for all volume raycasters
    StringOptionProperty maskingMode_;                 ///< What masking should be applied (thresholding, segmentation)
    StringOptionProperty gradientMode_;                ///< What type of calculation should be used for on-the-fly gradients
    StringOptionProperty classificationMode_;          ///< What type of transfer function should be used for classification
    StringOptionProperty shadeMode_;                   ///< What shading method should be applied
    StringOptionProperty compositingMode1_;             ///< What compositing mode should be applied
    StringOptionProperty compositingMode2_;             ///< What compositing mode should be applied
    StringOptionProperty compositingMode3_;             ///< What compositing mode should be applied

    IntProperty interactionCoarseness_;                ///< RenderPorts are resized to size_/interactionCoarseness_ in interactionmode
    FloatProperty interactionQuality_;                 ///< Used to reduce the sampling step size


    tgt::Shader* rescaleShader_;                       ///< Shader used by the rescaleRendering() method
    bool usePortgroup_;

    static const std::string loggerCat_; ///< category used in logging
};

} // namespace voreen

#endif // VRN_VOLUMERAYCASTER_H
