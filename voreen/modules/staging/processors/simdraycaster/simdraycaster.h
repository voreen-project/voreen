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

#ifndef VRN_SIMDRAYCASTER_H
#define VRN_SIMDRAYCASTER_H

#include "../../utils/simdraycaster/brickedvolumebase.h"
#include "../../utils/simdraycaster/performancemetric.h"

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/fontproperty.h"


#include <cstdint>
#include <memory>

namespace voreen {
class CameraInteractionHandler;
    
/**
 *
 */
class SIMDRayCaster : public VolumeRaycaster {
public:

    SIMDRayCaster();
    virtual ~SIMDRayCaster();

    /**
     * Allocation of resources. Called by voreen.
     */
    virtual void initialize();

    /**
     * Cleanup of resources. Called by voreen.
     */
    virtual void deinitialize();

    /**
     * The main method, that does the rendering. It is called by voreen each frame.
     */
    virtual void process();

    // Configuration for voreen
    virtual Processor* create() const { return new SIMDRayCaster(); }
    virtual std::string getCategory() const { return "Raycaster"; }
    virtual std::string getClassName() const { return "SIMDRayCaster"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }
protected:
    /**
     * Gets the transferfunction
     * \param size[out] The size of the transfer function, that is returned. 
     *                  This is written by the function
     * \return The transferfunction. This array has \p size elements
     */
    tgt::vec4* getTransFunc1D(size_t *size);

    /**
     * Gets the preintegrated transferfunction
     * \param[out] size The size of the transfer function, that is returned. 
     *                  This is written by the function.
     * \param samplingrate[in] The samplingrate used for preintegration.
     *                         This needs to equal the samplingrate used
     *                         for raycasting.
     * \return The transferfunction. This array has \p size * \p size elements
     */
    tgt::vec4* getPreintegrateTransFunc(size_t *size, float samplingrate);
    
    
    /**
     * Performs opacitiy correction and alpha premultiplication  to a transfer function.
     * \param[in,out] data The data that is used. This is modified!
     *                     It must use the default stepsize!
     * \param elements[in] The number of elements in the \p data  array
     * \param volumeDimension[in] The dimension of the volume used for the 
     *                            calculation of the default stepsize in voreen.
     * \param stepsize[in] the desired stepsize.
     */
    void precalcOpacityCorrectionAndAlphaPremultiplication(tgt::vec4* data, size_t elements, tgt::svec3 volumeDimension, float stepsize);

    /**
     * Performs alpha premultiplication.
     * \param data[in,out] The array of values modified.
     * \param elements[in] The number of elements
     *        in the \p data array
     */
    void premultiplyAlpha(tgt::vec4* data, size_t elements);

    /**
     * Outputs a image from a buffer on the cpu to the gpu framebuffer.
     * \param image[in] The array of RGBA colors
     * \param size[in] Size of the Framebuffer.
     */
    void blitImage(std::uint32_t *img, tgt::ivec2 size);

    /**
     * Creates a bricked volume
     * \param[in] VolumeRAM The volume data, that should be used.
     * \param[in] bricksize Size of the Bricks for the bricked volume
     *                      Each brick has this value cubed voxel-
     * \return The bricked volume.
     */
    BrickedVolumeBase * createBrickedVolume(const VolumeRAM * vol, int bricksize);

    /**
     * Adapts the processor to a new volume.
     */
    void adaptToNewVolume();

private:
    /**
     * Sets descriptions in the GUI
     */ 
    virtual void setDescriptions() {
        setDescription("Renders a Volume with SSE instructions on the CPU");
        performanceInfoPort_.setDescription("Outputs performance info as Text.");
        transferFunc_.setDescription("Transferfunction used for rendering.");
        preintegration_.setDescription("Enables the use of preintegration of the tranfer function for DVR.");
        mipRendering_.setDescription("Selects the mode of rendering.");
        pixelBlockConfiguration_.setDescription("Sets the block size of pixels processed together.");
        enableBrickedVolume_.setDescription("Enables the use of the bricked volume as a optimization for some datasets.");
        brickSize_.setDescription("The size of a block for the bricked volume.");
    }

    /**
     * Executed after process for performance metrics
     */
    virtual void afterProcess();

    /**
     * Tells voreen if the processor is ready
     */
    virtual bool isReady() const;

    TextPort performanceInfoPort_; ///< Outport for info on performance metrics
    CameraProperty camera_; ///< The camera of the scene.
    TransFunc1DKeysProperty transferFunc_;  ///< the property that controls the transfer function
    BoolProperty  preintegration_; ///< use preintegration for DVR rendering
    IntOptionProperty  mipRendering_; ///< use MIP rendering
    IntVec2Property pixelBlockConfiguration_; ///< configuration of blocks of pixels
    BoolProperty  enableBrickedVolume_; ///< use the bricked volume
    IntProperty brickSize_; ///< Blocksize for the bricked Volume

    GLuint resultTexure_; ///< Texture for the resulting images, reused for performance reasons

    BrickedVolumeBase* brickedVolume_; ///< The data for the bricked volume
    PerformanceMetric performanceMetric_; ///< Data for the performance metrics
};

}   // namespace

#endif  // VRN_SIMDRAYCASTER_H

