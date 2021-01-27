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

#ifndef VRN_SINGLEVOLUMERAYCASTER_H
#define VRN_SINGLEVOLUMERAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"

#include "voreen/core/properties/transfunc/transfunctypeproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/1d/1dgaussian/transfunc1dgaussianproperty.h"
#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * This is the standard raycaster within Voreen. It allows to generate three output
 * renderings, whereas only the first one provides depth values. There are several
 * shading and compositing modes available.
 *
 * @see CubeMeshProxyGeometry, MeshEntryExitPoints
 * @see MultiVolumeRaycaster
 */
class VRN_CORE_API SingleVolumeRaycaster : public VolumeRaycaster {
public:
    SingleVolumeRaycaster();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "SingleVolumeRaycaster"; }
    virtual std::string getCategory() const     { return "Raycasting"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    /// All inports and at least one outport need to be connected
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("This is the standard volume renderer in Voreen. It generates up to three output renderings, where the first one also provides depth values. Several shading and compositing modes are supported.\
<p>See CubeProxyGeometry, MeshEntryExitPoints.</p>");
    }

    /**
     * Recompiles the shader, if the invalidation level >= Processor::INVALID_PROGRAM.
     */
    virtual void beforeProcess();

    virtual void process();

    /**
     * Loads the shader and initializes the port group.
     */
    virtual void initialize();

    /**
     * Deinitializes the port group and disposes the shader.
     */
    virtual void deinitialize();

    /**
     * Adds compositing macros for the additional outports.
     *
     * @see VolumeRaycaster::generateHeader
     */
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    /// Rebuilds the loaded shader.
    virtual void compile();

private:
    /**
     * Handels all processor changes related to volume changes.
     */
    void volumeInportOnChange();
    /**
     * Adjusts the property visibility according to the current property settings.
     */
    void adjustPropertyVisibilities();

    TransFuncTypeProperty transFuncTypeProp_;       ///< property to switch between transfer functions
    TransFunc1DKeysProperty transFunc1DProp1_;       ///< the property that controls the 1d transfer function
    TransFunc1DKeysProperty transFunc1DProp2_;       ///< the property that controls the 1d transfer function
    TransFunc1DKeysProperty transFunc1DProp3_;       ///< the property that controls the 1d transfer function
    TransFunc1DKeysProperty transFunc1DProp4_;       ///< the property that controls the 1d transfer function
    TransFunc2DPrimitivesProperty transFunc2DProp1_; ///< the property that controls the 2d transfer function
    TransFunc2DPrimitivesProperty transFunc2DProp2_; ///< the property that controls the 2d transfer function
    TransFunc2DPrimitivesProperty transFunc2DProp3_; ///< the property that controls the 2d transfer function
    TransFunc2DPrimitivesProperty transFunc2DProp4_; ///< the property that controls the 2d transfer function

    TransFunc1DGaussianProperty transFunc1DGaussianProp1_;       ///< the property that controls the 1d gaussian transfer function
    TransFunc1DGaussianProperty transFunc1DGaussianProp2_;
    TransFunc1DGaussianProperty transFunc1DGaussianProp3_;
    TransFunc1DGaussianProperty transFunc1DGaussianProp4_;

    FloatVec3Property channelShift1_;
    FloatVec3Property channelShift2_;
    FloatVec3Property channelShift3_;
    FloatVec3Property channelShift4_;
    BoolProperty enableChannelShift_;

    ShaderProperty shaderProp_;                     ///< the shader property used by this raycaster.
    CameraProperty camera_;                         ///< the camera used for lighting calculations

    static const std::string loggerCat_; ///< category used in logging
};


} // namespace voreen

#endif // VRN_SINGLEVOLUMERAYCASTER_H
