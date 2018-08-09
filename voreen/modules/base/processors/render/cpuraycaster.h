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

#ifndef VRN_CPURAYCASTER_H
#define VRN_CPURAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/properties/transfunc/transfunctypeproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"

namespace voreen {

/**
 * This is a simple CPURaycaster.
 * The processor allows the use of pre-integration for 1D transfer functions.
 * OpenMP is used for multithreading, if the OpenMP module is activated.
 */
class VRN_CORE_API CPURaycaster : public VolumeRaycaster {
public:
    CPURaycaster();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "CPURaycaster"; }
    virtual std::string getCategory() const   { return "Raycasting";   }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }

    /// All ports except the gradient port need to be connected.
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Performs a simple ray casting on the CPU. Supports OpenMP for parallization (if OpenMP Module is activated) and Pre-Integration for 1D transfer functions.");
    }

    virtual void process();

    virtual void initialize();

private:

    enum CpuRaycasterClassificationMode {
        TF,
        PREINTEGRATED
    };

    /**
     * Performs the actual ray casting for a single ray,
     * which determined by the passed entry and exit points.
     */
    virtual tgt::vec4 directRendering(const tgt::vec3& first, const tgt::vec3& last, tgt::Texture* tfTexture, const VolumeRAM* volume, float samplingStepSize, const PreIntegrationTable* table, CpuRaycasterClassificationMode mode, tgt::mat4 textureToVoxelMatrix);
    tgt::vec4 apply1DTF(tgt::Texture* tfTexture, float intensity);
    tgt::vec4 apply2DTF(tgt::Texture* tfTexture, float intensity, float gradientMagnitude);

    VolumePort gradientVolumePort_;

    TransFuncTypeProperty tFuncType_;                   ///< the property that controls the transfer-function type
    TransFunc1DKeysProperty tFunc1DKeys_;               ///< the property that controls the 1d keys transfer-function
    TransFunc2DPrimitivesProperty tFunc2DPrimitives_;   ///< the property that controls the 2d primitives transfer-function
    IntOptionProperty texFilterMode_;                   ///< texture filtering mode to use for volume access

    IntOptionProperty preIntegrationTableSize_; ///< sets the width of the Pre-Integration table
};

} // namespace voreen

#endif // VRN_CPURAYCASTER_H
