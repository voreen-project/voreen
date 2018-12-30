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

#ifndef VRN_RAYCASTERCL_H
#define VRN_RAYCASTERCL_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "modules/opencl/utils/clwrapper.h"
#include "modules/opencl/processors/openclprocessor.h"
#include "modules/opencl/properties/openclproperty.h"

namespace voreen {

class VRN_CORE_API RaycasterCL : public cl::OpenCLProcessor<VolumeRenderer> {
public:
    RaycasterCL();
    virtual ~RaycasterCL();

    virtual std::string getClassName() const    { return "RaycasterCL"; }
    virtual std::string getCategory() const     { return "Raycasting"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isReady() const;

    virtual Processor* create() const;
    void portResized();

    virtual void initializeCL();
    virtual void deinitializeCL();

protected:
    virtual void setDescriptions() {
        setDescription("A simple raycaster using OpenCL");
    }

    virtual void beforeProcess();
    virtual void process();

    virtual void initialize();
    virtual void deinitialize();

    void updateTF();

    // ports
    VolumePort volumePort_;
    RenderPort entryPort_;
    RenderPort exitPort_;

    RenderPort outport_;

    /// Category used for logging.
    static const std::string loggerCat_;

private:

    FloatProperty samplingRate_;  ///< Sampling rate of the raycasting, specified relative to the size of one voxel
    TransFunc1DKeysProperty transferFunc_;
    OpenCLProperty openclProp_;

    //VolumeRAM_4xUInt8* volumeData_; //Currently not used

    //float samplingStepSize_; //Currently not used

    cl::OpenCL* opencl_;
    cl::Context* context_;
    cl::CommandQueue* queue_;

    cl::ImageObject3D* volumeTex_;
    cl::ImageObject2D* tfData_;
    cl::SharedTexture* entryTexCol_;
    cl::SharedTexture* exitTexCol_;
    cl::SharedTexture* outCol_;
};

} // namespace voreen

#endif
