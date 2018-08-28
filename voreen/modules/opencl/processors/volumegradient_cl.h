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

#ifndef VRN_VOLUMEGRADIENT_CL_H
#define VRN_VOLUMEGRADIENT_CL_H

#include "voreen/core/processors/volumeprocessor.h"
#include "modules/opencl/utils/clwrapper.h"
#include "modules/opencl/processors/openclprocessor.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include <string>

namespace voreen {

class Volume;

class VRN_CORE_API VolumeGradientCL : public cl::OpenCLProcessor<VolumeProcessor> {
public:
    VolumeGradientCL();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual std::string getClassName() const  { return "VolumeGradientCL";  }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

    virtual void initializeCL();
    virtual void deinitializeCL();

protected:
    virtual void setDescriptions() {
        setDescription("OpenCL-based version of the VolumeGradient processor.");
    }

    virtual void process();
    virtual void beforeProcess();

private:
    StringOptionProperty technique_;
    StringOptionProperty outputType_;
    BoolProperty copyIntensityChannel_;  ///< if true, a four channel volume is created including the original intensity channel

    VolumePort inport_;
    VolumePort outport_;

    cl::OpenCL* opencl_;
    cl::Context* context_;
    cl::CommandQueue* queue_;
    cl::Program* prog_;

    cl::ImageObject3D* volumeTex_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
