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

#ifndef VRN_SINGLEVOLUMESLICER_H
#define VRN_SINGLEVOLUMESLICER_H

#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/processors/volumeslicer.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * Slicing-based volume renderer using view plane aligned slices.
 */
class VRN_CORE_API SingleVolumeSlicer : public VolumeSlicer {
public:
    SingleVolumeSlicer();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "SingleVolumeSlicer"; }
    virtual std::string getCategory() const  { return "Slice Rendering";    }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;    }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("This class implements a slicing based volume renderer. "
            "The slices are view plane aligned and clipped against the volume on the GPU as described "
            "in 'Real-Time Volume Graphics' (1st Ed., Chapter 3.5.1). "
            "At the moment no correct depth values are calculated. Currently only 1D transfer functions are supported.");
    }

    virtual void process();
    virtual void initialize();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void compile();

private:

    RenderPort outport_;

    ShaderProperty shaderProp_;

    // interaction handlers
    CameraInteractionHandler cameraHandler_;
};


} // namespace voreen

#endif // VRN_SINGLEVOLUMESLICER_H
