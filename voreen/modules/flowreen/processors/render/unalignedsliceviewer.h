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

#ifndef VRN_UNALIGNEDSLICEVIEWER_H
#define VRN_UNALIGNEDSLICEVIEWER_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * Performs slice rendering of a single slice on an arbitrary plane.
 */
class VRN_CORE_API UnalignedSliceViewer : public VolumeRenderer {

public:
    UnalignedSliceViewer();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "UnalignedSliceViewer";  }
    virtual std::string getCategory() const     { return "Slice Rendering";       }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }

protected:

    virtual void setDescriptions() {
        setDescription("Performs slice rendering of a single slice on an arbitrary plane.");
    }

    virtual void initialize();
    virtual void deinitialize();

    virtual void beforeProcess();
    virtual void process();
    virtual void afterProcess();

    virtual void adjustPropertiesToInput();

    /// Generates the header for the shader depending on the choice of features to be used.
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    /// Recompiles the shader.
    bool rebuildShader();

private:

    VolumePort inport_;
    RenderPort outport_;

    FloatVec3Property normal_;
    FloatProperty distance_;
    FloatProperty samplingRate_;
    TransFunc1DKeysProperty transferFunc_;

    tgt::Shader* sliceShader_;

    static const std::string loggerCat_;

};

} // namespace

#endif // VRN_UNALIGNEDSLICEVIEWER_H
