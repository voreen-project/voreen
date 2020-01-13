/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_EEPGEOMETRYINTEGRATOR_H
#define VRN_EEPGEOMETRYINTEGRATOR_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"

namespace voreen {

class VRN_CORE_API EEPGeometryIntegrator : public ImageProcessor {

public:
    EEPGeometryIntegrator();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "EEPGeometryIntegrator"; }
    virtual std::string getCategory() const  { return "Entry-Exit Points";     }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;       }
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Modifies raycasting entry/exit points with regard to the input geometry rendering such that alpha compositing of the geometry rendering with the raycasted image gives correct results (for opaque geometry). The raycasted volume has to be passed to the volume inport.\
<p>See Compositor</p>");
    }

    virtual void beforeProcess();
    virtual void process();

    RenderPort inport0_;
    RenderPort inport1_;
    RenderPort geometryPort_;
    VolumePort volumeInport_;

    RenderPort entryPort_;
    RenderPort exitPort_;
    RenderPort tmpPort_;

    BoolProperty useFloatRenderTargets_;
    CameraProperty camera_;

private:
    /// category used in logging
    static const std::string loggerCat_;
};


} // namespace voreen

#endif
