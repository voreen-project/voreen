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

#ifndef VRN_EXPLOSIONCOMPOSITOR_H
#define VRN_EXPLOSIONCOMPOSITOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"

namespace voreen {

/**
 * Composites the partial renderings of an exploded view.
 *
 * @see ExplosionProxyGeometry
 */
class VRN_CORE_API ExplosionCompositor : public ImageProcessor {

public:
    ExplosionCompositor();
    ~ExplosionCompositor();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "ExplosionCompositor"; }
    virtual std::string getCategory() const     { return "Image Processing";    }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE;     }

protected:
    virtual void setDescriptions() {
        setDescription("Composites the partial renderings of an exploded view.\
<p>See ExplosionProxyGeometry</p>");
    }

    virtual void process();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void compile();

    void compositingModeChanged();

    StringOptionProperty compositingMode_;

    FloatProperty weightingFactor_;

    RenderPort inport0_;
    RenderPort interalPort_;
    RenderPort outport_;
    LoopPort loopOutport_;
};


} // namespace voreen

#endif //VRN_EXPLOSIONCOMPOSITOR_H
