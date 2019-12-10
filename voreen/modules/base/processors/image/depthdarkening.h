/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_DEPTHDARKENING_H
#define VRN_DEPTHDARKENING_H

#include "voreen/core/ports/renderport.h"
#include "voreen/core/processors/imageprocessorbypassable.h"
#include "voreen/core/properties/floatproperty.h"

namespace voreen {

/**
 * Performs depth darkening as
 * presented by Luft et al. in 2006.
 */
class VRN_CORE_API DepthDarkening : public ImageProcessorBypassable {
public:
    DepthDarkening();
    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "DepthDarkening"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }
    virtual Processor* create() const;
    void process();

protected:
    virtual void setDescriptions() {
        setDescription("Performs depth darkening as presented by Luft et al. at Siggraph 2006. The technique is based on an unsharping mask applied to the depth buffer. The kernel size of the used Gaussian filter is set to round(2.5*sigma). Lambda controls the strength of the effect.");
    }

    FloatProperty sigma_; ///< Standard deviation of the used Gaussian filter
    FloatProperty lambda_; ///< Controls the strength of the depth masking effect

    RenderPort inport_;
    RenderPort outport_;
    RenderPort privatePort_;
};


} // namespace voreen

#endif //VRN_DEPTHDARKENING_H
