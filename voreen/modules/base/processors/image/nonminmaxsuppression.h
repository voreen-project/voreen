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

#ifndef VRN_NONMINMAXSUPPRESSION_H
#define VRN_NONMINMAXSUPPRESSION_H

#include "voreen/core/processors/imageprocessorbypassable.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

/**
 * Performs a non-maximum suppression.
 */
class VRN_CORE_API NonMinMaxSuppression : public ImageProcessorBypassable {
public:
    NonMinMaxSuppression();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Image Processing"; }
    virtual std::string getClassName() const { return "NonMinMaxSuppression"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Performs a gray-value based non-maximum or non-minimum suppression, according to the selected mode. The kernel radius can be specified for both dimensions independently. For anisotropic non-min/max suppression, an additional sobel image has to be provided.");
    }

    void process();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void compile();

    RenderPort inportImage_;
    RenderPort inportSobel_;
    RenderPort outport_;

    StringOptionProperty mode_;
    BoolProperty isotropic_;
    IntVec2Property kernelRadius_;
};


} // namespace voreen

#endif //VRN_NONMAXSUPPRESSION_H
