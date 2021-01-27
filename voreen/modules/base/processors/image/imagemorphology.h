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

#ifndef VRN_IMAGEMORPHOLOGY_H
#define VRN_IMAGEMORPHOLOGY_H

#include "voreen/core/processors/imageprocessorbypassable.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

/**
 * Performs a dilation or erosion operation with variable kernel size.
 */
class VRN_CORE_API ImageMorphology : public ImageProcessorBypassable{
public:
    ImageMorphology();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Image Processing"; }
    virtual std::string getClassName() const { return "ImageMorphology"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Performs erosion or dilation operation with variable Kernel Size");
        modeProp_.setDescription("Sets the filter to either Erosion or Dilation mode");
        kernelRadius_.setDescription("Change radius of the kernel that is used for filtering");
    }

    virtual void process();

    RenderPort inport_;
    RenderPort outport_;

    IntProperty kernelRadius_;
    StringOptionProperty modeProp_;
    StringOptionProperty shapeProp_;
};


} // namespace voreen

#endif //VRN_IMAGEMORPHOLOGY_H
