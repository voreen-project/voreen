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

#ifndef VRN_CONVOLUTION_H
#define VRN_CONVOLUTION_H

#include "voreen/core/processors/imageprocessorbypassable.h"

#include "voreen/core/properties/intproperty.h"

namespace voreen {

class VRN_CORE_API Convolution : public ImageProcessorBypassable {
public:
    Convolution();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "Convolution"; }
    virtual std::string getCategory() const  { return "Image Processing"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Applies an arbitrary filter mask to an image. The mask is passed through the second RenderPort.");
    }

    void process();

    IntProperty filterSize_;

    RenderPort inport_;
    RenderPort filterPort_;
    RenderPort outport_;
};

} // namespace

#endif // VRN_CONVOLUTION_H
