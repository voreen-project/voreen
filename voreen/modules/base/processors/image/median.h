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

#ifndef VRN_MEDIAN_H
#define VRN_MEDIAN_H

#include "voreen/core/processors/imageprocessorbypassable.h"

#include "voreen/core/properties/intproperty.h"

namespace voreen {

/**
 * Performs a median filtering.
 */
class VRN_CORE_API Median : public ImageProcessorBypassable {
public:
    Median();
    virtual Processor* create() const {return new Median();}

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "Median"; }
    virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Median image filter with selectable kernel size.");
    }

    void process();

    IntProperty halfKernelDim_; ///< Half the kernel size of the box filter

    RenderPort inport_;
    RenderPort outport_;
};


} // namespace voreen

#endif //VRN_MEDIAN_H
