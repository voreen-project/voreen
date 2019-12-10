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

#ifndef VRN_CANNY_H
#define VRN_CANNY_H

#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

class Canny : public ImageProcessor {
public:
    Canny();
    virtual Processor* create() const;

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "Canny"; }
    virtual CodeState getCodeState() const { return CODE_STATE_OBSOLETE; }

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    void process();

    FloatProperty startThreshold_;
    FloatProperty runThreshold_;
    ColorProperty edgeColor_; ///< The color used for drawing the edges

    RenderPort inport_;
    RenderPort outport_;
    RenderPort privatePort_;
};


} // namespace voreen

#endif //VRN_CANNY_H
