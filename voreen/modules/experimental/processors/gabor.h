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

#ifndef VRN_GABOR_H
#define VRN_GABOR_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class Gabor : public ImageProcessor {
public:
    Gabor();
    virtual Processor* create() const;

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "Gabor"; }
    virtual CodeState getCodeState() const { return CODE_STATE_TESTING; }

    virtual void portResized(RenderPort* /*p*/, tgt::ivec2 /*newsize*/) {}

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    void process();

    FloatProperty angle_;
    FloatProperty wavelength_;
    FloatProperty offset_;
    FloatProperty sigma_;
    FloatProperty aspectRatio_;
    IntProperty resolution_;

    RenderPort outport_;
};

} // namespace

#endif // VRN_GABOR_H
