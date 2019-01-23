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

#ifndef VRN_IMAGESELECTOR_H
#define VRN_IMAGESELECTOR_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"

namespace voreen {

/**
 * Selects a single volume out of an input image sequence.
 */
class VRN_CORE_API ImageSelector : public RenderProcessor {

public:
    ImageSelector();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "ImageSelector";   }
    virtual std::string getCategory() const   { return "Input";           }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE; }
    virtual bool isUtility() const            { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Selects a single image from the input image sequence and renders it to its outport.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void adjustPropertiesToInput();

    ImageSequencePort inport_;   ///< Inport for the image sequence.
    RenderPort outport_;         ///< The render port the selected image is written to.

    IntProperty imageID_;        ///< id of the selected image
    IntVec2Property imageSize_;  ///< Read-only property providing the dimensions of the output image.
    MWheelNumPropInteractionHandler<int> wheelHandler_; ///< used for cycling through images by mouse wheel

    tgt::Shader* shader_;        ///< Shader program used for rendering the selected image to the outport.

    static const std::string loggerCat_;
};

} // namespace

#endif
