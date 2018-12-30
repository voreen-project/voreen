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

#ifndef VRN_COLORDEPTH_H
#define VRN_COLORDEPTH_H

#include "voreen/core/processors/imageprocessorbypassable.h"
#include "voreen/core/voreenapplication.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"

namespace voreen {

/**
 * Performs a color filtering which encodes depth information.
 *
 */
class VRN_CORE_API ColorDepth : public ImageProcessorBypassable {
public:
    ColorDepth();
    virtual ~ColorDepth();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "ColorDepth";       }
    virtual std::string getCategory() const  { return "Image Processing"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;  }

protected:
    virtual void setDescriptions() {
        setDescription("Performs a color adaptation based on the depth information. There are four modes available. 'Light-dark' computes a gray value depending on the depth value. The current color can be modulated or replaced with this gray value. Additionally, a chromadepth or pseudo chromadepth color coding can be applied [Ropinski et al., Smartgraphics 2006].");
    }

    virtual void initialize();

    void process();

    tgt::Texture* chromaDepthTex_; ///< The texture used for the chromadepth color coding
    IntOptionProperty colorMode_; ///< The color mode to choose (Light-dark (replace/modulate), chromadepth, pseudo chromadepth)

    FloatProperty factor_;  ///< Controls the influence of the depth values

    RenderPort inport_;
    RenderPort outport_;

    static const std::string loggerCat_;
};


} // namespace voreen

#endif //VRN_COLORDEPTHPP_H
