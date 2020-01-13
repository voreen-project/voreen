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

#ifndef VRN_SAMPLERENDERPROCESSOR_H
#define VRN_SAMPLERENDERPROCESSOR_H

//header of base class
#include "voreen/core/processors/renderprocessor.h"
//port headers
#include "voreen/core/ports/renderport.h"
//property headers
#include "voreen/core/properties/floatproperty.h"
//header of the used shader object
#include "tgt/shadermanager.h"

//use namespace voreen
namespace voreen {

/**
 * Sample render processor, which gray-scales an input image based on a user-defined parameter.
 * VRN_CORE_API is a macro needed for shared libs on windows (see voreencoreapi.h)
 */
class VRN_CORE_API SampleRenderProcessor : public RenderProcessor {

public:
    /**
     * Constructor
     */
    SampleRenderProcessor();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new SampleRenderProcessor();     }
    virtual std::string getClassName() const { return "SampleRenderProcessor";  }
    virtual std::string getCategory() const  { return "Image Processing";       }

protected:

    virtual void setDescriptions() { setDescription("Sample render processor that transforms its input image to gray scale."); }
    virtual void process();

    /**
     * Overwrites the base implementation of this function.
     * It is used to load the needed shader.
     * @see Processor
     */
    virtual void initialize();

    /**
     * Overwrites the base implementation of this function.
     * It is used to free the used shader.
     * @see Processor
     */
    virtual void deinitialize();

private:

    //-------------
    //  members
    //-------------
    RenderPort inport_;             ///< input of the image which should be modified
    RenderPort outport_;            ///< output of the modified image
    FloatProperty saturationProp_;  ///< property for the color saturation parameter

    tgt::Shader* shader_;           ///< GLSL shader object used in process()
};

} // namespace

#endif // VRN_SAMPLERENDERPROCESSOR_H
